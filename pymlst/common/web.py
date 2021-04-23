import os
import tempfile

import zipfile
import requests
import questionary

from bs4 import BeautifulSoup
from Bio import SeqIO


class StructureError(Exception):
    pass


def request(query):
    result = requests.get(query, timeout=300)
    result.raise_for_status()
    return result


def display_prompt(message, choices):
    style = questionary.Style([
        ('qmark', 'fg:#673ab7 bold'),
        ('question', 'bold'),
        ('answer', 'fg:#f44336 bold'),
        ('pointer', 'fg:#673ab7 bold'),
        ('highlighted', 'fg:#673ab7 bold'),
        ('selected', 'fg:#cc5454'),
        ('separator', 'fg:#cc5454'),
        ('instruction', ''),
        ('text', ''),
    ])

    return questionary.select(message,
                              choices,
                              style=style) \
                      .ask()


def is_mlst_scheme(url, description):
    desc_lower = description.lower()
    blacklist = ['cgmlst', 'wgmlst', 'extended mlst']
    for word in blacklist:
        if word in desc_lower:
            return False
    scheme_json = request(url).json()
    if 'profiles_csv' not in scheme_json:
        return False
    return len(scheme_json['loci']) < 10


def process_results(choices, query, prompt):
    choices_length = len(choices)
    if choices_length == 0:
        raise Exception('No result found for \'{}\'\n'.format(query))
    if choices_length == 1:
        return choices[0]
    if prompt:
        return display_prompt('({}) Results found'.format(choices_length),
                              choices)
    raise Exception('More than 1 result found for \'{}\'\n'.format(query))


def retrieve_mlst(query, prompt_enabled, mlst=''):
    url = 'https://rest.pubmlst.org/db'

    try:
        whole_base = request(url).json()
    except ValueError:
        raise StructureError()

    species = {}
    query_low = query.lower()

    try:
        for record in whole_base:
            if record['name'] == 'test':
                continue
            for database in record['databases']:
                des = database['description'].replace('sequence/profile definitions', '').lower()
                if query_low in des:
                    if database['name'].endswith('seqdef'):
                        species[des] = database['href']
    except KeyError:
        raise StructureError()

    species_choice = process_results(list(species.keys()), query, prompt_enabled)

    if species_choice is None:
        return None

    schemes_url = species[species_choice] + '/schemes'
    schemes_json = request(schemes_url).json()

    schemes = {}
    mlst_low = mlst.lower()

    try:
        for scheme in schemes_json['schemes']:
            if not is_mlst_scheme(scheme['scheme'], scheme['description']):
                continue
            des = scheme['description'].lower()
            if mlst_low in des:
                schemes[des] = scheme['scheme']
    except KeyError:
        raise StructureError()

    scheme_choice = process_results(list(schemes.keys()), mlst, prompt_enabled)

    return schemes[scheme_choice]


def retrieve_cgmlst(query, prompt_enabled):
    url = 'https://www.cgmlst.org/ncs'
    page = request(url)

    soup = BeautifulSoup(page.content, 'html.parser')

    table = soup.find('tbody')
    if table is None:
        raise StructureError()

    lines = table.find_all('a')

    addresses = {}
    query_low = query.lower()

    for line in lines:
        text = line.get_text()
        if 'cgMLST' not in text:
            continue
        name = text.replace('cgMLST', '').strip()
        if query_low in name.lower():
            link = line.get('href')
            if link is None:
                raise StructureError()
            addresses[name] = link

    choice = process_results(list(addresses.keys()), query, prompt_enabled)
    genome_url = addresses[choice]

    return genome_url + 'alleles'


def get_coregene_file(url, handle):
    with tempfile.TemporaryDirectory() as tmp_dir:
        zip_req = request(url)
        zip_tmp = os.path.join(tmp_dir, 'tmp.zip')
        open(zip_tmp, 'wb').write(zip_req.content)

        fas_tmp = os.path.join(tmp_dir, 'fas')
        os.mkdir(fas_tmp)
        with zipfile.ZipFile(zip_tmp) as z_file:
            z_file.extractall(fas_tmp)
        skipped = []
        for fasta in os.listdir(fas_tmp):
            try:
                iterator = next(SeqIO.parse(os.path.join(fas_tmp, fasta), 'fasta'))
            except (StopIteration, ValueError, TypeError):
                skipped.append(fasta)
                continue
            handle.write('> ' + fasta.replace('.fasta', '') + '\n')
            handle.write(str(iterator.seq) + '\n')
        return skipped


def clean_csv(csv_content, locus_nb):
    lines = csv_content.split('\n')
    header = lines[0].split('\t')
    diff = len(header) - (locus_nb + 1)
    if diff > 0:
        lines[0] = '\t'.join(header[0:-diff])
    return '\n'.join(lines)


def get_mlst_files(directory, url):
    mlst_scheme = request(url).json()

    # Downloading the locus files in a directory :
    locus_dir = os.path.join(directory, 'locus')
    os.mkdir(locus_dir)
    for loci in mlst_scheme['loci']:
        name = loci.split('/')[-1]
        loci_fasta = request(loci + '/alleles_fasta')
        loci_file_name = os.path.join(locus_dir, name + '.fasta')
        open(loci_file_name, 'wb').write(loci_fasta.content)

    # Downloading the profiles CSV + removing last header column :
    profiles_url = url + '/profiles_csv'
    profiles = request(profiles_url)
    with open(os.path.join(directory, 'profiles.csv'), 'wt') as profiles_dir:
        profiles_dir.write(clean_csv(profiles.text, len(mlst_scheme['loci'])))
