import os
import tempfile

import zipfile
import requests
import questionary

from bs4 import BeautifulSoup
from Bio import SeqIO


class StructureException(Exception):
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


def prompt_mlst(query, prompt_enabled, mlst=None):
    url = 'https://rest.pubmlst.org/db'
    whole_base = request(url).json()

    species_choices = []
    species = {}
    query_low = query.lower()

    for record in whole_base:
        if record['name'] == 'test':
            continue
        for database in record['databases']:
            des = database['description'].replace('sequence/profile definitions', '').strip().lower()
            if query_low in des:
                if database['name'].endswith('seqdef'):
                    species_choices.append({'name': des})
                    species[des] = database['href']

    species_length = len(species_choices)
    if species_length == 0:
        raise Exception('No match found for \'' + query + '\'\n')
    if species_length == 1:
        species_choice = species_choices[0]['name']
    elif not prompt_enabled:
        raise Exception('More than 1 match found for \'' + query + '\'\n')
    else:
        species_choice = display_prompt('(' + str(species_length) + ') Matching species found',
                                        species_choices)
        if species_choice is None:
            return None

    schemes_url = species[species_choice] + '/schemes'
    schemes_json = request(schemes_url).json()

    schemes_choices = []
    schemes = {}
    matching_mlst = None

    for scheme in schemes_json['schemes']:
        if is_mlst_scheme(scheme['scheme'], scheme['description']):
            if mlst is None:
                schemes_choices.append({'name': scheme['description']})
                schemes[scheme['description']] = scheme['scheme']
            else:
                if mlst.lower() in scheme['description'].lower():
                    if matching_mlst is None:
                        matching_mlst = scheme['scheme']
                    else:
                        raise Exception('More than 1 MLST scheme found matching \'' + mlst + '\'\n')

    if mlst is not None:
        if matching_mlst is None:
            raise Exception('No MLST scheme found matching \'' + mlst + '\'')
        return matching_mlst

    schemes_length = len(schemes)
    if schemes_length == 1:
        return schemes[schemes_choices[0]['name']]
    elif schemes_length == 0:
        raise Exception('No MLST scheme found')
    elif not prompt_enabled:
        raise Exception('More than 1 MLST scheme found')

    scheme_choice = display_prompt('Choose an MLST scheme', schemes_choices)
    if scheme_choice is None:
        return None

    return schemes[scheme_choice]


def prompt_cgmlst(query, prompt_enabled):
    url = 'https://www.cgmlst.org/ncs'
    page = request(url)

    soup = BeautifulSoup(page.content, 'html.parser')

    table = soup.find('tbody')
    if table is None:
        raise StructureException()

    lines = table.find_all('a')

    choices = []
    addresses = {}
    query_low = query.lower()

    for line in lines:
        text = line.get_text()
        if 'cgMLST' not in text:
            continue
        name = text.replace('cgMLST', '').strip()
        if query_low in name.lower():
            choices.append({'name': name})
            link = line.get('href')
            if link is None:
                raise StructureException()
            addresses[name] = link

    choices_length = len(choices)
    if choices_length == 1:
        genome_url = addresses[choices[0]['name']]
    elif choices_length == 0:
        raise Exception('No coregenome found')
    elif not prompt_enabled:
        raise Exception('More than 1 coregenome found')
    else:
        choice = display_prompt('Select a coregenome', choices)
        if choice is None:
            return ''
        genome_url = addresses[choice]

    return genome_url + 'alleles'


def build_coregene(url, handle):
    with tempfile.TemporaryDirectory() as tmp_dir:
        zip_req = request(url)
        zip_tmp = tmp_dir + '/tmp.zip'
        open(zip_tmp, 'wb').write(zip_req.content)

        fas_tmp = tmp_dir + '/fas'
        os.mkdir(fas_tmp)
        with zipfile.ZipFile(zip_tmp) as z_file:
            z_file.extractall(fas_tmp)
        skipped = []
        for fasta in os.listdir(fas_tmp):
            try:
                iterator = next(SeqIO.parse(fas_tmp + '/' + fasta, 'fasta'))
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
    locus_dir = directory + '/locus'
    os.mkdir(locus_dir)
    for loci in mlst_scheme['loci']:
        name = loci.split('/')[-1]
        loci_fasta = request(loci + '/alleles_fasta')
        loci_file_name = locus_dir + '/' + name + '.fasta'
        open(loci_file_name, 'wb').write(loci_fasta.content)

    # Downloading the profiles CSV + removing last header column :
    profiles_url = url + '/profiles_csv'
    profiles = request(profiles_url)
    with open(directory + '/profiles.csv', 'wt') as profiles_dir:
        profiles_dir.write(clean_csv(profiles.text, len(mlst_scheme['loci'])))
