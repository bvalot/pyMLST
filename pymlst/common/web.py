import os
import tempfile
import logging

import zipfile
import requests
import questionary

from bs4 import BeautifulSoup
from Bio import SeqIO

import urllib3
urllib3.disable_warnings()

from pymlst.common.exceptions import PyMLSTError

PUBMLST_URL = 'https://rest.pubmlst.org/db'
CGMLST_URL = 'https://www.cgmlst.org/ncs'


class PyMLSTWebError(PyMLSTError):
    pass


class StructureError(PyMLSTWebError):
    pass


def request(query):
    result = requests.get(query, timeout=600, verify=False)
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
        raise PyMLSTWebError('No result found for \'{}\'\n'.format(query))
    if choices_length == 1:
        logging.info("One element found : {}".format(choices[0]))
        return choices[0]
    if prompt:
        logging.info("{} elements found, please choose one:".format(str(len(choices))))       
        return display_prompt('({}) Results found'.format(choices_length),
                              choices)
    raise PyMLSTWebError('More than 1 result found for \'{}\'\n'.format(query))


def get_mlst_species(query):
    """Gets MLST species from pubmlst.org.

    :param query: A sub-string to filter species names.
    :return: A Dictionary with species name in Key and URL in Value.
    """
    try:
        whole_base = request(PUBMLST_URL).json()
    except ValueError as error:
        raise StructureError() from error

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
    except KeyError as error:
        raise StructureError() from error

    return species


def get_mlst_schemes(species_url, query):
    """Gets schemes profiles from PubMLST for a given species URL.

    :param species_url: The species URL (see get_mlst_species()).
    :param query: A sub-string to filter schemes names.
    :return: A Dictionary with schemes name in Key and URL in Value.
    """
    schemes_url = species_url + '/schemes'
    schemes_json = request(schemes_url).json()

    schemes = {}
    query_low = query.lower()

    try:
        for scheme in schemes_json['schemes']:
            if not is_mlst_scheme(scheme['scheme'], scheme['description']):
                continue
            des = scheme['description'].lower()
            if query_low in des:
                schemes[des] = scheme['scheme']
    except KeyError as error:
        raise StructureError() from error

    return schemes


def retrieve_mlst(query, prompt_enabled, mlst=''):
    """Retrieves MLST data, prompts user if necessary and if possible.

    :param query: A sub-string to filter species names.
    :param prompt_enabled: Whether or not to prompt user for actions.
                           If disabled and many choices are possible,
                           will raise an Exception.
    :param mlst: A sub-string to filter schemes names.
    :return: A scheme URL.
    """
    species = get_mlst_species(query)
    species_choice = process_results(list(species.keys()), query, prompt_enabled)
    if species_choice is None:
        return None

    species_url = species[species_choice]

    schemes = get_mlst_schemes(species_url, mlst)
    scheme_choice = process_results(list(schemes.keys()), mlst, prompt_enabled)
    if scheme_choice is None:
        return None

    return schemes[scheme_choice]


def get_cgmlst_species(query):
    """Gets cgMLST species from cgmlst.org.

    :param query: A sub-string to filter species names.
    :return: A Dictionary with species name in Key and download URL in Value.
    """
    page = request(CGMLST_URL)

    soup = BeautifulSoup(page.content, 'html.parser')

    table = soup.find('tbody')
    if table is None:
        raise StructureError()

    lines = table.find_all('a')

    species = {}
    query_low = query.lower()

    for line in lines:
        text = line.get_text()
        if 'cgMLST' not in text:
            continue
        name = text.replace('cgMLST', '').strip()
        if query_low in name.lower():
            url = line.get('href')
            if url is None:
                raise StructureError()
            species[name] = url

    return species


def retrieve_cgmlst(query, prompt_enabled):
    """Retrieves cgMLST data, prompts user if necessary and if possible.

    :param query: A sub-string to filter species names.
    :param prompt_enabled: Whether or not to prompt user for actions.
                           If disabled and many choices are possible,
                           will raise an Exception.
    :return: A species download URL.
    """
    species = get_cgmlst_species(query)
    choice = process_results(list(species.keys()), query, prompt_enabled)
    if choice is None:
        return None

    species_url = species[choice]

    return species_url


def get_cgmlst_file(url, handle):
    """Download cgMLST data and use them to initialize a fasta file.

    :param url: The download URL.
    :param handle: The file handle.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        url += 'alleles'
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


def get_mlst_files(url, directory):
    """Download MLST data and puts them in the given directory.

    :param url: The scheme URL.
    :param directory: The directory.
    """
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
