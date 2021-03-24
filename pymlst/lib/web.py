import requests
import zipfile
import tempfile
import os

from bs4 import BeautifulSoup
from Bio import SeqIO

from PyInquirer import style_from_dict, Token, prompt


def display_prompt(name, message, choices):
    style = style_from_dict({
        Token.Separator: '#cc5454',
        Token.QuestionMark: '#673ab7 bold',
        Token.Selected: '#cc5454',  # default
        Token.Pointer: '#673ab7 bold',
        Token.Instruction: '',  # default
        Token.Answer: '#f44336 bold',
        Token.Question: '',
    })

    questions = [
        {
            'type': 'list',
            'message': message,
            'name': name,
            'choices': choices,
            'validate': lambda answer: 'You must choose an option' \
                if len(answer) == 0 else True
        }
    ]

    try:
        return prompt(questions, style=style)[name]
    except KeyError:
        return None


def prompt_mlst(query, prompt_enabled, mlst=None):
    url = 'https://rest.pubmlst.org/db'
    whole_base = requests.get(url).json()

    species_choices = []
    species = {}

    for record in whole_base:
        if record['name'] == 'test':
            continue
        for db in record['databases']:
            des = db['description'].replace('sequence/profile definitions', '').strip().lower()
            if query.lower() in des:
                if db['name'].endswith('seqdef'):
                    species_choices.append({'name': des})
                    species[des] = db['href']

    species_length = len(species_choices)
    if species_length == 0:
        raise Exception('No match found for \'' + query + '\'\n')
    elif species_length == 1:
        species_choice = species_choices[0]['name']
    elif not prompt_enabled:
        raise Exception('More than 1 match found for \'' + query + '\'\n')
    else:
        species_choice = display_prompt('species',
                                     '(' + str(species_length) + ') Matching species found',
                                     species_choices)
        if species_choice is None:
            return None

    schemes_url = species[species_choice] + '/schemes'
    schemes_json = requests.get(schemes_url).json()

    schemes_choices = []
    schemes = {}

    for s in schemes_json['schemes']:
        if mlst is not None:
            if mlst.lower() in s['description'].lower():
                return s['scheme']

        if s['description'].startswith('MLST'):
            schemes_choices.append({'name': s['description']})
            schemes[s['description']] = s['scheme']

    if mlst is not None:
        raise Exception('No MLST schema was found matching \'' + mlst + '\'')

    schemes_length = len(schemes)
    if schemes_length == 1:
        return schemes[schemes_choices[0]['name']]
    elif schemes_length == 0:
        raise Exception('No MLST scheme was found')
    elif not prompt_enabled:
        raise Exception('Multiple MLST schemes were found')

    scheme_choice = display_prompt('scheme', 'Choose an MLST scheme', schemes_choices)
    if scheme_choice is None:
        return None

    return schemes[scheme_choice]


def prompt_cgmlst():
    url = 'https://www.cgmlst.org/ncs'
    page = requests.get(url)

    soup = BeautifulSoup(page.content, 'html.parser')

    table = soup.find('tbody')
    lines = table.find_all('a')

    choices = []
    addresses = {}

    for l in lines:
        name = l.get_text().replace('cgMLST', '').strip()
        choices.append({'name': name})
        addresses[name] = l.get('href')

    style = style_from_dict({
        Token.Separator: '#cc5454',
        Token.QuestionMark: '#673ab7 bold',
        Token.Selected: '#cc5454',  # default
        Token.Pointer: '#673ab7 bold',
        Token.Instruction: '',  # default
        Token.Answer: '#f44336 bold',
        Token.Question: '',
    })

    questions = [
        {
            'type': 'list',
            'message': 'Select a coregenome',
            'name': 'coregenome',
            'choices': choices,
            'validate': lambda answer: 'You must choose a coregenome.' \
                if len(answer) == 0 else True
        }
    ]

    answer = prompt(questions, style=style)

    try:
        dll_url = addresses.get(answer['coregenome']) + 'alleles'
    except KeyError:
        return ''

    return dll_url


def build_coregene(url, handle):
    with tempfile.TemporaryDirectory() as tmp_dir:
        zip_req = requests.get(url)
        zip_tmp = tmp_dir + '/tmp.zip'
        open(zip_tmp, 'wb').write(zip_req.content)

        fas_tmp = tmp_dir + '/fas'
        os.mkdir(fas_tmp)
        with zipfile.ZipFile(zip_tmp) as z:
            z.extractall(fas_tmp)

        for fasta in os.listdir(fas_tmp):
            try:
                it = next(SeqIO.parse(fas_tmp + '/' + fasta, 'fasta'))
            except StopIteration:
                continue
            handle.write('> ' + fasta.replace('.fasta', '') + '\n')
            handle.write(str(it.seq) + '\n')


def clean_csv(csv_content, locus_nb):
    lines = csv_content.split('\n')
    header = lines[0].split('\t')
    diff = len(header) - (locus_nb + 1)
    if diff > 0:
        lines[0] = '\t'.join(header[0:-diff])
        print('cleaned: ', header[-diff])
    return '\n'.join(lines)


def get_mlst_files(directory, url):
    mlst_scheme = requests.get(url).json()

    # Downloading the locus files in a directory :
    locus_dir = directory + '/locus'
    os.mkdir(locus_dir)
    for loci in mlst_scheme['loci']:
        name = loci.split('/')[-1]
        loci_fasta = requests.get(loci + '/alleles_fasta')
        loci_file_name = locus_dir + '/' + name + '.fasta'
        open(loci_file_name, 'wb').write(loci_fasta.content)

    # Downloading the profiles CSV + removing last header column :
    profiles_url = url + '/profiles_csv'
    profiles = requests.get(profiles_url)
    with open(directory + '/profiles.csv', 'wt') as p:
        p.write(clean_csv(profiles.text, len(mlst_scheme['loci'])))
