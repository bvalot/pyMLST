import requests
import zipfile
import tempfile
import os

from bs4 import BeautifulSoup
from Bio import SeqIO

from PyInquirer import style_from_dict, Token, prompt


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


def get_mlst_files(species, directory):
    base_url = 'https://rest.pubmlst.org/db/'
    url = base_url + '/pubmlst_' + species + '_seqdef/schemes/1'
    scheme = requests.get(url).json()
    locus_dir = directory + '/locus'
    os.mkdir(locus_dir)
    for loci in scheme['loci']:
        name = loci.split('/')[-1]
        loci_fasta = requests.get(loci + '/alleles_fasta')
        loci_file_name = locus_dir + '/' + name + '.fasta'
        open(loci_file_name, 'wb').write(loci_fasta.content)
    profiles_url = url + '/profiles_csv'
    profiles = requests.get(profiles_url)
    lines = profiles.text.split('\n')
    lines[0] = '\t'.join(lines[0].split('\t')[0:-1])
    with open(directory + '/profiles.csv', 'wt') as p:
        p.write('\n'.join(lines))
