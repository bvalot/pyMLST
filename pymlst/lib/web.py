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

    return addresses.get(answer['coregenome']) + 'alleles'


def build_coregene(url, handle):
    with tempfile.TemporaryDirectory() as tmp_dir:
        zip_req = requests.get(url, allow_redirects=True)
        zip_tmp = tmp_dir + '/tmp.zip'
        open(zip_tmp, 'wb').write(zip_req.content)

        fas_tmp = tmp_dir + '/fas'
        os.mkdir(fas_tmp)
        with zipfile.ZipFile(zip_tmp) as z:
            z.extractall(fas_tmp)

        for fasta in os.listdir(fas_tmp):
            it = next(SeqIO.parse(fas_tmp + '/' + fasta, 'fasta'))
            handle.write('> ' + fasta.replace('.fasta', '') + '\n')
            handle.write(str(it.seq) + '\n')
