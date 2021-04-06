from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
import logging


def read_genome(genome):
    seqs = {}
    for seq in SeqIO.parse(genome, 'fasta'):
        seqs[seq.id] = seq
    return seqs


def strip_file(file):
    found = []
    if file is not None:
        for line in file.readLines():
            found.append(line.rstrip('\n'))
    return []


def compar_seqs(seqs):
    count = 0
    dim = len(seqs[0])
    for j in range(0, len(seqs[0])):
        d = set([s[j] for s in seqs])
        if '-' in d:
            d.remove('-')
        if len(d) > 1:
            count += 1
    return count


def write_count(count, texte):
    if count:
        count.write(texte)


def validate_sequence(sequence):
    try:
        sequence.translate(cds=True, table=11)
    except TranslationError:
        return False
    else:
        return True


def create_logger():
    logging.basicConfig(
        level=logging.INFO,
        format='[%(levelname)s: %(asctime)s] %(message)s')
    return logging.getLogger('mlst')