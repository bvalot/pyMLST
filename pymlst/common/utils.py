from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
import logging


def records_to_dict(records):
    print('records: ' + str(type(records)))
    seq_dict = {}
    for seq in records:
        print('seq: ' + str(type(seq)))
        seq_dict[seq.id] = seq
    return seq_dict


def read_genome(handle):
    records = SeqIO.parse(handle, 'fasta')
    return records_to_dict(records)


def write_genome(genome_dict, handle):
    for seq_id, seq in genome_dict.items():
        handle.write('> ' + seq_id + '\n'
                     + seq + '\n')


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