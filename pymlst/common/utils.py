import logging

from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError


def records_to_dict(records):
    seq_dict = {}
    for seq in records:
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
    for index in range(0, len(seqs[0])):
        seqs_char = set([s[index] for s in seqs])
        if '-' in seqs_char:
            seqs_char.remove('-')
        if len(seqs_char) > 1:
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


def create_logger(verbose=False):
    if verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(
        level=level,
        format='[%(levelname)s: %(asctime)s] %(message)s')


def clean_kwargs(kwargs):
    """Removes kwargs with None values produced by Click.

    Because of the way the Click library binds
    every arguments and options to kwargs entries,
    when a user doesn't specify an option, its name
    is bound to None in the kwargs dictionary.

    By removing the None entries we can pass the kwargs directly
    to the API core functions without overriding the default values.
    """
    for key, value in kwargs.copy().items():
        if value is None:
            kwargs.pop(key)
    return kwargs

