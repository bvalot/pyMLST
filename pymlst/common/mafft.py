import tempfile
import logging
import sys
import subprocess

from Bio import AlignIO

from pymlst import config
from pymlst.common import utils, exceptions


def align(genes):
    path = config.get_binary_path('mafft')
    if not path:
        raise exceptions.BinaryNotFound('MAFFT binary was not found')
    with tempfile.NamedTemporaryFile(mode='w+t', suffix='.fasta') as tmp:
        utils.write_genome(genes, tmp)
        tmp.flush()
        p = subprocess.Popen([path, "--auto", tmp.name], \
                             stdout=subprocess.PIPE, \
                             stderr=subprocess.PIPE, \
                             encoding=sys.stdout.encoding)
        records = AlignIO.parse(p.stdout, "fasta")
        try:
            alignments = next(records)
        except StopIteration:
            logging.error("MAFFT doesn't finish correctly\n" + \
                          p.stderr.read())
            return {}
        return utils.records_to_dict(alignments)


def __first_aligned_position(sequence):
    position = 0
    for char in sequence:
        if char != '-':
            return position
        position += 1
    return -1


def get_aligned_area(query, target):
    alignments = align({'query': query, 'target': target})
    if len(alignments) != 2:
        return None, None
    q_align = alignments['query']
    q_len = len(q_align)
    start_index = __first_aligned_position(q_align)
    if start_index == -1:
        return None, None
    end_index = q_len - __first_aligned_position(reversed(q_align))
    return start_index, end_index
