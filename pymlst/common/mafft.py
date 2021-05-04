import tempfile
from io import StringIO

from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline

from pymlst.common import binaries
from pymlst.common import utils


def align(genes):
    path = binaries.get_binary_path('mafft')
    if not path:
        raise Exception('Unable to locate the Mafft executable\n')
    with tempfile.NamedTemporaryFile(mode='w+t', suffix='.fasta') as tmp:
        utils.write_genome(genes, tmp)
        tmp.flush()
        mafft_cmd = MafftCommandline(path, input=tmp.name, quiet=True)
        stdout, _ = mafft_cmd()
        records = AlignIO.parse(StringIO(stdout), "fasta")
        try:
            alignments = next(records)
        except StopIteration:
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
