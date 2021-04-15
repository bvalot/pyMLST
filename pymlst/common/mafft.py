import tempfile
from io import StringIO

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MafftCommandline

from pymlst.common import binaries
from pymlst.common.utils import records_to_dict, write_genome


def align(genes):
    path = binaries.get_binary_path('mafft')
    if not path:
        raise Exception('Unable to locate the Mafft executable\n')
    with tempfile.NamedTemporaryFile(mode='w+t', suffix='.fasta') as tmp:
        write_genome(genes, tmp)
        tmp.flush()
        mafft_cmd = MafftCommandline(path, input=tmp.name, quiet=True)
        stdout, stderr = mafft_cmd()
        records = AlignIO.parse(StringIO(stdout), "fasta")
        try:
            alignments = next(records)
        except StopIteration:
            return {}
        return records_to_dict(alignments)


def _first_aligned_position(sequence):
    position = 0
    for c in sequence:
        if c != '-':
            return position
        position += 1
    return -1


def get_aligned_area(query, target):
    alignments = align({'query': query, 'target': target})
    if not alignments:
        return None
    q_align = alignments['query']
    q_len = len(q_align)
    start_index = _first_aligned_position(q_align)
    if start_index == -1:
        return None
    end_index = q_len - _first_aligned_position(reversed(q_align))
    return start_index, end_index

