import time
import os

import pymlst
from pymlst.wg.extractors import TableExtractor, SequenceExtractor

genome_path = '/home/abordy/workspace/data/FATgenome/genome'
db_path = '/home/abordy/workspace/data/FATgenome/database.db'
script = '/home/abordy/workspace/pyMLST/mlst_add_strain.py'


low = 20
up = 40


def with_reopen():
    for file_name in os.listdir(genome_path)[low:up]:
        file_path = os.path.join(genome_path, file_name)
        with pymlst.open_wg(db_path) as mlst:
            mlst.add_strain(open(file_path))


def with_keepopen():
    i = 0
    with pymlst.open_wg(db_path) as mlst:
        for file_name in os.listdir(genome_path)[low:up]:
            i += 1
            if i % 50 == 0:
                mlst.commit()
            file_path = os.path.join(genome_path, file_name)
            mlst.add_strain(open(file_path))


def with_old():
    for file_name in os.listdir(genome_path)[low:up]:
        file_path = os.path.join(genome_path, file_name)
        os.system('python {} {} {}'.format(script, file_path, db_path))


def calculate_resemblance(seq1, seq2):
    alike = 0
    for s1, s2 in zip(seq1, seq2):
        if s1 == s2:
            alike += 1
    return alike / len(seq1)


if __name__ == '__main__':
    start = time.time()

    with_keepopen()
    # with pymlst.open_wg('/home/abordy/workspace/data/database_BIG.db') as mlst:
    #     mlst.extract(TableExtractor(export='stat'))

    elapsed = time.time() - start
    print('Elapsed: {}s'.format(elapsed))
