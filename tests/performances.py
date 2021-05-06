import time
import os

import pymlst


genome_path = '/home/abordy/workspace/data/FATgenome/genome'
db_path = '/home/abordy/workspace/data/FATgenome/database.db'
script = '/home/abordy/workspace/pyMLST/mlst_add_strain.py'


low = 20
up = 40


def with_reopen():
    for file_name in os.listdir(genome_path):
        file_path = os.path.join(genome_path, file_name)
        with pymlst.open_wg(db_path) as mlst:
            mlst.add_strain(open(file_path))


def with_keepopen():
    i = 0
    with pymlst.open_wg(db_path) as mlst:
        for file_name in os.listdir(genome_path):
            i += 1
            if i % 50 == 0:
                mlst.commit()
            file_path = os.path.join(genome_path, file_name)
            mlst.add_strain(open(file_path))


def with_old():
    for file_name in os.listdir(genome_path):
        file_path = os.path.join(genome_path, file_name)
        os.system('python {} {} {}'.format(script, file_path, db_path))


def calculate_resemblance(seq1, seq2):
    alike = 0
    for s1, s2 in zip(seq1, seq2):
        if s1 == s2:
            alike += 1
    return alike / len(seq1)


def process_similarities(gene):
    with pymlst.open_wg(db_path) as mlst:
        sequences = mlst.database.get_sequences_by_gene(gene)
        ref_seq = sequences[0][1]
        different = set()
        for seq in sequences:
            resemblance = calculate_resemblance(ref_seq, seq[1])
            print('{} %'.format(resemblance * 100))
            print('Seq: {}'.format(seq[1]))
            different.add(resemblance)
        print('There is {} different sequences stored'.format(len(different)))


if __name__ == '__main__':
    start = time.time()

    process_similarities('PA4819')

    elapsed = time.time() - start
    print('Elapsed: {}s'.format(elapsed))
