import os
import pytest
from sqlalchemy import select, exists
from sqlalchemy.sql.functions import count
from sqlalchemy.sql.operators import in_op as in_

import pymlst
from pymlst.wg.core import DatabaseWG


data_path = os.path.join(os.path.dirname(__file__), 'data')
wg_path = os.path.join(data_path, 'wg')


def fasta(name):
    return open(os.path.join(wg_path, name + '.fasta'))


@pytest.fixture()
def wg():
    with pymlst.open_wg() as wg_mlst:
        yield wg_mlst


@pytest.fixture()
def db():
    db = DatabaseWG(None, 'ref')
    yield db
    db.close()


@pytest.fixture()
def db_simple(db):
    seqs = [
        ('g1', 'AAA'),
        ('g1', 'ATA'),
        ('g2', 'TTT'),
        ('g3', 'CCC'),
    ]
    for gene, seq in seqs:
        _, seq_id = db.add_sequence(seq)
        db.add_mlst('A', gene, seq_id)
    db.add_mlst('A', 'g4', 4)  # g3 and g4 both on CCC
    return db


@pytest.fixture()
def db_many(db):
    seqs = [
        ('A', 'g1', 'AAA'),
        ('A', 'g2', 'ATA'),
        ('A', 'g3', 'ATT'),
        ('A', 'g4', 'CCC'),
        ('B', 'g1', 'AAT'),
        ('B', 'g2', 'ATA'),
        ('B', 'g3', 'TTT'),
        ('B', 'g4', 'CAC'),
        ('B', 'g5', 'GGG'),
        ('C', 'g1', 'AAA'),
        ('C', 'g3', 'TTT'),
        ('C', 'g4', 'CAA'),
        ('D', 'g4', 'CAC'),
    ]
    for strain, gene, seq in seqs:
        _, seq_id = db.add_sequence(seq)
        db.add_mlst(strain, gene, seq_id)
    return db


@pytest.fixture()
def db_ref(db_many):
    seqs = [
        ('g1', 'AAA'),
        ('g2', 'ATA'),
        ('g3', 'TTT'),
        ('g4', 'CCC'),
        ('g5', 'GGG'),
    ]
    for gene, seq in seqs:
        _, seq_id = db_many.add_sequence(seq)
        db_many.add_mlst('ref', gene, seq_id)
    return db_many


def test_add_sequence(db):
    inserted, key = db.add_sequence('AAA')
    assert inserted
    seq = db.connection.execute(
        select([db.sequences.c.sequence])
        .where(db.sequences.c.id == key)
    ).fetchone()
    assert seq.sequence == 'AAA'


def test_add_sequence_exist(db):
    db.add_sequence('AAA')
    inserted, key = db.add_sequence('AAA')
    assert not inserted
    seq = db.connection.execute(
        select([db.sequences.c.id])
            .where(db.sequences.c.sequence == 'AAA')
    ).fetchall()
    assert len(seq) == 1
    assert seq[0].id == key


def test_add_mlst(db):
    db.add_mlst('A', 'g1', 0)
    res = db.connection.execute(
        select([db.mlst])
    ).fetchall()
    assert len(res) == 1
    assert (res[0].souche == 'A'
            and res[0].gene == 'g1'
            and res[0].seqid == 0)


def test_concatenate_gene(db):
    db.add_mlst('A', 'g1', 0)
    db.concatenate_gene(0, 'g2')
    res = db.connection.execute(
        select([db.mlst.c.gene])
    ).fetchall()
    assert len(res) == 1
    assert res[0].gene == 'g1;g2'


def test_get_sequences_by_gene(db_simple):
    seqs = db_simple.get_sequences_by_gene('g1')
    assert len(seqs) == 2
    assert (seqs[0].sequence == 'AAA'
            and seqs[1].sequence == 'ATA')


def test_remove_sequences(db_simple):
    db_simple.remove_sequences([1, 3])
    mlst_c = db_simple.connection.execute(
        select([count(db_simple.mlst.c.id)])
    ).scalar()
    assert mlst_c == 3
    seq_c = db_simple.connection.execute(
        select([count(db_simple.sequences.c.id)])
    ).scalar()
    assert seq_c == 2
    seq_e = db_simple.connection.execute(
        select([db_simple.sequences])
        .where(in_(db_simple.sequences.c.id, [1, 3]))
    ).fetchone()
    assert seq_e is None
    mlst_e = db_simple.connection.execute(
        select([db_simple.mlst])
        .where(in_(db_simple.mlst.c.seqid, [1, 3]))
    ).fetchone()
    assert mlst_e is None


def test_remove_gene(db_simple):
    db_simple.remove_gene('g1')
    mlst_e = db_simple.connection.execute(
        select([db_simple.mlst])
        .where(db_simple.mlst.c.gene == 'g1')
    ).fetchone()
    assert mlst_e is None


def test_remove_orphan_sequences(db_simple):
    db_simple.remove_gene('g2')
    db_simple.remove_gene('g4')
    db_simple.remove_orphan_sequences([3, 4])
    seq_1_e = db_simple.connection.execute(
        select([db_simple.sequences])
        .where(db_simple.sequences.c.id == 3)
    ).fetchone()
    assert seq_1_e is None
    seq_2_e = db_simple.connection.execute(
        select([db_simple.sequences])
        .where(db_simple.sequences.c.id == 4)
    ).fetchone()
    assert seq_2_e is not None


def test_remove_strain(db_simple):
    db_simple.remove_strain('A')
    mlst_e = db_simple.connection.execute(
        select([db_simple.mlst])
        .where(db_simple.mlst.c.souche == 'A')
    ).fetchone()
    assert mlst_e is None


def test_get_gene_sequences_ids(db_simple):
    ids = db_simple.get_gene_sequences_ids('g1')
    assert ids == {1, 2}


def test_get_strain_sequences_ids(db_many):
    ids = db_many.get_strain_sequences_ids('A')
    assert ids == {1, 2, 3, 4}


def test_contains_souche(db_many):
    assert db_many.contains_souche('B')
    db_many.connection.execute(
        db_many.mlst.delete()
        .where(db_many.mlst.c.souche == 'B'))
    assert not db_many.contains_souche('B')


def test_get_gene_sequences_many_strains(db_ref):
    g1_seq = db_ref.get_gene_sequences('g1')
    assert g1_seq == [
        [1, ['A', 'C'], 'AAA'],
        [5, ['B'], 'AAT']]
    g2_seq = db_ref.get_gene_sequences('g2')
    assert g2_seq == [
        [2, ['A', 'B'], 'ATA']]


def test_get_gene_sequences_one_strain_duplicated_gene(db_simple):
    g1_seq = db_simple.get_gene_sequences('g1')
    assert g1_seq == [
        [1, ['A'], 'AAA'],
        [2, ['A'], 'ATA']]


def test_get_duplicated_genes(db_simple):
    dupli = db_simple.get_duplicated_genes()
    assert dupli == {'g1'}


def test_get_all_strains(db_ref):
    strains = db_ref.get_all_strains()
    assert strains == ['A', 'B', 'C', 'D']


def test_get_core_genes(db_ref):
    genes = db_ref.get_core_genes()
    assert genes == ['g1', 'g2', 'g3', 'g4', 'g5']


def test_count_sequences_per_gene(db_ref):
    seq_c = db_ref.count_sequences_per_gene()
    assert seq_c == {
        'g1': 2,
        'g2': 1,
        'g3': 2,
        'g4': 3,
        'g5': 1
    }


def test_count_souches_per_gene(db_ref):
    str_c = db_ref.count_souches_per_gene()
    assert str_c == {
        'g1': 3,
        'g2': 2,
        'g3': 3,
        'g4': 4,
        'g5': 1
    }


def test_count_genes_per_souche(db_ref):
    gene_c = db_ref.count_genes_per_souche(['g1', 'g2', 'g3', 'g4', 'g5'])
    assert gene_c == {
        'A': 4,
        'B': 5,
        'C': 3,
        'D': 1,
        'ref': 5
    }


def test_count_sequences(db_ref):
    seq_c = db_ref.count_sequences()
    assert seq_c == 9


def test_get_strains_distances(db_ref):
    distances = db_ref.get_strains_distances(['g1', 'g2', 'g3', 'g4', 'g5'])
    assert distances == {
        'A': {
            'B': 3,
            'C': 2,
            'D': 1,
        },
        'B': {
            'A': 3,
            'C': 2,
        },
        'C': {
            'A': 2,
            'B': 2,
            'D': 1,
        },
        'D': {
            'A': 1,
            'C': 1,
        },
    }


def test_get_mlst(db_ref):
    mlst = db_ref.get_mlst(['g1', 'g2', 'g3', 'g4', 'g5'])
    assert mlst == {
        'g1': {
            'A': '1',
            'B': '5',
            'C': '1',
        },
        'g2': {
            'A': '2',
            'B': '2',
        },
        'g3': {
            'A': '3',
            'B': '6',
            'C': '6',
        },
        'g4': {
            'A': '4',
            'B': '7',
            'C': '9',
            'D': '7',
        },
        'g5': {
            'B': '8',
        },
    }
