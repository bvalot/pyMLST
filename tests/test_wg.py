import os
import pytest
from sqlalchemy import select, exists
from sqlalchemy.sql.functions import count
from sqlalchemy.sql.operators import in_op as in_

import pymlst
from pymlst.common import exceptions
from pymlst.wg import model
from pymlst.wg.core import DatabaseWG, DuplicationHandling

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
    try:
        yield db
    finally:
        db.close()


@pytest.fixture()
def db_simple(db):
    seqs = [
        ('g1', 'AAA'),
        ('g1', 'ATA'),
        ('g2', 'TTT'),
        ('g3', 'CCC'),
        ('g4', 'CCC'),
    ]
    for gene, seq in seqs:
        db.add_genome(gene, 'A', seq)
    return db


@pytest.fixture()
def db_many(db):
    seqs_ref = [
        ('g1', 'AAA'),
        ('g2', 'ATA'),
        ('g3', 'TTT'),
        ('g4', 'CCC'),
        ('g5', 'GGG'),
    ]
    for gene, seq in seqs_ref:
        db.add_core_genome(gene, seq)
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
        db.add_genome(gene, strain, seq)
    return db


def test_add_genome(db):
    db.add_genome('g1', 'A', 'AAA')
    seq = db.connection.execute(
        select([model.sequences])
    ).fetchall()
    assert len(seq) == 1
    assert seq[0].sequence == 'AAA'
    mlst = db.connection.execute(
        select([model.mlst])
    ).fetchall()
    assert len(mlst) == 1
    assert (mlst[0].gene == 'g1'
            and mlst[0].souche == 'A'
            and mlst[0].seqid == seq[0].id)


def test_add_core_genome(db):
    added = db.add_core_genome('g1', 'AAA')
    assert added
    seq = db.connection.execute(
        select([model.sequences])
    ).fetchone()
    assert seq.sequence == 'AAA'
    mlst = db.connection.execute(
        select([model.mlst])
    ).fetchone()
    assert mlst.souche == db.ref == 'ref'
    assert mlst.gene == 'g1' and mlst.seqid == seq.id


def test_add_core_genome_exist_no_duplication_handle(db):
    db.add_core_genome('g1', 'AAA')
    with pytest.raises(exceptions.DuplicatedGeneSequence):
        db.add_core_genome('g2', 'AAA')


def test_add_core_genome_exist_concatenate_handle(db):
    db.add_core_genome('g1', 'AAA')
    added = db.add_core_genome('g2', 'AAA', DuplicationHandling.CONCATENATE)
    assert not added
    seq = db.connection.execute(
         select([model.sequences.c.id])
         .where(model.sequences.c.sequence == 'AAA')
    ).fetchall()
    assert len(seq) == 1
    mlst = db.connection.execute(
        select([model.mlst.c.gene])
    ).fetchall()
    assert len(mlst) == 1
    assert mlst[0].gene == 'g1;g2'


def test_add_core_genome_exist_remove_handle(db):
    db.add_core_genome('g1', 'AAA')
    added = db.add_core_genome('g2', 'AAA', DuplicationHandling.REMOVE)
    assert not added
    seq = db.connection.execute(
         select([model.sequences])
    ).fetchall()
    assert len(seq) == 0
    mlst = db.connection.execute(
        select([model.mlst])
    ).fetchall()
    assert len(mlst) == 0


def test_add_core_genome_gene_exist(db):
    db.add_core_genome('g1', 'AAA')
    with pytest.raises(exceptions.DuplicatedGeneName):
        db.add_core_genome('g1', 'AAT')


def test_add_genome_with_invalid_gene_name(db):
    with pytest.raises(exceptions.InvalidGeneName):
        db.add_core_genome('g1;', 'AAA')


def test_get_core_genome(db):
    db.add_core_genome('g1', 'AAA')
    db.add_core_genome('g2', 'TTT')
    db.add_genome('g3', 'A', 'CCC')
    core_genome = db.core_genome
    assert core_genome == {
        'g1': 'AAA',
        'g2': 'TTT',
    }


def test_remove_gene(db_simple):
    db_simple.remove_gene('g1')
    mlst_e = db_simple.connection.execute(
        select([model.mlst])
        .where(model.mlst.c.gene == 'g1')
    ).fetchone()
    assert mlst_e is None
    seq_e = db_simple.connection.execute(
        select([model.sequences])
        .where(in_(model.sequences.c.sequence,
                   ['AAA', 'ATA']))
    ).fetchone()
    assert seq_e is None


def test_remove_gene_sequence_still_referenced(db_simple):
    db_simple.remove_gene('g3')
    seq_e = db_simple.connection.execute(
        select([model.sequences])
        .where(model.sequences.c.sequence == 'CCC')
    ).fetchone()
    assert seq_e is not None


def test_remove_gene_from_core_genome_dict(db):
    db.add_core_genome('g1', 'AAA')
    assert 'g1' in db.core_genome
    db.remove_gene('g1')
    assert 'g1' not in db.core_genome


def test_remove_strain(db_many):
    db_many.remove_strain('B')
    mlst_e = db_many.connection.execute(
        select([model.mlst])
        .where(model.mlst.c.souche == 'B')
    ).fetchone()
    assert mlst_e is None
    seq_c = db_many.connection.execute(
        select([count(model.sequences.c.id)])
    ).fetchone()
    assert seq_c[0] == 8  # Removed 1 sequence only


def test_remove_reference_strain_attempt(db):
    db.add_core_genome('g1', 'AAA')
    with pytest.raises(exceptions.ReferenceStrainRemoval):
        db.remove_strain('ref')


def test_contains_souche(db_many):
    assert db_many.contains_souche('B')
    db_many.connection.execute(
        model.mlst.delete()
        .where(model.mlst.c.souche == 'B'))
    assert not db_many.contains_souche('B')


def test_get_gene_sequences_many_strains(db_many):
    g1_seq = db_many.get_gene_sequences('g1')
    assert g1_seq == [
        [1, ['A', 'C'], 'AAA'],
        [7, ['B'], 'AAT'],
    ]
    g2_seq = db_many.get_gene_sequences('g2')
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


def test_get_all_strains(db_many):
    strains = db_many.get_all_strains()
    assert strains == ['A', 'B', 'C', 'D']


def test_get_core_genes(db_many):
    genes = db_many.get_core_genes()
    assert genes == ['g1', 'g2', 'g3', 'g4', 'g5']


def test_count_sequences_per_gene(db_many):
    seq_c = db_many.count_sequences_per_gene()
    assert seq_c == {
        'g1': 2,
        'g2': 1,
        'g3': 2,
        'g4': 3,
        'g5': 1
    }


def test_count_souches_per_gene(db_many):
    str_c = db_many.count_souches_per_gene()
    assert str_c == {
        'g1': 3,
        'g2': 2,
        'g3': 3,
        'g4': 4,
        'g5': 1
    }


def test_count_genes_per_souche(db_many):
    gene_c = db_many.count_genes_per_souche(['g1', 'g2', 'g3', 'g4', 'g5'])
    assert gene_c == {
        'A': 4,
        'B': 5,
        'C': 3,
        'D': 1,
        'ref': 5
    }


def test_count_sequences(db_many):
    seq_c = db_many.count_sequences()
    assert seq_c == 9


def test_get_strains_distances(db_many):
    distances = db_many.get_strains_distances(['g1', 'g2', 'g3', 'g4', 'g5'])
    assert distances == {
        'A': {
            'A': 0,
            'B': 3,
            'C': 2,
            'D': 1,
        },
        'B': {
            'A': 3,
            'B': 0,
            'C': 2,
            'D': 0,
        },
        'C': {
            'A': 2,
            'B': 2,
            'C': 0,
            'D': 1,
        },
        'D': {
            'A': 1,
            'B': 0,
            'C': 1,
            'D': 0,
        },
    }


def test_get_mlst(db_many):
    mlst = db_many.get_mlst(['g1', 'g2', 'g3', 'g4', 'g5'])
    assert mlst == {
        'g1': {
            'A': '1',
            'B': '7',
            'C': '1',
        },
        'g2': {
            'A': '2',
            'B': '2',
        },
        'g3': {
            'A': '6',
            'B': '3',
            'C': '3',
        },
        'g4': {
            'A': '4',
            'B': '8',
            'C': '9',
            'D': '8',
        },
        'g5': {
            'B': '5',
        },
    }
