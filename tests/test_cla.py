import pytest
from sqlalchemy import select

import pymlst
from pymlst.cla import model
from pymlst.cla.core import DatabaseCLA
from pymlst.common import exceptions


@pytest.fixture()
def cla():
    with pymlst.open_cla() as cla_mlst:
        yield cla_mlst


@pytest.fixture()
def db():
    db = DatabaseCLA(None, 1)
    try:
        yield db
    finally:
        db.close()


@pytest.fixture()
def db_many(db):
    seqs = [  # gene, seq, allele
        ('g1', 'AAA', 1),
        ('g2', 'ATA', 1),
        ('g3', 'TTT', 1),
        ('g4', 'CCC', 1),
        ('g5', 'CCT', 1),
        ('g1', 'AAT', 2),
        ('g2', 'ATT', 2),
        ('g3', 'TCT', 2),
        ('g4', 'ACC', 2),
        ('g5', 'CTC', 2),
    ]
    for gene, seq, allele in seqs:
        db.add_sequence(seq, gene, allele)
    mlst = [  # st, gene, allele
        (1, 'g1', 1),
        (1, 'g2', 2),
        (1, 'g3', 1),
        (1, 'g4', 2),
        (1, 'g5', 2),
        (2, 'g1', 1),
        (2, 'g2', 1),
        (2, 'g3', 2),
        (2, 'g4', 1),
        (2, 'g5', 1),
        (3, 'g1', 2),
        (3, 'g2', 1),
        (3, 'g3', 1),
        (3, 'g4', 1),
        (3, 'g5', 2),
    ]
    for st, gene, allele in mlst:
        db.add_mlst(st, gene, allele)
    return db


def test_add_sequence(db):
    db.add_sequence('AAA', 'g1', 2)
    seq = db.connection.execute(
        select(model.sequences)
    ).fetchall()
    assert len(seq) == 1
    assert (seq[0].sequence == 'AAA'
            and seq[0].gene == 'g1'
            and seq[0].allele == 2)


def test_add_mlst(db):
    db.add_sequence('AAA', 'g1', 2)
    db.add_mlst(5, 'g1', 2)
    mlst = db.connection.execute(
        select(model.mlst)
    ).fetchall()
    assert len(mlst) == 1
    assert (mlst[0].st == 5
            and mlst[0].gene == 'g1'
            and mlst[0].allele == 2)
    assert len(db.core_genome) == 0


# def test_add_mlst_no_sequence(db):
#     db.add_sequence('AAA', 'g1', 1)
#     with pytest.raises(exceptions.AlleleSequenceNotFound):
#         db.add_mlst(5, 'g1', 2)


def test_add_mlst_reference(db):
    db.add_sequence('AAA', 'g1', 1)
    db.add_mlst(5, 'g1', 1)
    assert len(db.core_genome) == 1
    assert db.core_genome['g1'] == 'AAA'


def test_get_genes_by_allele(db_many):
    genes = db_many.get_genes_by_allele(2)
    assert genes == {
        'g1': 'AAT',
        'g2': 'ATT',
        'g3': 'TCT',
        'g4': 'ACC',
        'g5': 'CTC',
    }


def test_get_allele_by_sequence_and_gene(db_many):
    allele = db_many.get_allele_by_sequence_and_gene('AAT', 'g1')
    assert allele == 2


def test_get_allele_by_sequence_and_gene_none(db_many):
    allele = db_many.get_allele_by_sequence_and_gene('AAT', 'g2')
    assert allele is None


def test_get_st_by_gene_and_allele(db_many):
    st = db_many.get_st_by_gene_and_allele('g3', 1)
    assert st == [1, 3]
    st = db_many.get_st_by_gene_and_allele('g2', 2)
    assert st == [1]


def test_get_st_by_gene_and_allele_none(db_many):
    st = db_many.get_st_by_gene_and_allele('g5', 3)
    assert st == []


def test_get_sequence_by_gene_and_allele(db_many):
    seq = db_many.get_sequence_by_gene_and_allele('g3', 2)
    assert seq == 'TCT'


def test_get_sequence_by_gene_and_allele_none(db_many):
    seq = db_many.get_sequence_by_gene_and_allele('g3', 6)
    assert seq is None
