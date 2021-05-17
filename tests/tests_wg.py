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
    ]
    for strain, gene, seq in seqs:
        _, seq_id = db.add_sequence(seq)
        db.add_mlst(strain, gene, seq_id)
    return db


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
        select([exists(db_simple.sequences)])
        .where(in_(db_simple.sequences.c.id, [1, 3]))
    ).scalar()
    assert not seq_e
    mlst_e = db_simple.connection.execute(
        select([exists(db_simple.mlst)])
        .where(in_(db_simple.mlst.c.seqid, [1, 3]))
    ).scalar()
    assert not mlst_e


def test_remove_gene(db_simple):
    db_simple.remove_gene('g1')
    mlst_e = db_simple.connection.execute(
        select([exists(db_simple.mlst)])
        .where(db_simple.mlst.c.gene == 'g1')
    ).scalar()
    assert not mlst_e


def test_remove_orphan_sequences(db_simple):
    db_simple.remove_gene('g2')
    db_simple.remove_gene('g4')
    db_simple.remove_orphan_sequences([3, 4])
    seq_1_e = db_simple.connection.execute(
        select([exists(db_simple.sequences)])
        .where(db_simple.sequences.c.id == 3)
    ).scalar()
    assert not seq_1_e
    seq_2_e = db_simple.connection.execute(
        select([exists(db_simple.sequences)])
        .where(db_simple.sequences.c.id == 4)
    ).scalar()
    assert seq_2_e


def test_remove_strain(db_simple):
    db_simple.remove_strain('A')
    mlst_e = db_simple.connection.execute(
        select([exists(db_simple.mlst.c.id)])
        .where(db_simple.mlst.c.souche == 'A')
    ).scalar()
    assert not mlst_e


def test_get_gene_sequences_ids(db_simple):
    ids = db_simple.get_gene_sequences_ids('g1')
    assert ids == [1, 2]


def test_get_strain_sequences_ids(db_many):
    ids = db_many.get_strain_sequences_ids('A')
    assert ids == {1, 2, 3, 4}


def test_get_sequence_by_gene_and_souche(db_many):
    seq = db_many.get_sequence_by_gene_and_souche('g1', 'B')
    assert seq.sequence == 'AAT'
