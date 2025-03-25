import pytest
from sqlalchemy import select

import pymlst
from pymlst.pytyper import model
from pymlst.pytyper.method import FIM, SPA, CLMT
from pymlst.pytyper.core import DatabaseTyper, TypingResult, FimH, Spa, Clmt
from pymlst.common import exceptions


@pytest.fixture()
def fim():
    with pymlst.open_typer(FIM) as fim_typer:
        yield fim_typer

@pytest.fixture()
def spa():
    with pymlst.open_typer(SPA) as spa_typer:
        yield spa_typer

@pytest.fixture()
def clmt():
    with pymlst.open_typer(CLMT) as clmt_typer:
        yield clmt_typer


@pytest.fixture()
def db():
    db = DatabaseTyper(None)
    try:
        yield db
    finally:
        db.close()

@pytest.fixture()
def result():
    res = TypingResult('sample1', FIM)
    yield res

@pytest.fixture()
def db_many(db):
    seqs = [  # seq, typing, allele
        ('AAA', 'fim', 'fimH1'),
        ('ATA', 'fim', 'fimH2'),
        ('TTT', 'fim', 'fimH3'),
        ('CCC', 'spa', '01'),
        ('CCT', 'spa', '02'),
        ('AAT', 'spa', '03'),
        ('ATT', 'clmt', 'arpA'),
        ('TCT', 'clmt', 'chuA'),
        ('ACC', 'clmt', 'yjaA'),
        ('CTC', 'clmt', 'TspE4.C2'),
    ]
    for seq, method, allele in seqs:
        db.add_sequence(seq, method, allele)
    sts = [  # st, typing, allele
        ('fimH1', 'fim', 'fimH1'),
        ('fimH2', 'fim', 'fimH2'),
        ('fimH3', 'fim', 'fimH3'),
        ('t1', 'spa', '01-02-02-01'),
        ('t2', 'spa', '02-02-03-01'),
        ('t3', 'spa', '01-01-02-01'),
        ('A', 'clmt', 'arpA|+,chuA|-,yjaA|-,TspE4.C2|-'),
        ('B1', 'clmt', 'arpA|+,chuA|-,yjaA|-,TspE4.C2|+'),
        ('B2', 'clmt', 'arpA|-,chuA|+,yjaA|+,TspE4.C2|+'),
        ('D|E', 'clmt', 'arpA|+,chuA|+,yjaA|-,TspE4.C2|-'),
    ]
    for st, method, allele in sts:
        db.add_st(st, method, allele)
    return db

def test_check_db(db_many):
    res = db_many.check_db(FIM)
    assert(res) == True

def test_check_new_db(db):
    res = db.check_db(CLMT)
    assert(res) == False

def test_add_sequence(db):
    db.add_sequence('AAA', FIM, '02')    
    seq = db.connection.execute(
        select(model.typerSeq)
    ).fetchall()
    assert len(seq) == 1
    assert (seq[0].sequence == 'AAA'
            and seq[0].typing == FIM
            and seq[0].allele == '02')

def test_add_st(db):
    db.add_st('fimH1', FIM, 'fimH1')    
    st = db.connection.execute(
        select(model.typerSt)
    ).fetchall()
    assert len(st) == 1
    assert (st[0].st == 'fimH1'
            and st[0].typing == FIM
            and st[0].allele == 'fimH1')

def test_get_sequences(db_many):
    seqs = db_many.get_sequences(FIM)
    seqs2 = db_many.get_sequences(CLMT)
    assert len(seqs) == 3
    assert len(seqs2) == 4
    assert seqs[1] == ('fimH2', 'ATA')
    assert seqs2[0] == ('arpA','ATT')

def test_get_sequence_by_allele(db_many):
    seq = db_many.get_sequence_by_allele(FIM, 'fimH1')
    assert seq == 'AAA'
    with pytest.raises(exceptions.AlleleSequenceNotFound):
        db_many.get_sequence_by_allele(SPA, '04')

def test_get_allele_by_sequence(db_many):
    al = db_many.get_allele_by_sequence(SPA, 'CCT')
    assert al == '02'
    al2 = db_many.get_allele_by_sequence(FIM, 'GCG')
    assert al2 == 'New'

def test_get_st(db_many):
    st = db_many.get_st(FIM, 'fimH2')
    assert st == 'fimH2'
    st2 = db_many.get_st(CLMT, '02')
    assert st2 == ''
    
def test_pyTyper_instance(fim, spa, clmt):
    assert isinstance(fim, FimH)
    assert isinstance(spa, Spa)
    assert isinstance(clmt, Clmt)

def test_pyTyper_check_input(fim):
    a = fim.check_input(0.9, 0.9)
    with pytest.raises(exceptions.BadCoverageRange):
        fim.check_input(0.9, 18)
    with pytest.raises(exceptions.BadIdentityRange):
        fim.check_input(12, 0.2)

def test_typingResult_full(result):
    result.set_allele('12')
    result.set_st('t1235')
    result.set_notes('Some informations')
    assert str(result) == 'sample1 fim t1235 12'

def test_typingResult_empty(result):
    assert str(result) == 'sample1 fim  '
