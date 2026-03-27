#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2019 Benoit Valot
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

from Bio.Data.CodonTable import TranslationError

from pymlst.common import mafft

class Psl:
    """A simple Psl class"""
    def __init__(self, pslline):
        pslelement = pslline.rstrip("\n").split("\t")
        if len(pslelement) != 21:
            raise Exception("Psl line have not 21 elements:\n"+pslline)
        self.pslelement = pslelement
        self.chro = pslelement[13]
        self.start = int(pslelement[15])
        self.end = int(pslelement[16])
        self.strand = pslelement[8]
        self.rstart = int(pslelement[11])
        self.rend = int(pslelement[12])
        self.rtotal = int(pslelement[10])
        self.coverage = (float(self.rend) - self.rstart)/self.rtotal

    def gene_id(self):
        return self.pslelement[9]

    def get_sequence(self, seq):
        if self.strand == '+':
            return seq[self.start:self.end]
        return seq[self.start:self.end].reverse_complement()

    def get_aligned_sequence(self, seq, coregene):
        if self.strand == '+':
            expand_start = self.rstart > 0
            expand_end = self.rend < self.rtotal
        else:
            expand_start = self.rend < self.rtotal
            expand_end = self.rstart > 0

        if expand_start:
            start = self.start - 36
            if start < 0:
                start = 0
        else:
            start = self.start

        if expand_end:
            end = self.end + 36
            if end > len(seq):
                end = len(seq)
        else:
            end = self.end

        target = seq[start:end]
        if self.strand != '+':
            target = target.reverse_complement()

        al_start, al_end = mafft.get_aligned_area(coregene, str(target))
        if al_start is not None:
            return target[al_start:al_end]

        return None

