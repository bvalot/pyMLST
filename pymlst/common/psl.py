#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2019 Benoit Valot
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

from Bio.Data.CodonTable import TranslationError

from pymlst.common import mafft


def test_cds(seq):
    try:
        seq.translate(table="Bacterial", cds=True)
    except TranslationError:
        return False
    else:
        return True

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

    # def searchCorrect(self):
    #     if int(self.pslelement[11]) != 0:
    #         diff = int(self.pslelement[11])
    #         if self.strand == "+":
    #             self.start = self.start - diff
    #         else:
    #             self.end = self.end + diff
    #     elif int(self.pslelement[10]) != int(self.pslelement[12]):
    #         diff = int(self.pslelement[10]) - int(self.pslelement[12])
    #         if self.strand == "+":
    #             self.end = self.end + diff
    #         else:
    #             self.start = self.start - diff
    #     self.coverage = 1
    #
    # def searchCorrectCDS(self, seq, coverage):
    #     prot = self.get_sequence(seq)
    #     ##modifs start and stop not create
    #     if prot.startswith("M") is False and prot.endswith("*") is False:
    #         return False
    #     windows = int((1-coverage)*self.rtotal)
    #     if prot.startswith("M") is False:
    #         return self.__searchCDS(seq, True, False, windows, 0)
    #     elif prot.endswith("*") is False:
    #         return self.__searchCDS(seq, False, True, windows, 0)
    #     else:
    #         raise Exception("A problem of start/stop  for gene " + self.gene_id())

    # def searchPartialCDS(self, seq, coverage):
    #     ##modifs start and stop not create
    #     if self.rstart !=0 and self.rend != self.rtotal:
    #         return False
    #     windows = int((1-coverage)*self.rtotal)
    #     if self.rstart !=0:
    #         diff = self.rstart
    #         return self.__searchCDS(seq, True, False, windows, diff)
    #     elif self.rend != self.rtotal:
    #         diff = self.rtotal - self.rend
    #         return self.__searchCDS(seq, False, True, windows, diff)
    #     else:
    #         raise Exception("A problem of start/stop for gene " + self.gene_id())

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

    # def __searchCDS(self, seq, start, stop, windows, diff):
    #     ##correct windows/diff multiple of 3
    #     windows = windows - windows%3
    #     diff = diff - diff%3
    #     ##modifs start and stop not create
    #     if start and stop:
    #         return False
    #     ##modifs start
    #     if start:
    #         ##modulo = (self.end-self.start)%3
    #         if self.strand == "+":
    #             theoStart = self.__getTheoricStart(diff)
    #             val = [i for i in range(theoStart+windows, theoStart-windows, -3) \
    #                    if test_cds(seq.seq[i:self.end], False)]
    #             if len(val)==1:
    #                 self.start=val[0]
    #                 return True
    #             elif len(val) >1:
    #                 best = self.__getBest(val)
    #                 self.logger.info("Choosing best start for gene " + self.gene_id() + " " \
    #                                  + str(best) + " " + str(val))
    #                 self.start = best
    #                 return True
    #             else:
    #                 return False
    #         else:
    #             theoEnd = self.__getTheoricEnd(diff)
    #             val = [i for i in range(theoEnd-windows, theoEnd+windows, 3) \
    #                    if test_cds(seq.seq[self.start:i], True)]
    #             if len(val) == 1:
    #                 self.end = val[0]
    #                 return True
    #             elif len(val) >1:
    #                 best = self.__getBest(val)
    #                 self.logger.info("Choosing best start for gene " + self.gene_id() + " " \
    #                                  + str(best) + " " + str(val))
    #                 self.end = best
    #                 return True
    #             else:
    #                 return False
    #     ##modifs end
    #     elif stop:
    #         ##modulo = (self.end-self.start)%3
    #         if self.strand == "+":
    #             theoEnd = self.__getTheoricEnd(diff)
    #             val = [i for i in range(theoEnd-windows, theoEnd+windows, 3) \
    #                    if test_cds(seq.seq[self.start:i], False)]
    #             if len(val) == 1:
    #                 self.end = val[0]
    #                 return True
    #             else:
    #                 return False
    #         else:
    #             theoStart = self.__getTheoricStart(diff)
    #             val = [i for i in range(theoStart+windows, theoStart-windows, -3) \
    #                    if test_cds(seq.seq[i:self.end], True)]
    #             if len(val) == 1:
    #                 self.start = val[0]
    #                 return True
    #             else:
    #                 return False
    #
    # def __getTheoricStart(self, diff):
    #     modulo = (self.end-self.start)%3
    #     return self.start + modulo - diff
    #
    # def __getTheoricEnd(self, diff):
    #     modulo = (self.end-self.start)%3
    #     return self.end - modulo + diff
    #
    # def __getBest(self, val):
    #     best = val[0]
    #     for v in val[1:]:
    #         if self.strand == "+":
    #             if abs(abs(self.end - v) - self.rtotal) < abs(abs(self.end - best) - self.rtotal):
    #                 best = v
    #         else:
    #             if abs(abs(v - self.start) - self.rtotal)
    #             < abs(abs(best - self.start) - self.rtotal):
    #                 best = v
    #     return best
