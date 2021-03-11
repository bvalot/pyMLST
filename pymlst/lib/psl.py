#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2019 Benoit Valot
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

import sys

def testCDS(seq, reverse):
    try:
        if reverse:
            seq.reverse_complement().translate(table="Bacterial", cds=True)
        else:
            seq.translate(table="Bacterial", cds=True)
    except:
        return False
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
        # if self.coverage !=1 and self.coverage>=0.95:
        #     self.correct()
        
    def geneId(self):
        return self.pslelement[9]

    def getSequence(self, seq):
        if self.strand =="+":
            sequence = seq.seq[self.start:self.end]
        else:
            sequence = seq.seq[self.start:self.end].reverse_complement()
        # ##Verify sequence correct
        # if len(sequence) != (self.end-self.start):
        #     raise Exception("Gene " + self.geneId() + " incomplete\n")
        return sequence

    def searchCorrect(self):
        if int(self.pslelement[11]) != 0:
            diff = int(self.pslelement[11])
            if self.strand == "+":
                self.start = self.start - diff
            else:
                self.end = self.end + diff
        elif int(self.pslelement[10]) != int(self.pslelement[12]):
            diff = int(self.pslelement[10]) - int(self.pslelement[12])
            if self.strand == "+":
                self.end = self.end + diff
            else:
                self.start = self.start - diff   
        self.coverage = 1
    
    def searchCorrectCDS(self, seq, coverage):
        prot = self.getSequence(seq)
        ##modifs start and stop not create
        if prot.startswith("M") is False and prot.endswith("*") is False:
            return False
        windows = int((1-coverage)*self.rtotal)
        if prot.startswith("M") is False:
            return self.__searchCDS(seq, True, False, windows, 0)
        elif prot.endswith("*") is False:
            return self.__searchCDS(seq, False, True, windows, 0)
        else:
            raise Exception("A problem of start/stop  for gene " + self.geneId())

    def searchPartialCDS(self, seq, coverage):
        ##modifs start and stop not create
        if self.rstart !=0 and self.rend != self.rtotal:
            return False
        windows = int((1-coverage)*self.rtotal)
        if self.rstart !=0:
            diff = self.rstart
            return self.__searchCDS(seq, True, False, windows, diff)
        elif self.rend != self.rtotal:
            diff = self.rtotal - self.rend
            return self.__searchCDS(seq, False, True, windows, diff)
        else:
            raise Exception("A problem of start/stop for gene " + self.geneId())
    
    def __searchCDS(self, seq, start, stop, windows, diff):
        ##correct windows/diff multiple of 3
        windows = windows - windows%3
        diff = diff - diff%3
        ##modifs start and stop not create
        if start and stop:
            return False
        ##modifs start
        if start:
            ##modulo = (self.end-self.start)%3
            if self.strand == "+":
                theoStart = self.__getTheoricStart(diff)
                val = [i for i in range(theoStart+windows, theoStart-windows, -3) \
                       if testCDS(seq.seq[i:self.end], False)]
                if len(val)==1:
                    self.start=val[0]
                    return True
                elif len(val) >1:
                    best = self.__getBest(val)
                    sys.stderr.write("Choice best start for gene " + self.geneId() + " " \
                                     + str(best) + " " + str(val) + "\n")
                    self.start = best
                    return True
                else:
                    return False
            else:
                theoEnd = self.__getTheoricEnd(diff)
                val = [i for i in range(theoEnd-windows, theoEnd+windows, 3) \
                       if testCDS(seq.seq[self.start:i], True)]
                if len(val) == 1:
                    self.end = val[0]
                    return True
                elif len(val) >1:
                    best = self.__getBest(val)
                    sys.stderr.write("Choice best start for gene " + self.geneId() + " " \
                                     + str(best) + " " + str(val) + "\n")
                    self.end = best
                    return True
                else:
                    return False
        ##modifs end
        elif stop:
            ##modulo = (self.end-self.start)%3
            if self.strand == "+":
                theoEnd = self.__getTheoricEnd(diff)
                val = [i for i in range(theoEnd-windows, theoEnd+windows, 3) \
                       if testCDS(seq.seq[self.start:i], False)]
                if len(val) == 1:
                    self.end = val[0]
                    return True
                else:
                    return False
            else:
                theoStart = self.__getTheoricStart(diff)
                val = [i for i in range(theoStart+windows, theoStart-windows, -3) \
                       if testCDS(seq.seq[i:self.end], True)]
                if len(val) == 1:
                    self.start = val[0]
                    return True
                else:
                    return False

    def __getTheoricStart(self, diff):
        modulo = (self.end-self.start)%3
        return self.start + modulo - diff

    def __getTheoricEnd(self, diff):
        modulo = (self.end-self.start)%3
        return self.end - modulo + diff
        
    def __getBest(self, val):
        best = val[0]
        for v in val[1:]:
            if self.strand == "+":
                if abs(abs(self.end - v) - self.rtotal) < abs(abs(self.end - best) - self.rtotal):
                    best = v
            else:
                if abs(abs(v - self.start) - self.rtotal) < abs(abs(best - self.start) - self.rtotal):
                    best = v
        return best

