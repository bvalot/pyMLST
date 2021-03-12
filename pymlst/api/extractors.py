import os
import subprocess
import tempfile

from pymlst.api.wgmlst import Extractor

def run_mafft(path, tmpfile):
    command = [path, '--quiet', tmpfile.name]
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    genes = {}
    ids = None
    seq = ""
    for line in iter(proc.stdout.readline, ''):
        if ids is None and line[0] == '>':
            ids = line.lstrip('>').rstrip("\n")
        elif ids is not None and line[0] == '>':
            genes[int(ids)] = seq.upper()
            ids = line.lstrip('>').rstrip("\n")
            seq = ""
        elif ids is not None:
            seq += line.rstrip("\n")
        else:
            raise Exception("A problem occurred while running mafft" + str(line))
    if seq != "":
        genes[int(ids)] = seq.upper()

    return genes


def write_tmp_seqs(tmpfile, seqs):
    tmp = open(tmpfile.name, 'w+t')
    for s in seqs:
        tmp.write(">"+str(s[0])+"\n"+s[2]+"\n")
    tmp.close()


def add_sequence_strain(seqs, strains, sequences):
    """Add a sequence to multi-align, take the first gene in case of repetition"""
    size = 0
    if len(seqs) > 0:
        size = len(seqs[0][2])
    for s in strains:
        seq = [i[2] for i in seqs if s in i[1]]
        if len(seq) == 0:
            sequences.get(s).append('-' * size)
        elif len(seq) == 1:
            sequences.get(s).append(seq[0])
        else:
            raise Exception("repeated gene must be excluded in order to align export\n")


class SequenceExtractor(Extractor):

    def __init__(self, list, align, realign, mincover):
        self.list = list
        self.align = align
        self.realign = realign
        self.mincover = mincover
        self.mafft_path = '/usr/bin/mafft'

    def extract(self, base, ref, output, logger):
        tmpfile = tempfile.NamedTemporaryFile(mode='w+t', suffix='.fasta', delete=False)
        tmpfile.close()
        try:
            # Minimun number of strain
            strains = [i[0] for i in base.get_different_souches(ref)]
            if self.mincover < 1 or self.mincover > len(strains):
                raise Exception("Mincover must be between 1 to number of strains : " + str(len(strains)))

            #  Coregene
            coregene = []
            if self.list is not None:
                coregene = [l.rstrip("\n") for l in iter(self.list.readline, '')]
            else:
                coregene = [l[0] for l in base.get_genes_coverages(ref) if l[1] >= self.mincover]

            logger.info("Number of gene to analyse : " + str(len(coregene)) + "\n")

            if self.align is False:
                # no multialign
                for g in coregene:
                    seqs = base.get_gene_sequences(g, ref)
                    for seq in seqs:
                        output.write(">" + g + "|" + str(seq[0]) + " "
                                     + ";".join(seq[1]) + "\n")
                        output.write(seq[2] + "\n")
            else:
                # multialign
                # search duplicate
                dupli = base.get_duplicated_genes(ref)
                for d in dupli:
                    logger.info('Duplicated: ', d)

                sequences = {s: [] for s in strains}
                for i, g in enumerate(coregene):
                    logger.info(str(i + 1) + "/" + str(len(coregene)) + " | " + g + "     ")
                    if g in dupli:
                        logger.info("No: Repeat gene\n")
                        continue
                    seqs = base.get_gene_sequences(g, ref)
                    size = set()
                    for seq in seqs:
                        size.add(len(seq[2]))
                    if len(size) == 1 and self.realign is False:
                        logger.info("Direct")
                        add_sequence_strain(seqs, strains, sequences)
                    else:
                        logger.info("Align")
                        write_tmp_seqs(tmpfile, seqs)
                        corrseqs = run_mafft(self.mafft_path, tmpfile)
                        for seq in seqs:
                            seq[2] = corrseqs.get(seq[0])
                        add_sequence_strain(seqs, strains, sequences)
                    logger.info("\n")

                # output align result
                for s in strains:
                    output.write('>' + s + "\n")
                    output.write("\n".join(sequences.get(s)) + "\n")
        finally:
            if os.path.exists(tmpfile.name):
                os.remove(tmpfile.name)
