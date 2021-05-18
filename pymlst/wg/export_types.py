import logging

import pandas as pd
import numpy as np

from pymlst.wg.extractors import ExportType


class StrainExport(ExportType):
    def export(self, data, base, output):
        if data.count is False:
            output.write("\n".join(data.strains) + "\n")
        else:
            tmp = base.count_genes_per_souche(data.valid_schema)
            for strain in data.strains:
                output.write(strain + "\t" + str(tmp.get(strain)) + "\n")

    @staticmethod
    def name():
        return 'strain'


class GeneExport(ExportType):
    def export(self, data, base, output):
        output.write("\n".join(sorted(data.valid_schema)) + "\n")

    @staticmethod
    def name():
        return 'gene'


class DistanceExport(ExportType):
    def export(self, data, base, output):
        if data.duplicate:
            logging.info("WARNINGS : Calculate distance between strains ",
                         "using duplicate genes could reported bad result.")
        output.write(str(len(data.strains)) + "\n")
        distance = base.get_strains_distances(data.valid_schema)
        for strain in data.strains:
            output.write(strain + "\t")
            dist = [str(distance.get(strain, {}).get(s2, 0)) for s2 in data.strains]
            output.write("\t".join(dist) + "\n")

    @staticmethod
    def name():
        return 'distance'


class MlstExport(ExportType):
    def export(self, data, base, output):
        output.write("GeneId\t" + "\t".join(data.strains) + "\n")
        mlst = base.get_mlst(data.valid_schema)
        for gene in data.valid_schema:
            towrite = [gene]
            mlstg = mlst.get(gene, {})
            for strain in data.strains:
                towrite.append(mlstg.get(strain, ""))
            output.write("\t".join(towrite) + "\n")

    @staticmethod
    def name():
        return 'mlst'


class GrapetreeExport(ExportType):
    def export(self, data, base, output):
        mlst = base.get_mlst(data.valid_schema)
        strains = pd.DataFrame(columns=["#GeneId"] + data.strains)
        for gene in data.valid_schema:
            row = {"#GeneId": gene}
            mlstg = mlst.get(gene, {})
            for strain in data.strains:
                row[strain] = mlstg.get(strain, np.NaN)
            strains = strains.append(row, ignore_index=True)
        strains = strains.set_index('#GeneId')
        strains = strains.transpose()
        strains = strains.fillna(-1).astype(int)
        strains.to_csv(output, sep='\t')

    @staticmethod
    def name():
        return 'grapetree'


class StatExport(ExportType):
    def export(self, data, base, output):
        output.write("Strains\t" + str(len(data.strains)) + "\n")
        output.write("Coregenes\t" + str(len(data.all_genes)) + "\n")
        output.write("Sequences\t" + str(base.count_sequences()) + "\n")

    @staticmethod
    def name():
        return 'stat'
