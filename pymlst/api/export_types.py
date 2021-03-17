import pandas as pd
import numpy as np

from pymlst.api.extractors import ExportType


class StrainExport(ExportType):
    def export(self, data, base, output, logger):
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
    def export(self, data, base, output, logger):
        output.write("\n".join(sorted(data.valid_schema)) + "\n")

    @staticmethod
    def name():
        return 'gene'


class DistanceExport(ExportType):
    def export(self, data, base, output, logger):
        if data.duplicate is False:
            logger.info("WARNINGS : Calculate distance between strains " +
                        "using duplicate genes could reported bad result\n")
        output.write(str(len(data.strains)) + "\n")
        distance = base.get_strains_distances(data.ref, data.valid_schema)
        for s1 in data.strains:
            output.write(s1 + "\t")
            c = [str(distance.get(s1, {}).get(s2, 0)) for s2 in data.strains]
            output.write("\t".join(c) + "\n")

    @staticmethod
    def name():
        return 'distance'


class MlstExport(ExportType):
    def export(self, data, base, output, logger):
        output.write("GeneId\t" + "\t".join(data.strains) + "\n")
        mlst = base.get_mlst(data.ref, data.valid_schema)
        for g in data.valid_schema:
            towrite = [g]
            mlstg = mlst.get(g, {})
            for s in data.strains:
                towrite.append(mlstg.get(s, ""))
            output.write("\t".join(towrite) + "\n")

    @staticmethod
    def name():
        return 'mlst'


class GrapetreeExport(ExportType):
    def export(self, data, base, output, logger):
        mlst = base.get_mlst(data.ref, data.valid_schema)
        df = pd.DataFrame(columns=["#GeneId"] + data.strains)
        for g in data.valid_schema:
            row = {"#GeneId": g}
            mlstg = mlst.get(g, {})
            for s in data.strains:
                row[s] = mlstg.get(s, np.NaN)
            df = df.append(row, ignore_index=True)
        df = df.set_index('#GeneId')
        df = df.transpose()
        df = df.fillna(-1).astype(int)
        df.to_csv(output, sep='\t')

    @staticmethod
    def name():
        return 'grapetree'


class StatExport(ExportType):
    def export(self, data, base, output, logger):
        output.write("Strains\t" + str(len(data.strains)) + "\n")
        output.write("Coregenes\t" + str(len(data.all_genes)) + "\n")
        output.write("Sequences\t" + str(base.get_sequences_number(data.ref)) + "\n")

    @staticmethod
    def name():
        return 'stat'
