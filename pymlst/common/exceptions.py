class PyMLSTError(Exception):
    pass


class GeneError(PyMLSTError):
    pass


class DuplicatedGeneSequence(GeneError):
    pass


class DuplicatedGeneName(GeneError):
    pass


class InvalidGeneName(GeneError):
    pass


class AlleleSequenceNotFound(PyMLSTError):
    pass


class WrongBaseType(PyMLSTError):
    pass


class ReferenceStrainRemoval(PyMLSTError):
    pass


class BadIdentityRange(PyMLSTError):
    pass

class BadCoverageRange(PyMLSTError):
    pass


class BinaryNotFound(PyMLSTError):
    pass


class StrainAlreadyPresent(PyMLSTError):
    pass


class ChromosomeNotFound(PyMLSTError):
    pass


class CoreGenomePathNotFound(PyMLSTError):
    pass


class NothingToRemove(PyMLSTError):
    pass


class UndefinedExportType(PyMLSTError):
    pass

class EmptyDatabase(PyMLSTError):
    pass
