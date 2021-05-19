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


class BadIdentityRange(PyMLSTError):
    pass


class BinaryNotFound(PyMLSTError):
    pass


class StrainAlreadyPresent(PyMLSTError):
    pass


class ChromosomeNotFound(PyMLSTError):
    pass


class CoreGenomePathNotFound(PyMLSTError):
    pass
