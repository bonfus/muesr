class SampleException(Exception):
    """
    Exceptions connected with problems with the Sample definition.
    """
    pass
    
class SymmetryError(SampleException):
    """
    Symmetry needed but not defined.
    """
    pass

class CellError(SampleException):
    """
    Lattice structure needed but not defined.
    """
    pass
    
class MagDefError(SampleException):
    """
    Magnetic structure not defined or incompatible.
    """
    pass

class MuonError(SampleException):
    """
    Muon position needed but not defined.
    """
    pass

