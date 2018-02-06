from .procedure import IntermediateProcedureStep
from ..misc.log import logprocedure

class ShearCatEqualSubdivision(object):

    def __init__(self, num_input_bins):
        pass

class ShearCatRebinning(IntermediateProcedureStep):
    """
    This step rebins the shear catalogs - it usually increases or decreases 
    (but does not completely consolidate) the number of tomographic bins.

    In principle, this should also support discarding unecessary tomographic 
    information.

    For now, this is implemented with the assumption that individual bins will 
    be subdivided. In the future, if more complex rebinning is required, we can 
    refactor.
    """

    def __init__(self):
        pass
