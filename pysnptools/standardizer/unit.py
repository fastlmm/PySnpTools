import numpy as np
import scipy as sp
import logging
import pysnptools.standardizer as stdizer

# attempt to import wrapped plink parser
WRAPPED_PLINK_PARSER_PRESENT = True
try:
    import pysnptools.snpreader.wrap_plink_parser as wrap_plink_parser
except Exception:
    WRAPPED_PLINK_PARSER_PRESENT = False



class Unit(object):  #IStandardizer
    """The specification for unit standardization"""
    def __init__(self):
        pass

    def standardize(self, snps, blocksize=None, force_python_only=False):
        l = self.lambdaFactory(snps, blocksize=blocksize, force_python_only=force_python_only)
        return stdizer.standardize_with_lambda(snps, l, blocksize)

    def __repr__(self):
        return "{0}()".format(self.__class__.__name__)

    def lambdaFactory(self, snps, blocksize=None, force_python_only=False):
        if not force_python_only:
            hideSNCWarning = logging.getLogger().getEffectiveLevel() > logging.WARNING # Only show the C code's SNC warning if the logging level is set to WARNING or below

            if snps.dtype == np.float64:
                if snps.flags['F_CONTIGUOUS'] and (snps.flags["OWNDATA"] or snps.base.nbytes == snps.nbytes):
                    #!!! LATER: set snps to np.load(D:\Source\carlk_snpreader\tests\temp.npz.npz) and then run this. It fails without error. Could it be that standardizing twice, sometimes causes this?
                    return lambda s, hideSNCWarning=hideSNCWarning : wrap_plink_parser.standardizedoubleFAAA(s,False,float("NaN"),float("NaN"),hideSNCWarning)
                elif snps.flags['C_CONTIGUOUS'] and (snps.flags["OWNDATA"] or snps.base.nbytes == snps.nbytes) and blocksize is None:
                    return lambda s, hideSNCWarning=hideSNCWarning : wrap_plink_parser.standardizedoubleCAAA(s,False,float("NaN"),float("NaN"),hideSNCWarning)
                else:
                    logging.info("Array is not contiguous, so will standardize with python only instead of C++")
            elif snps.dtype == np.float32:
                if snps.flags['F_CONTIGUOUS'] and (snps.flags["OWNDATA"] or snps.base.nbytes == snps.nbytes):
                    return lambda s, hideSNCWarning=hideSNCWarning: wrap_plink_parser.standardizefloatFAAA(s,False,float("NaN"),float("NaN"),hideSNCWarning)
                elif snps.flags['C_CONTIGUOUS'] and (snps.flags["OWNDATA"] or snps.base.nbytes == snps.nbytes) and blocksize is None:
                    return lambda s, hideSNCWarning=hideSNCWarning: wrap_plink_parser.standardizefloatCAAA(s,False,float("NaN"),float("NaN"),hideSNCWarning)
                else:
                    logging.info("Array is not contiguous, so will standardize with python only instead of C++")
            else:
                logging.info("Array type is not float64 or float32, so will standardize with python only instead of C++")

        return lambda s: stdizer.standardize_unit_python(s)
