import numpy as np
import logging
import pysnptools.util as pstutil
from bed_reader import get_num_threads, standardize_f64, standardize_f32


class Standardizer(object):
    '''
    A Standardizer is a class such as :class:`.Unit` and :class:`.Beta` to be used by the :meth:`.SnpData.standardize` and :meth:`.SnpReader.read_kernel` method to standardize SNP data.

    :Example:

    Read and standardize SNP data.

    >>> from pysnptools.standardizer import Unit
    >>> from pysnptools.snpreader import Bed
    >>> from pysnptools.util import example_file # Download and return local file name
    >>> bedfile = example_file("tests/datasets/all_chr.maf0.001.N300.*","*.bed")
    >>> snpdata1 = Bed(bedfile,count_A1=False).read().standardize(Unit())
    >>> print('{0:.6f}'.format(snpdata1.val[0,0]))
    0.229416

    Create a kernel from SNP data on disk.

    >>> bedfile2 = example_file("pysnptools/examples/toydata.5chrom.*","*.bed")
    >>> kerneldata = Bed(bedfile2,count_A1=False).read_kernel(Unit())
    >>> print('{0:.6f}'.format(kerneldata.val[0,0]))
    9923.069928


    Can also return a constant SNP standardizer that can be applied to other :class:`.SnpData`.

    >>> snp_whole = Bed(bedfile,count_A1=False)
    >>> train_idx, test_idx = range(10,snp_whole.iid_count), range(0,10) #test on the first 10, train on the rest
    >>> snp_train, trained_standardizer = Unit().standardize(snp_whole[train_idx,:].read(),return_trained=True)
    >>> print('{0:.6f}'.format(snp_train.val[0,0]))
    0.233550
    >>> print(trained_standardizer.stats[0,:]) #The mean and stddev of the 1st SNP on the training data # '...' for a possible space character
    [...1.94827586  0.22146953]
    >>> snp_test = snp_whole[test_idx,:].read().standardize(trained_standardizer)
    >>> snp_test.val[0,0]
    0.23354968324845735

    Standardize any Numpy array.

    >>> val = Bed(bedfile,count_A1=False).read().val
    >>> print(val[0,0])
    2.0
    >>> val = Unit().standardize(val)
    >>> print('{0:.6f}'.format(val[0,0]))
    0.229416

    Details of Methods & Properties:
    '''

    def __init__(self):
        super(Standardizer, self).__init__()

    def standardize(self, snps, block_size=None, return_trained=False, force_python_only=False, num_threads=None):
        '''
        Applies standardization, in place, to :class:`.SnpData` (or a NumPy array). For convenience also returns the :class:`.SnpData` (or a NumPy array).

        :param snps: SNP values to standardize
        :type snps: :class:`.SnpData` (or a NumPy array)

        :param block_size: *Not used*
        :type block_size: None

        :param return_trained: If true, returns a second value containing a constant standardizer trained on this data.
        :type return_trained: bool

        :param force_python_only: optional -- If False (default), may use outside library code. If True, requests that the read
            be done without outside library code.
        :type force_python_only: bool

        :param num_threads: optional -- The number of threads with which to standardize data. Defaults to all available
            processors. Can also be set with these environment variables (listed in priority order):
            'PST_NUM_THREADS', 'NUM_THREADS', 'MKL_NUM_THREADS'.
        :type num_threads: None or int

        :rtype: :class:`.SnpData` (or a NumPy array), (optional) constant :class:`.Standardizer`

        '''
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)
        raise NotImplementedError("subclass {0} needs to implement method '.standardize'".format(self.__class__.__name__))

    @staticmethod
    #changes snps in place
    def _standardize_unit_and_beta(snps, is_beta, a, b, apply_in_place, use_stats, stats, num_threads, force_python_only=False):
        '''
        When snps is a cupy ndarray, will use cupy to compute new stats for unit. (Other paths are not defined for cupy)
        '''
        xp = pstutil.get_array_module(snps)

        assert snps.flags["C_CONTIGUOUS"] or snps.flags["F_CONTIGUOUS"], "Expect snps to be order 'C' or order 'F'"

        #Make sure stats is the same type as snps. Because we might be creating a new array, we return it
        if stats is None:
            stats = xp.empty([snps.shape[1],2],dtype=snps.dtype,order="F" if snps.flags["F_CONTIGUOUS"] else "C")
        elif not (
             stats.dtype == snps.dtype   #stats must have the same dtype as snps
             and (stats.flags["OWNDATA"]) # stats must own its data
             and (snps.flags["C_CONTIGUOUS"] and stats.flags["C_CONTIGUOUS"]) or (snps.flags["F_CONTIGUOUS"] and stats.flags["F_CONTIGUOUS"]) #stats must have the same order as snps
             ):
            stats = xp.array(stats,dtype=snps.dtype,order="F" if snps.flags["F_CONTIGUOUS"] else "C")
        assert stats.shape == (snps.shape[1],2), "stats must have size [sid_count,2]"

        if not force_python_only and xp is np:
            num_threads = get_num_threads(num_threads)

            if snps.dtype == np.float64:
                if (snps.flags['F_CONTIGUOUS'] or snps.flags['C_CONTIGUOUS']) and (snps.flags["OWNDATA"] or snps.base.nbytes == snps.nbytes): #!!create a method called is_single_segment
                    standardize_f64(snps,is_beta,a,b,apply_in_place,use_stats,stats,num_threads)
                    return stats
                else:
                    logging.info("Array is not contiguous, so will standardize with python only instead of C++")
            elif snps.dtype == np.float32:
                if (snps.flags['F_CONTIGUOUS'] or snps.flags['C_CONTIGUOUS']) and (snps.flags["OWNDATA"] or snps.base.nbytes == snps.nbytes):
                    standardize_f32(snps,is_beta,a,b,apply_in_place,use_stats,stats,num_threads)
                    return stats
                else:
                    logging.info("Array is not contiguous, so will standardize with python only instead of C++")
            else:
                logging.info("Array type is not float64 or float32, so will standardize with python only instead of C++")

        import pysnptools.standardizer as stdizer
        if is_beta:
            Standardizer._standardize_beta_python(snps, a, b, apply_in_place, use_stats=use_stats, stats=stats)
            return stats
        else:
            Standardizer._standardize_unit_python(snps, apply_in_place, use_stats=use_stats, stats=stats)
            return stats

    @staticmethod
    def _standardize_unit_python(snps,apply_in_place,use_stats,stats):
        '''
        Standardize snps to zero-mean and unit variance.

        Will work with both numpy and cupy ndarray.
        '''
        assert snps.dtype in [np.float64,np.float32], "snps must be a float in order to standardize in place."
        xp = pstutil.get_array_module(snps)

        imissX = xp.isnan(snps)

        if use_stats:
            snp_mean = stats[:,0]
            snp_std = stats[:,1]
        else:
            snp_std = xp.nanstd(snps, axis=0)
            snp_mean = xp.nanmean(snps, axis=0)
            # avoid div by 0 when standardizing
            #Don't need this warning because SNCs are still meaning full in QQ plots because they should be thought of as SNPs without enough data.
            #logging.warn("A least one snps has only one value, that is, its standard deviation is zero")
            snp_std[snp_std == 0.0] = xp.inf #We make the stdev infinity so that applying as a trained_standardizer will turn any input to 0. Thus if a variable has no variation in the training data, then it will be set to 0 in test data, too. 
            stats[:,0] = snp_mean
            stats[:,1] = snp_std

        if apply_in_place:
            snps -= snp_mean
            snps /= snp_std
            snps[imissX] = 0

    @property
    def is_constant(self):
        '''
        Tell if a standardizer's statistics are constant. For example, :class:`.Unit_Trained` contains a pre-specified mean and stddev for each SNP, so :attr:`Unit_Trained.is_constant` is True.
        :class:`.Unit`, on the other hand, measures the mean and stddev for each SNP, so :attr:`Unit.is_constant` is False.

        :rtype: bool
        '''
        return False        

    @staticmethod
    def _standardize_beta_python(snps, betaA, betaB, apply_in_place, use_stats, stats):
        '''
        standardize snps with Beta prior
        '''
        assert snps.dtype in [np.float64,np.float32], "snps must be a float in order to standardize in place."

        imissX = np.isnan(snps)
        snp_sum =  np.nansum(snps,axis=0)
        n_obs_sum = (~imissX).sum(0)
    
        if use_stats:
            snp_mean = stats[:,0]
            snp_std = stats[:,1]
        else:
            snp_mean = (snp_sum*1.0)/n_obs_sum
            snp_std = np.sqrt(np.nansum((snps-snp_mean)**2, axis=0)/n_obs_sum)
            if 0.0 in snp_std:
                #Don't need this warning because SNCs are still meaning full in QQ plots because they should be thought of as SNPs without enough data.
                #logging.warn("A least one snps has only one value, that is, its standard deviation is zero")
                snp_std[snp_std==0] = np.inf
            stats[:,0] = snp_mean
            stats[:,1] = snp_std

        if apply_in_place:
            maf = snp_mean/2.0
            maf[maf>0.5]=1.0 - maf[maf>0.5]

            # avoid div by 0 when standardizing
            import scipy.stats as st
            maf_beta = st.beta.pdf(maf, betaA, betaB)
            #print("BetaPdf[{0},{1},{2}]={3}".format(maf,betaA,betaB,maf_beta))
            snps -= snp_mean
            snps*=maf_beta
            snps[imissX] = 0.0
            if use_stats: #If we're applying to test data, set any variables with to 0 if they have no variation in the training data.
                snps[:,snp_std==np.inf] = 0.0

    def _merge_trained(self, trained_list):
        raise Exception("Not defined")

class _CannotBeTrained(Standardizer):

    def __init__(self, name):
        super(_CannotBeTrained, self).__init__()
        self.name=name

    def __repr__(self): 
        return "{0}({1})".format(self.__class__.__name__,self.name)

    def standardize(self, snps, block_size=None, return_trained=False, force_python_only=False, num_threads=None):
        if block_size is not None:
            warnings.warn("block_size is deprecated (and not needed, since standardization is in-place", DeprecationWarning)
        raise Exception("Standardizer '{0}' cannot be trained",self)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS|doctest.NORMALIZE_WHITESPACE)
