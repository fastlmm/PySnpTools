import os
import numpy as np
import logging
from pysnptools.snpreader import SnpReader
from pysnptools.util import Pairs as UtilPairs
import doctest
import unittest
from pysnptools.kernelreader import SnpKernel
from pysnptools.standardizer import Unit
from fastlmm.association import single_snp #!!!is this wanted?
from pysnptools.util.mapreduce1.runner import LocalMultiProc
import multiprocessing      
import datetime

class Pairs(SnpReader):
    
    def __init__(self, snpreader0,snpreader1=None, do_standardize=True,sid_materialize_limit=1000*1000,_include_single_times_single=False): #!!!cmk could add option to change snp separator and another to encode chrom, etc in the snp name
        super(Pairs, self).__init__()
        self._ran_once = False
        self.snpreader0 = snpreader0
        self.snpreader1 = snpreader1
        self.do_standardize = do_standardize
        self.sid_materialize_limit = sid_materialize_limit
        self._include_single_times_single=_include_single_times_single
        self._utilpairs = UtilPairs(snpreader0.sid,snpreader1.sid if snpreader1 is not None else snpreader0.sid,include_singles=_include_single_times_single,duplicates_ok=False)

    def __repr__(self):
        part2 = '' if self.snpreader1 is None else ',{0}'.format(self.snpreader1)
        return "{0}({1}{2})".format(self.__class__.__name__,self.snpreader0,part2)

    @property
    def row(self):
        """*same as* :attr:`iid`
        """
        if not hasattr(self,"_row"):
            self._row = self.snpreader0.row
            assert self.snpreader1 is None or np.array_equal(self._row,self.snpreader1.row), "Expect snpreaders to have the same iids in the same order"
        return self._row

    @property
    def col(self):
        """*same as* :attr:`sid`
        """
        if not hasattr(self,"_col"):
            assert self.col_count < self.sid_materialize_limit, '{:,} is too many sids to materialize'.format(self.col_count)
            snpreader1 = self.snpreader1 if self.snpreader1 is not None else self.snpreader0
            #!!!cmkself.index0_list = self.snpreader0.sid_to_index((sid0 for sid0,sid1 in self._utilpairs[:])) #!!!cmk can we do without these?
            #!!!cmkself.index1_list = snpreader1.sid_to_index(sid1 for sid0,sid1 in self._utilpairs[:])#!!!cmk can we do without these?
            self._col = np.array(list(','.join(pair) for pair in self._utilpairs[:]))
        return self._col

    @property
    def col_count(self):
        return len(self._utilpairs)

    @property
    def col_property(self):
        """*same as* :attr:`pos`
        """
        if not hasattr(self,"_col_property"):
            self._col_property = np.zeros([self.sid_count,3],dtype=np.int64)
        return self._col_property

    def copyinputs(self, copier):
        # doesn't need to self.run_once() because only uses original inputs !!!is this true?
        self.snpreader0.copyinputs(copier)
        if self.snpreader1 is not None:
            self.snpreader1.copyinputs(copier)

    def run_once(self):
        if self._ran_once:
            return

        self._ran_once = True
        self.col

    def _read(self, iid_index_or_none, sid_index_or_none, order, dtype, force_python_only, view_ok):
        self.run_once()

        iid_count_in = self.iid_count
        sid_count_in = self.sid_count

        if iid_index_or_none is not None:
            iid_count_out = len(iid_index_or_none)
            iid_index_out = iid_index_or_none
        else:
            iid_count_out = iid_count_in
            iid_index_out = range(iid_count_in)

        if sid_index_or_none is not None:
            sid_count_out = len(sid_index_or_none)
            sid_index_out = sid_index_or_none
        else:
            sid_count_out = sid_count_in
            sid_index_out = splice(None)#cmk test this

        pair_array = np.array(list(self._utilpairs[sid_index_out])) #!!!cmk make more efficient with npfromiter?
        snpreader1 = self.snpreader1 if self.snpreader1 is not None else self.snpreader0
        sid_index_inner_0 = self.snpreader0.sid_to_index(pair_array[:,0])
        sid_index_inner_1 = snpreader1.sid_to_index(pair_array[:,1])


        if self.snpreader1 is None:
            sid_index_inner_01 = np.unique(np.r_[sid_index_inner_0,sid_index_inner_1]) #Index of every snp of interest
            inner_01 = self.snpreader0[iid_index_or_none,sid_index_inner_01].read(order=order, dtype=dtype, force_python_only=force_python_only, view_ok=True) #read every val of interest
            val_inner_01 = inner_01.standardize().val if self.do_standardize else inner_01.val

            sid_index_inner_01_reverse = {v:i for i,v in enumerate(sid_index_inner_01)} #Dictionary of snp_index to position in sid_index_inner_01
            sid_index_inner_0_in_val = np.array([sid_index_inner_01_reverse[i] for i in sid_index_inner_0])  #Replace snp_index0 with column # in val_inner_01
            sid_index_inner_1_in_val = np.array([sid_index_inner_01_reverse[i] for i in sid_index_inner_1])  #Replace snp_index1 with column # in val_inner_01
            val_inner_0 = val_inner_01[:,sid_index_inner_0_in_val] #Extract the vals for the left snps of interest
            val_inner_1 = val_inner_01[:,sid_index_inner_1_in_val]#Extract the vals for the right snps of interest
        else:
            inner_0 = self.snpreader0[iid_index_or_none,sid_index_inner_0].read(order=order, dtype=dtype, force_python_only=force_python_only, view_ok=True) #read every val of interest
            inner_1 = self.snpreader1[iid_index_or_none,sid_index_inner_1].read(order=order, dtype=dtype, force_python_only=force_python_only, view_ok=True) #read every val of interest
            val_inner_0 = inner_0.standardize().val if self.do_standardize else inner_0.val
            val_inner_1 = inner_1.standardize().val if self.do_standardize else inner_1.val
        val = val_inner_0*val_inner_1 #Element multiplication creates the vals for the pairs
        return val

#!!!cmk keep these?

def split_on_sids(snpreader,part_count):
    sid_count = snpreader.sid_count
    start = 0
    for part_index in range(1,part_count+1): #!!!cmk make work better in Python 2
        end = part_index*sid_count//part_count
        yield snpreader[:,start:end]
        start=end

def epi_reml(pair_snps,pheno,covar=None,kernel_snps=None,output_dir='results',part_count=33,runner=None,override=False):
    part_list = list(split_on_sids(pair_snps,part_count))
    part_pair_count = (part_count*part_count+part_count)/2
    part_pair_index = -1
    print("part_pair_count={0:,}".format(part_pair_count))
    K0 = SnpKernel(kernel_snps or pair_snps,standardizer=Unit()).read() #Precompute the similarity
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    start_time = datetime.datetime.now()
    for i in range(part_count): #!!!cmk python 2
        part_i = part_list[i]        
        for j in range(i,part_count):
            part_pair_index+=1
            pairs = Pairs(part_i) if i==j else Pairs(part_i,part_list[j])
            print("Looking at pair {0},{1} which is {2} of {3}".format(i,j,part_pair_index,part_pair_count))
            output_file = '{0}/result.{1}.{2}.tsv'.format(output_dir,part_pair_index,part_pair_count)
            if override or not os.path.exists(output_file):
                result_df_ij = single_snp(pairs, K0=K0, pheno=pheno, covar=covar, leave_out_one_chrom=False, count_A1=True,runner=runner)
                result_df_ij.to_csv(output_file, sep="\t", index=False)
                print(result_df_ij[:1])
                time_so_far = datetime.datetime.now()-start_time
                total_time_estimate = time_so_far*part_pair_count/(part_pair_index+1)
                print(total_time_estimate)


    
class TestPairs(unittest.TestCase):

    def test_run1(self):
        from pysnptools.snpreader import Bed
        old_dir = os.getcwd()
        os.chdir(r'D:\OneDrive\programs\epireml\epireml') #!!!cmk make this work without this
        runner = None
        #runner = LocalMultiProc(multiprocessing.cpu_count(),just_one_process=False)

        bed_original = Bed('../syndata.bed',count_A1=False)        
        pheno= '../pheno.txt'
        covar= '../cov.txt' #set to None, if none
        output_dir = r'm:/deldir/results.run1'
        part_count = 2

        epi_reml(bed_original[:,-20:],pheno,covar=covar,output_dir=output_dir,part_count=part_count,runner=runner,override=True)        
        #!!!cmk check answer?
        os.chdir(old_dir)

    def test_run2(self):
        from pysnptools.snpreader import Bed
        old_dir = os.getcwd()
        os.chdir(r'D:\OneDrive\programs\epireml\epireml') #!!!cmk make this work without this
        runner = None
        #runner = LocalMultiProc(multiprocessing.cpu_count(),just_one_process=False)
        bed_original = Bed(r'../syndata.bed',count_A1=False) #Read only the first 10 SNPs
        pheno= r'../pheno.txt'
        covar= r'../cov.txt'
        output_dir = r'm:/deldir/results.run2'
        part_count = 1

        epi_reml(bed_original[:,:20],pheno,kernel_snps=bed_original,covar=covar,output_dir=output_dir,part_count=part_count,runner=runner,override=True)

        #!!!cmk check answer? against M:\deldir\refresults.run1 and 2
        os.chdir(old_dir)

def getTestSuite():
    """
    set up composite test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPairs))
    return test_suite


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    TestPairs().test_run1()
    TestPairs().test_run2()

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=True)
    r.run(suites)

    result = doctest.testmod()
    assert result.failed == 0, "failed doc test: " + __file__


    print("done")

