import numpy as np
import logging
import os
import unittest
import doctest
from pysnptools.distreader import DistReader

class DistGen(DistReader):
    '''
    A :class:`.DistReader` that generates -- on the fly -- deterministic "random" distributions over SNP data..

    See :class:`.DistReader` for general examples of using DistReader.

       **Tip:** The 'gen' in `DistGen` stands for generate, because it generates random (genetic) data. The 'gen' in :class:`Bgen` stands for 'genetic'. *

    **Constructor:**
        :Parameters: * **seed** (*number*) -- The random seed that (along with *block_size*) determines the SNP values.
        :Parameters: * **iid_count** (*number*) --  the number of iids (number of individuals)
        :Parameters: * **sid_count** (*number*) --  the number of sids (number of SNPs)
        :Parameters: * **chrom_count** (*number*) --  the number of chromosomes to generate. Default is 22. Must be 22 or less.
        :Parameters: * **cache_file** (*string*) -- (default None) If provided, tells where to cache the common iid, sid, and pos information. Using it can save time.
        :Parameters: * **block_size** (*number* or None) -- (default None) Tells how many SNP to generate at once. The default picks a *block_size* such that *iid_count \* *block_size* is about 100,000.
        
        :Example:

        >>> from pysnptools.distreader import DistGen
        >>> #Prepare to generate data for 1000 individuals and 1,000,000 SNPs
        >>> dist_gen = DistGen(seed=332,iid_count=1000,sid_count=1000*1000)
        >>> print(dist_gen.iid_count,dist_gen.sid_count)
        1000 1000000
        >>> dist_data = dist_gen[:,200*1000:201*1000].read() #Generate data for all users and for SNPs 200K to 201K
        >>> print(dist_data.val[2,2], dist_data.iid_count, dist_data.sid_count)
        [0.2760024  0.19624834 0.52774925] 1000 1000

    '''

    def __init__(self, seed, iid_count, sid_count, chrom_count=22, cache_file=None, block_size=None):#LATER fails on 450000x1
        self._ran_once = False
        self._cache_file = cache_file

        self._seed = seed
        self._iid_count = iid_count
        self._sid_count = sid_count
        self._chrom_count = chrom_count
        self._block_size = block_size or max((100*1000)//max(1,iid_count),1)

        if cache_file is not None:
            if not os.path.exists(cache_file):
                import pysnptools.util as pstutil
                pstutil.create_directory_if_necessary(cache_file)
                self._run_once()
                np.savez(cache_file, _row=self._row, _col=self._col, _col_property=self._col_property)
            else:
                with np.load(cache_file,allow_pickle=True) as data:
                    self._row = data['_row']
                    assert len(self._row) == iid_count, "The iid in the cache file has a different length than iid_count"
                    self._col = data['_col']
                    assert len(self._col) == sid_count, "The sid in the cache file has a different length than sid_count"
                    self._col_property = data['_col_property']
                    assert len(self._col_property) == sid_count, "The pos in the cache file has a different length than sid_count"
                    self._ran_once = True


    def __repr__(self): 
        return "{0}(seed={1},iid_count={2},sid_count={3},chrom_count={4},block_size={5},cache_file={6})".format(self.__class__.__name__,self._seed,self._iid_count,self._sid_count,self._chrom_count,self._block_size,self._cache_file)

    @property
    def row(self):
        self._run_once()
        return self._row

    @property
    def col(self):
        self._run_once()
        return self._col

    @property
    def col_property(self):
        self._run_once()
        return self._col_property

    _chrom_size = np.array([263,255,214,203,194,183,171,155,145,144,144,143,114,109,106,98,92,85,67,72,50,56],dtype=np.int64)*int(1e6) #The approximate size of human chromosomes in base pairs

    def _run_once(self):
        if (self._ran_once):
            return
        self._ran_once = True
        self._row = np.array([('0','iid_{0}'.format(i)) for i in range(self._iid_count)])
        self._col = np.array(['sid_{0}'.format(i) for i in range(self._sid_count)])
        self._col_property = np.zeros(((self._sid_count),3)) #Must be zero (not empty) because will leave _col_property[:,1] unchanged

        chrom_total = DistGen._chrom_size[:self._chrom_count].sum()
        step = chrom_total // self._sid_count
        start = 0
        chrom_size_so_far = 0
        for chrom_index in range(self._chrom_count):
            chrom_size_so_far += DistGen._chrom_size[chrom_index]
            stop = chrom_size_so_far * self._sid_count // chrom_total
            self._col_property[start:stop,0] = chrom_index+1
            random_increment = np.random.RandomState(chrom_index ^ 2292).randint(step)
            self._col_property[start:stop,2] = np.arange(0,stop-start)*step+random_increment+1
            #print(chrom_index+1,start,stop,self._sid_count)
            start = stop

    def copyinputs(self, copier):
        self._run_once()
        copier.input(self._cache_file)

    # Most _read's support only indexlists or None, but this one supports Slices, too.
    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok, num_threads):
        self._run_once()
        import pysnptools.util as pstutil

        if order == 'A':
            order = 'F'
        dtype = np.dtype(dtype)

        row_index_count = len(row_index_or_none) if row_index_or_none is not None else self._iid_count # turn to a count of the index positions e.g. all of them
        col_index = col_index_or_none if col_index_or_none is not None else np.arange(self._sid_count) # turn to an array of index positions, e.g. 0,1,200,2200,10
        batch_index = col_index // self._block_size  #find the batch index of each index position, e.g. 0,0,0,2,0
        val = np.empty((row_index_count,len(col_index),3),order=order, dtype=dtype) #allocate memory for result
        list_batch_index = list(set(batch_index))
        for ii,i in enumerate(list_batch_index):  #for each distinct batch index, generate dists #!!!fix up snpgen this way, too with ii
            i = int(i) # convert np.uint64 to python int to avoid uint*int->float
            #LATER logging.info("working on distgen batch {0} of {1}".format(ii,len(list_batch_index))) #!!!why does this produce messages like 'working on distgen batch 8 of 2'?
            start = i*self._block_size  #e.g. 0 (then 2000)
            stop = start + self._block_size #e.g. 1000, then 3000
            batch_val = self._get_val(start,stop,dtype) # generate whole batch
            a = (batch_index==i) #e.g. [True,True,True,False,True], then [False,False,False,True,False]
            b = col_index[a]-start #e.g.  0,1,200,10, then 200
            val[:,a,:] = batch_val[:,b,:] if row_index_or_none is None else pstutil.sub_matrix(batch_val, row_index_or_none, b, num_threads=num_threads)

        return val

    def _get_dist(self):
        w = np.array([-0.6482249 , -8.49790398])
    
        def dist_fit(x):
            log_y = w[0] * np.log(x) + w[1]
            return np.exp(log_y)
    
        x_sample = np.logspace(np.log10(.1/self._iid_count),np.log10(.5),100,base=10) #discretize from 1e-7 (assuming iid_count=1e6) to .5, logarithmically
        y_sample = np.array([dist_fit(x) for x in x_sample])                    #Find the relative weight of each point
        dist = y_sample/y_sample.sum()
        return x_sample, dist



    @staticmethod
    def _get_sid(sid_start, sid_stop):
        sid = ["sid_{0}".format(i) for i in range(sid_start,sid_stop)]
        return sid

    def _get_val(self, sid_start, sid_stop,dtype):
        missing_rate = .218
        block_size = sid_stop-sid_start

        np.random.seed(self._seed+sid_start)
        val = np.random.random((self.iid_count,block_size,3))
        val /= val.sum(axis=2,keepdims=True)  #make probabilities sum to 1

        missing = np.random.rand(self.iid_count,block_size)<missing_rate
        val[missing,:] = np.nan

        return val


class TestDistGen(unittest.TestCase):     

    def test1(self):
        logging.info("in TestDistGen test1")
        seed = 0
        distgen = DistGen(seed=seed,iid_count=1000,sid_count=1000*1000,block_size=1000)
        distdata = distgen[:,[0,1,200,2200,10]].read()
        distdata2 = distgen[:,[0,1,200,2200,10]].read()
        assert(distdata.allclose(distdata2))

        from pysnptools.distreader import DistNpz
        ref = DistNpz( os.path.dirname(os.path.realpath(__file__))+'/../../tests/datasets/distgen.dist.npz').read()
        assert(distdata.allclose(ref,equal_nan=True))

        cache_file = 'tempdir/cache_file_test1.npz'
        os.remove(cache_file) if os.path.exists(cache_file) else None
        distgen3 = DistGen(seed=seed,iid_count=1000,sid_count=1000*1000,block_size=1000,cache_file=cache_file)
        distdata3 = distgen3[::10,[0,1,200,2200,10]].read()
        assert(distdata3.allclose(distdata2[::10,:].read()))
        distgen4 = DistGen(seed=seed,iid_count=1000,sid_count=1000*1000,block_size=1000,cache_file=cache_file)
        distdata4 = distgen4[::10,[0,1,200,2200,10]].read()
        assert(distdata4.allclose(distdata2[::10,:].read()))

    def test_doctest(self):
        import pysnptools.distreader.distgen as gen_mod
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(gen_mod)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__


def getTestSuite():
    """
    set up composite test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestDistGen))
    return test_suite


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=True)
    ret = r.run(suites)
    assert ret.wasSuccessful()


    result = doctest.testmod()
    assert result.failed == 0, "failed doc test: " + __file__

