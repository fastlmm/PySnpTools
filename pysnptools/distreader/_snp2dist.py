import numpy as np
from itertools import *
import pandas as pd
import logging
import time
from pysnptools.distreader import DistReader
from pysnptools.snpreader import SnpReader
import pysnptools.standardizer as stdizer

class _Snp2Dist(DistReader):

    def __init__(self, snpreader, max_weight=2.0, block_size=None):
        super(_Snp2Dist, self).__init__()

        self.snpreader = snpreader
        self.max_weight = max_weight
        self.block_size = block_size

    @property
    def row(self):
        return self.snpreader.iid

    @property
    def col(self):
        return self.snpreader.col

    def __repr__(self):
        return self._internal_repr()

    def _internal_repr(self): #!!! merge this with __repr__
        s = "{0}.as_dist(".format(self.snpreader)
        if self.block_size is not None:
            s += "block_size={0}".format(self.block_size)
        s += ")"
        return s

    def copyinputs(self, copier):
        #Doesn't need run_once
        copier.input(self.snpreader)

    def _snpval_to_distval(self, snpval, order, dtype):
        if order=='A':
            order='F'
        dtype = np.dtype(dtype)

        distval = np.zeros([snpval.shape[0],snpval.shape[1],3],dtype=dtype,order=order)
        distval = distval.reshape(-1,distval.shape[-1],order=order)
        factor = 2.0/self.max_weight
        for count in range(distval.shape[-1]):
            bool_array = snpval.reshape(-1,order=order)*factor==count
            distval[bool_array,count]=1
        distval[distval.sum(axis=-1)!=1,:]=np.nan #replace any dist that doesn't sum to 1 with nan
        distval = distval.reshape([snpval.shape[0],snpval.shape[1],distval.shape[-1]],order=order)
        return distval

    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok, num_threads):
        assert row_index_or_none is None and col_index_or_none is None #real assert because indexing should already be pushed to the inner snpreader
        dtype = np.dtype(dtype)

        #Do all-at-once (not in blocks) if 1. No block size is given or 2. The #ofSNPs < Min(block_size,iid_count)
        if self.block_size is None or (self.sid_count <= self.block_size or self.sid_count <= self.iid_count):
            snpdata,_ = SnpReader._as_snpdata(self.snpreader,dtype=dtype,order=order,force_python_only=force_python_only,standardizer=stdizer.Identity(),num_threads=num_threads)
            val = self._snpval_to_distval(snpdata.val,order,dtype)

            has_right_order = order="A" or (order=="C" and val.flags["C_CONTIGUOUS"]) or (order=="F" and val.flags["F_CONTIGUOUS"])
            assert has_right_order
            return val
        else: #Do in blocks
            t0 = time.time()
            if order=='A':
                order = 'F'
            val = np.zeros([self.iid_count,self.sid_count,3],dtype=dtype,order=order)#LATER use empty or fillnan???

            logging.info("reading {0} value data in blocks of {1} SNPs and finding distribution (for {2} individuals)".format(self.sid_count, self.block_size, self.iid_count))
            ct = 0
            ts = time.time()

            for start in range(0, self.sid_count, self.block_size):
                ct += self.block_size
                snpdata = self.snpreader[:,start:start+self.block_size].read(order=order,dtype=dtype,force_python_only=force_python_only,view_ok=True,num_threads=num_threads) # a view is always OK, because we'll allocate memory in the next step
                val[:,start:start+self.block_size] = self._snpval_to_distval(snpdata.val,order,dtype)
                if ct % self.block_size==0:
                    diff = time.time()-ts
                    if diff > 5:
                       logging.info("read %s SNPs in %.2f seconds" % (ct, diff))

            t1 = time.time()
            logging.info("%.2f seconds elapsed" % (t1-t0))

            return val


    def __getitem__(self, iid_indexer_and_snp_indexer):
        row_index_or_none, col_index_or_none = iid_indexer_and_snp_indexer
        return _Snp2Dist(self.snpreader[row_index_or_none,col_index_or_none],max_weight=self.max_weight,block_size=self.block_size)


    @property
    def sid(self):
        '''The :attr:`SnpReader.sid` property of the SNP data.
        '''
        return self.snpreader.sid

    @property
    def sid_count(self):
        '''The :attr:`SnpReader.sid_count` property of the SNP data.
        '''
        return self.snpreader.sid_count

    @property
    def pos(self):
        '''The :attr:`SnpReader.pos` property of the SNP data.
        '''
        return self.snpreader.pos


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
