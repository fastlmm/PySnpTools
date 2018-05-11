import numpy as np
import subprocess, sys, os.path
from itertools import *
import pandas as pd
import logging
from snpreader import SnpReader
from snpdata import SnpData
import math
import warnings
from pysnptools.pstreader import PstData
from pgenlib import PgenReader


class PGen(SnpReader):
    '''
    A :class:`.SnpReader` for random-access reads of pgen/Bim/Fam files from disk.

    See :class:`.SnpReader` for details and examples.

    The format is described in http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml.

    **Constructor:**
        :Parameters: * **filename** (*string*) -- The \*.bed file to read. The '.bed' suffix is optional. The related \*.bim and \*.fam files will also be read.
                     * **count_A1** (*bool*) -- Tells if it should count the number of A1 alleles (the PLINK standard) or the number of A2 alleles.
                                                False is the current default, but in the future the default will change to True.

                     *The following options are never needed, but can be used to avoid reading large '.fam' and '.bim' files when their information is already known.*

                     * **iid** (an array of strings) -- The :attr:`.SnpReader.iid` information. If not given, reads info from '.fam' file.
                     * **sid** (an array of strings) -- The :attr:`.SnpReader.sid` information. If not given, reads info from '.bim' file.
                     * **pos** (optional, an array of strings) -- The :attr:`.SnpReader.pos` information.  If not given, reads info from '.bim' file.
                     * **skip_format_check** (*bool*) -- If False (default), will check that '.bed' file has expected starting bytes.

    **Methods beyond** :class:`.SnpReader`
    '''

    def __init__(self, filename, count_A1=None, iid=None, sid=None, pos=None, skip_format_check=False): #!!!document these new optionals. they are here
        super(PGen, self).__init__()

        self._ran_once = False
        self._file_pointer = None

        self.filename = filename
        if count_A1 is None:
             warnings.warn("'count_A1' was not set. For now it will default to 'False', but in the future it will default to 'True'", FutureWarning)
             count_A1 = False
        self.count_A1 =count_A1
        self.skip_format_check = skip_format_check
        if iid is not None:
            self._row = PstData._fixup_input(iid,empty_creator=lambda ignore:np.empty([0,2],dtype=str))
        if sid is not None:
            self._col = PstData._fixup_input(sid,empty_creator=lambda ignore:np.empty([0],dtype=str))
        if pos is not None:
            self._col_property = PstData._fixup_input(pos,count=len(self._col),empty_creator=lambda count:np.array([[np.nan, np.nan, np.nan]]*count))

    def __repr__(self): 
        return "{0}('{1}',count_A1={2})".format(self.__class__.__name__,self.filename,self.count_A1)

    @property
    def row(self):
        """*same as* :attr:`iid`
        """
        if not hasattr(self,"_row"):
            self._row = SnpReader._read_fam(self.filename,remove_suffix="pgen")
        return self._row

    @property
    def col(self):
        """*same as* :attr:`sid`
        """
        if not hasattr(self,"_col"):
            self._col, self._col_property = SnpReader._read_map_or_bim(self.filename,remove_suffix="pgen", add_suffix="bim", max_filesize=1E9)
        return self._col

    @property
    def col_property(self):
        """*same as* :attr:`pos`
        """
        if not hasattr(self,"_col_property"):
            self._col, self._col_property = SnpReader._read_map_or_bim(self.filename,remove_suffix="pgen", add_suffix="bim", max_filesize=1E9)
        return self._col_property

    def _open_pgen(self, iid_index_or_none=None):
        pgenfile = SnpReader._name_of_other_file(self.filename,"pgen","pgen")
        self._filepointer = PgenReader(pgenfile, raw_sample_ct = None, variant_ct=None, sample_subset=iid_index_or_none)#iid_index_or_none)
        
    def _close_pgen(self):
        self.__del__()
        self.file_pointer = None

    def _run_once(self):
        if self._ran_once:
            return
        self._ran_once = True

        if not hasattr(self,"_row"):
            self._row = SnpReader._read_fam(self.filename,remove_suffix="pgen")

        if not hasattr(self,"_col") or not hasattr(self,"_col_property"):
            self._col, self._col_property = SnpReader._read_map_or_bim(self.filename,remove_suffix="pgen", add_suffix="bim")
        self._assert_iid_sid_pos()

        if not self.skip_format_check:
            self._open_pgen()
            self._close_pgen()

    def __del__(self):
        if hasattr(self,'_filepointer') and self._filepointer is not None:  # we need to test this because Python doesn't guarantee that __init__ was fully run
            self._filepointer.close()
            self._filepointer = None

    def _read(self, iid_index_or_none, sid_index_or_none, order, dtype, force_python_only, view_ok, out_buffer):
        if force_python_only:
            raise NotImplementedError("No Python only mode")
        self._run_once()

        if order=='A':
            order='F'
        
        assert not hasattr(self, 'ind_used'), "A SnpReader should not have a 'ind_used' attribute"

        iid_count_in = self.iid_count
        sid_count_in = self.sid_count

        if iid_index_or_none is not None:
            iid_count_out = len(iid_index_or_none)
            iid_index_out = np.array(iid_index_or_none, dtype=np.uint32)
        else:
            iid_count_out = iid_count_in
            iid_index_out = None

        if sid_index_or_none is not None:
            sid_count_out = len(sid_index_or_none)
            sid_index_out = sid_index_or_none
        else:
            sid_count_out = sid_count_in
            sid_index_out = range(sid_count_in)

        out_compatible = (order=="F") and ((dtype == np.int8) or (dtype == np.int32) or (dtype == np.int64))
        
        if out_buffer is not None:
            assert out_buffer.dtype == dtype, "Wrong type"
            assert out_buffer.flags["F_CONTIGUOUS"] == (order=='F'), "wrong order"
            assert out_buffer.shape[0] == iid_count_out, "insufficient first dimension"
            assert out_buffer.shape[1] == sid_count_out, "insufficient second dimension"
            
        if out_compatible and (out_buffer is not None):
            val = out_buffer.T
            dtype_val = out_buffer.dtype
        elif out_compatible:
            dtype_val = dtype
            val = np.empty((sid_count_out, iid_count_out), order='C', dtype=dtype_val)
            out_buffer = val.T  
        else:
            dtype_val = np.int8
            val = np.empty((sid_count_out, iid_count_out), order='C', dtype=dtype_val)

        pgen_fn = SnpReader._name_of_other_file(self.filename,"pgen","pgen")

        self._open_pgen(iid_index_out)
        self._filepointer.read_list(np.array(sid_index_out, dtype=np.uint32), val)
        self._close_pgen()

        if (out_buffer is None) and out_compatible:
            out_buffer = val.T
        elif (out_buffer is not None) and (not out_compatible):
            out_buffer[:] = val.T
        elif (out_buffer is None) and (not out_compatible):
            out_buffer = np.array(val.T,  dtype=dtype, order=order)
        if (out_buffer.dtype == np.float64) or (out_buffer.dtype == np.float32):
            out_buffer[out_buffer==-9.0] = np.nan
        return out_buffer # [:iid_count_out,:sid_count_out]


if __name__ == "__main__":
    pgenfile = "/Users/clippert/Source/cerebro/CodeBase/Src/Pgenlib/example_data/qc-export-20160229-NODEV-VarMaj_Hg38_Mm.pgen"
    reader = PGen(filename=pgenfile, count_A1=None, iid=None, sid=None, pos=None, skip_format_check=False)
    l=reader[:,10000:20000].read()

    if 1:
        import time
        stepsize = 10000
        t0 = time.time()
        start = 0
        stop = 0
        for i in xrange(100):
            start = stop
            stop = start + stepsize
            ll = reader[:,start:stop].read(out_buffer=None)
        t1 = time.time()
        diff = t1 - t0

        print "time unbuffered: ", diff

        t0_ = time.time()
        start = 0
        stop = 0
        val = np.zeros((l.val.shape[0], stepsize), dtype=np.float64, order="F")
        for i in xrange(100):
            start = stop
            stop = start + stepsize
            ll = reader[:,start:stop].read(out_buffer=val)
        t1_ = time.time()
        diff_ = t1_- t0_
        print "time   buffered: ", diff_

        print l.val
