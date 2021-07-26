import numpy as np
from itertools import *
import pandas as pd
import logging
from bed_reader import open_bed, to_bed
from pysnptools.snpreader import SnpReader
from pysnptools.snpreader import SnpData
import math
import warnings
from pysnptools.pstreader import PstData


plink_chrom_map = {"X": 23, "Y": 24, "XY": 25, "MT": 26}

reverse_plink_chrom_map = {23:"X", 24:"Y", 25:"XY", 26:"MT"}

class Bed(SnpReader):
    """
    A :class:`.SnpReader` for random-access reads of Bed/Bim/Fam files from disk.

    See :class:`.SnpReader` for details and examples.

    The format is described `here <http://zzz.bwh.harvard.edu/plink/binary.shtml>`__.

    **Constructor:**
        :Parameters: * **filename** (*string*) -- The \*.bed file to read. The '.bed' suffix is optional. The related \*.bim and \*.fam files will also be read.
                     * **count_A1** (*bool*) -- Tells if it should count the number of A1 alleles (the PLINK standard) or the number of A2 alleles. False is the current default, but in the future the default will change to True.

                     *The following three options are never needed, but can be used to avoid reading large '.fam' and '.bim' files when their information is already known.*

                     * **iid** (an array of strings) -- The :attr:`SnpReader.iid` information. If not given, reads info from '.fam' file.
                     * **sid** (an array of strings) -- The :attr:`SnpReader.sid` information. If not given, reads info from '.bim' file.
                     * **pos** (optional, an array of strings) -- The :attr:`SnpReader.pos` information.  If not given, reads info from '.bim' file.
                     * **num_threads** (optinal, int) -- The number of threads with which to read data. Defaults to all available processors.
                            Can also be set with these environment variables (listed in priority order):
                            'PST_NUM_THREADS', 'NUM_THREADS', 'MKL_NUM_THREADS'.
                     * **skip_format_check** (*bool*) -- By default (False), checks that the '.bed' file starts with the expected bytes
                            the first time any file ('.bed', '.fam', or '.bim') is opened.
                     * **fam_filename** (optional, *string*) -- The \*.fam file to read. Defaults to the bed filename with the suffix replaced.
                     * **bim_filename** (optional, *string*) -- The \*.bim file to read. Defaults to the bed filename with the suffix replaced.
                     * **chrom_map** (optional, *dictionary*) -- A dictionary from non-numeric chromosome names to numbers. Defaults to the PLINK
                            mapping, namely, :data:`plink_chrom_map`.

    **Constants**

    * plink_chrom_map
        A dictionary that 'X' to 23, 'Y' to 24, 'XY' to 25, and 'MT' to 26.

    * reverse_plink_chrom_map
        A dictionary that maps 23 to 'X', 24 to 'Y', 25 to 'XY', and 26 to 'MT'

    **Methods beyond** :class:`.SnpReader`

        The :meth:`.SnpReader.read` method returns a :class:`SnpData` with a :attr:`SnpData.val` ndarray. By default, this ndarray is
        numpy.float32. Optionally, it can be numpy.float16. For :class:`Bed`, however, it can also be numpy.int8 with missing values
        represented by -127.

        When reading, any chromosome and position values of 0 (the PLINK standard for missing) will be represented in :attr:`pos` as NaN.

        :Example:

        >>> from pysnptools.snpreader import Bed
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bedfile = example_file("tests/datasets/distributed_bed_test1_X.*","*.bed")
        >>> snp_on_disk = Bed(bedfile, count_A1=False)
        >>> snpdata1 = snp_on_disk.read() # Read into the default, an float64 ndarray
        >>> snpdata1.val.dtype
        dtype('float64')
        >>> snpdata1 = snp_on_disk.read(dtype='int8',_require_float32_64=False) #Read into an 'int8' ndarray.
        >>> snpdata1.val.dtype
        dtype('int8')

    """

    def __init__(
        self,
        filename,
        count_A1=None,
        iid=None,
        sid=None,
        pos=None,
        num_threads=None,
        skip_format_check=False,
        fam_filename=None,
        bim_filename=None,
        chrom_map = plink_chrom_map,
    ):
        super(Bed, self).__init__()

        self._ran_once = False

        self.filename = SnpReader._name_of_other_file(filename,remove_suffix="bed", add_suffix="bed")
        self.fam_filename = fam_filename
        self.bim_filename = bim_filename
        if count_A1 is None:
            warnings.warn(
                "'count_A1' was not set. For now it will default to 'False', but in the future it will default to 'True'",
                FutureWarning,
            )
            count_A1 = False
        self.count_A1 = count_A1
        self._skip_format_check = skip_format_check
        self._original_iid = iid
        self._original_sid = sid
        self._original_pos = pos
        self._num_threads = num_threads
        self._open_bed = None
        self.chrom_map = chrom_map

    def __repr__(self):
        return "{0}('{1}',count_A1={2})".format(
            self.__class__.__name__, self.filename, self.count_A1
        )

    def _open_bed_if_needed(self):
        if self._open_bed is not None:
            return

        properties = {
            "father": None,
            "mother": None,
            "sex": None,
            "pheno": None,
            "allele_1": None,
            "allele_2": None,
        }
        if self._original_iid is not None:
            properties["fid"] = self._original_iid[:, 0]
            properties["iid"] = self._original_iid[:, 1]
        if self._original_sid is not None:
            properties["sid"] = self._original_sid
        if self._original_pos is not None:
            properties["chromosome"] = self._original_pos[:, 0]
            properties["cm_position"] = self._original_pos[:, 1]
            properties["bp_position"] = self._original_pos[:, 2]

        self._open_bed = open_bed(
            self.filename,
            properties=properties,
            skip_format_check=self._skip_format_check,
            count_A1=self.count_A1,
            num_threads=self._num_threads,
            fam_filepath=self.fam_filename,
            bim_filepath=self.bim_filename,
        )

    @property
    def row(self):
        """*same as* :attr:`iid`
        """
        if not hasattr(self, "_row"):
            self._open_bed_if_needed()

            self._row = np.array(
                [self._open_bed.fid, self._open_bed.iid]
            ).T  # LATER: could copy in batches or use concatenate

        return self._row

    @property
    def col(self):
        """*same as* :attr:`sid`
        """
        if not hasattr(self, "_col"):
            self._open_bed_if_needed()
            self._col = self._open_bed.sid
        return self._col

    @property
    def col_property(self):
        """*same as* :attr:`pos`
        """
        if not hasattr(self, "_col_property"):
            self._open_bed_if_needed()

            chromosome = self._open_bed.chromosome
            intersection = self.chrom_map.keys() & chromosome
            for key in intersection:
                chromosome[chromosome==key]=self.chrom_map[key]

            self._col_property = np.array(
                [
                    self._open_bed.chromosome.astype("float"),
                    self._open_bed.cm_position,
                    self._open_bed.bp_position,
                ]
            ).T  # LATER: Could copy in batches to use less memory
            self._col_property[
                self._col_property == 0
            ] = (
                np.nan
            )

        return self._col_property

    def _run_once(self):
        if self._ran_once:
            return
        self._ran_once = True

        self.row
        self.col
        self.col_property
        assert self._open_bed is not None  # real assert

        self._assert_iid_sid_pos(check_val=False)

    def __del__(self):
        pass

    def copyinputs(self, copier):
        # doesn't need to self.run_once() because only uses original inputs
        copier.input(
            SnpReader._name_of_other_file(
                self.filename, remove_suffix="bed", add_suffix="bed"
            )
        )
        copier.input(
            SnpReader._name_of_other_file(
                self.filename, remove_suffix="bed", add_suffix="bim"
            )
        )
        copier.input(
            SnpReader._name_of_other_file(
                self.filename, remove_suffix="bed", add_suffix="fam"
            )
        )

    @staticmethod
    def write(
        filename,
        snpdata,
        count_A1=False,
        force_python_only=False,
        _require_float32_64=True,
        num_threads=None,
        reverse_chrom_map = {},
    ):
        """Writes a :class:`SnpData` to Bed format and returns the :class:`.Bed`.

        :param filename: the name of the file to create
        :type filename: string
        :param snpdata: The in-memory data that should be written to disk.
        :type snpdata: :class:`SnpData`
        :param count_A1: Tells if it should count the number of A1 alleles (the PLINK standard) or the number of A2 alleles. False is the current default, but in the future the default will change to True.
        :type count_A1: bool
        :param force_python_only: Defaults to False. Tells to use Python code rather than faster Rust code. (Might be useful for debugging).
        :type force_python_only: bool
        :param _require_float32_64: Defaults to True. Requires that snpdata's dtype is float32 or float64. (False is useful for writing int8 data.)
        :type _require_float32_64: bool
        :param num_threads: Maximum number of threads to use. Currently ignored and number of threads is 1.
        :type num_threads: int
        :param reverse_chrom_map: Dictionary from chromosome number to chromsome string to write in the \*.bim file. Defaults to empty dictionary.
        :type reverse_chrom_map: dictionary
        :rtype: :class:`.Bed`

        Any :attr:`pos` values of NaN will be written as 0, the PLINK standard for missing chromosome and position values. 

        >>> from pysnptools.snpreader import Pheno, Bed
        >>> import pysnptools.util as pstutil
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bed_fn = example_file("pysnptools/examples/toydata.5chrom.bed")
        >>> snpdata = Bed(bed_fn)[:,::2].read() # Read every-other SNP
        >>> pstutil.create_directory_if_necessary("tempdir/everyother.bed")
        >>> Bed.write("tempdir/everyother.bed",snpdata,count_A1=False)   # Write data in Bed format
        Bed('tempdir/everyother.bed',count_A1=False)
        >>> # Can write from an int8 array, too.
        >>> snpdata_int = SnpData(val=np.int_(snpdata.val).astype('int8'),iid=snpdata.iid,sid=snpdata.sid,pos=snpdata.pos,_require_float32_64=False)
        >>> snpdata_int.val.dtype
        dtype('int8')
        >>> Bed.write("tempdir/everyother.bed",snpdata_int,count_A1=False,_require_float32_64=False)
        Bed('tempdir/everyother.bed',count_A1=False)
        """

        if isinstance(filename, SnpData) and isinstance(
            snpdata, str
        ):  # For backwards compatibility, reverse inputs if necessary
            warnings.warn(
                "write statement should have filename before data to write",
                DeprecationWarning,
            )
            filename, snpdata = snpdata, filename

        if count_A1 is None:
            warnings.warn(
                "'count_A1' was not set. For now it will default to 'False', but in the future it will default to 'True'",
                FutureWarning,
            )
            count_A1 = False

        filename = SnpReader._name_of_other_file(filename,remove_suffix="bed", add_suffix="bed")

        chromosome = snpdata.pos[:, 0]
        intersection = reverse_chrom_map.keys() & chromosome
        if len(intersection) > 0:
            chromosome = chromosome.astype("object")
            for key in intersection:
                chromosome[chromosome==key]=reverse_chrom_map[key]

        to_bed(
            filename,
            val=snpdata.val,
            properties={
                "fid": snpdata.iid[:, 0],
                "iid": snpdata.iid[:, 1],
                "sid": snpdata.sid,
                "chromosome": chromosome,
                "cm_position": snpdata.pos[:, 1],
                "bp_position": snpdata.pos[:, 2],
            },
            count_A1=count_A1,
            force_python_only=force_python_only,
            num_threads=num_threads,
        )

        return Bed(filename, count_A1=count_A1)

    def _read(
        self,
        iid_index_or_none,
        sid_index_or_none,
        order,
        dtype,
        force_python_only,
        view_ok,
        num_threads,
    ):
        self._run_once()

        if order == "A":
            order = "F"

        assert not hasattr(
            self, "ind_used"
        ), "A SnpReader should not have a 'ind_used' attribute"

        val = self._open_bed.read(
            index=(iid_index_or_none, sid_index_or_none),
            order=order,
            dtype=dtype,
            force_python_only=force_python_only,
            num_threads=num_threads
        )

        return val


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    import os

    if False:  # Look for example Bed files with missing data
        from pysnptools.util._example_file import pysnptools_hashdown
        from pysnptools.util import example_file

        for file in pysnptools_hashdown.walk():
            if file.endswith(".bed"):
                print(file + "?")
                bed_file = None
                try:
                    bed_file = example_file(file[:-4] + ".*", "*.bed")
                except Exception:
                    pass
                if bed_file is not None:
                    bed = Bed(bed_file)
                    snpdata = bed[:1000, :1000].read()
                    if not np.all(snpdata.val == snpdata.val):
                        print(bed_file + "!")

    if False:
        from pysnptools.snpreader import Bed
        from pysnptools.util import example_file  # Download and return local file name

        # bed_file = example_file('doc/ipynb/all.*','*.bed')
        bed_file = r"F:\backup\carlk4d\data\carlk\cachebio\genetics\onemil\id1000000.sid_1000000.seed0.byiid\iid990000to1000000.bed"
        bed = Bed(bed_file, count_A1=False)
        snpdata1 = bed[:, :1000].read()
        snpdata2 = bed[:, :1000].read(dtype="int8", _require_float32_64=False)
        print(snpdata2)
        snpdata3 = bed[:, :1000].read(
            dtype="int8", order="C", _require_float32_64=False
        )
        print(snpdata3)
        snpdata3.val = snpdata3.val.astype("float32")
        snpdata3.val.dtype

    if False:
        from pysnptools.snpreader import Bed, SnpGen

        iid_count = 487409
        sid_count = 5000
        sid_count_max = 5765294
        sid_batch_size = 50

        sid_batch_count = -(sid_count // -sid_batch_size)
        sid_batch_count_max = -(sid_count_max // -sid_batch_size)
        snpgen = SnpGen(seed=234, iid_count=iid_count, sid_count=sid_count_max)

        for batch_index in range(sid_batch_count):
            sid_index_start = batch_index * sid_batch_size
            sid_index_end = (batch_index + 1) * sid_batch_size  # what about rounding
            filename = r"d:\deldir\rand\fakeukC{0}x{1}-{2}.bed".format(
                iid_count, sid_index_start, sid_index_end
            )
            if not os.path.exists(filename):
                Bed.write(
                    filename + ".temp", snpgen[:, sid_index_start:sid_index_end].read()
                )
                os.rename(filename + ".temp", filename)

    if False:
        from pysnptools.snpreader import Pheno, Bed

        filename = r"m:\deldir\New folder (4)\all_chr.maf0.001.N300.bed"
        iid_count = 300
        iid = [["0", "iid_{0}".format(iid_index)] for iid_index in range(iid_count)]
        bed = Bed(filename, iid=iid, count_A1=False)
        print(bed.iid_count)

    if False:
        from pysnptools.util import example_file

        pheno_fn = example_file("pysnptools/examples/toydata.phe")

    if False:
        from pysnptools.snpreader import Pheno, Bed
        import pysnptools.util as pstutil
        import os

        print(os.getcwd())
        snpdata = Pheno("../examples/toydata.phe").read()  # Read data from Pheno format
        pstutil.create_directory_if_necessary("tempdir/toydata.5chrom.bed")
        Bed.write(
            "tempdir/toydata.5chrom.bed", snpdata, count_A1=False
        )  # Write data in Bed format

    import doctest

    doctest.testmod(optionflags=doctest.ELLIPSIS)
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
