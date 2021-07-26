from pysnptools.pstreader import PstNpz
from pysnptools.snpreader import SnpReader
import logging
import scipy as np
import warnings

class SnpNpz(PstNpz,SnpReader):
    '''
    A :class:`.SnpReader` for reading \*.snp.npz files from disk.

    See :class:`.SnpReader` for general examples of using SnpReaders.

    The general NPZ format is described `here <http://docs.scipy.org/doc/numpy/reference/generated/numpy.savez.html>`__. The SnpNpz format stores
    val, iid, sid, and pos information in NPZ format.
   
    **Constructor:**
        :Parameters: * **filename** (*string*) -- The SnpNpz file to read.

        :Example:

        >>> from pysnptools.snpreader import SnpNpz
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> npz_file = example_file('pysnptools/examples/toydata10.snp.npz')
        >>> data_on_disk = SnpNpz(npz_file)
        >>> print((data_on_disk.iid_count, data_on_disk.sid_count))
        (500, 10)

    **Methods beyond** :class:`.SnpReader`

    '''

    def __init__(self, *args, **kwargs):
        super(SnpNpz, self).__init__(*args, **kwargs)

    @property
    def row(self):
        self._run_once()
        if self._row.dtype.type is not np.str_:
            self._row = np.array(self._row,dtype='str')
        return self._row

    @property
    def col(self):
        self._run_once()
        if self._col.dtype.type is not np.str_:
            self._col = np.array(self._col,dtype='str')
        return self._col


    @staticmethod
    def write(filename, snpdata):
        """Writes a :class:`SnpData` to SnpNpz format and returns the :class:`.SnpNpz`

        :param filename: the name of the file to create
        :type filename: string
        :param snpdata: The in-memory data that should be written to disk.
        :type snpdata: :class:`SnpData`
        :rtype: :class:`.SnpNpz`

        >>> from pysnptools.snpreader import SnpNpz, Bed
        >>> import pysnptools.util as pstutil
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bed_file = example_file("pysnptools/examples/toydata.5chrom.*","*.bed")
        >>> snpdata = Bed(bed_file,count_A1=False)[:,:10].read()     # Read first 10 snps from Bed format
        >>> pstutil.create_directory_if_necessary("tempdir/toydata10.snp.npz")
        >>> SnpNpz.write("tempdir/toydata10.snp.npz",snpdata)          # Write data in SnpNpz format
        SnpNpz('tempdir/toydata10.snp.npz')
        """
        row_ascii = np.array(snpdata.row,dtype='S') #!!! would be nice to avoid this copy when not needed.
        col_ascii = np.array(snpdata.col,dtype='S') #!!! would be nice to avoid this copy when not needed.
        np.savez(filename, row=row_ascii, col=col_ascii, row_property=snpdata.row_property, col_property=snpdata.col_property,val=snpdata.val)
        logging.debug("Done writing " + filename)
        return SnpNpz(filename)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
