import logging
import scipy as np
from pysnptools.snpreader.snpreader import SnpReader
from pysnptools.pstreader import PstHdf5
import warnings

class SnpHdf5(PstHdf5,SnpReader):
    '''
    A :class:`.SnpReader` for reading \*.snp.hdf5 files from disk.

    See :class:`.SnpReader` for general examples of using SnpReaders.

    The general HDF5 format is described in http://www.hdfgroup.org/HDF5/. The SnpHdf5 format stores
    val, iid, sid, and pos information in Hdf5 format.
   
    **Constructor:**
        :Parameters: * **filename** (*string*) -- The SnpHdf5 file to read.

        :Example:

        >>> from pysnptools.snpreader import SnpHdf5
        >>> data_on_disk = SnpHdf5('../examples/toydata.snpmajor.snp.hdf5')
        >>> print data_on_disk.iid_count, data_on_disk.sid_count
        500 10000

    **Methods beyond** :class:`.SnpReader`

    '''
    def __init__(self, *args, **kwargs):
        super(SnpHdf5, self).__init__(*args, **kwargs)

    @staticmethod
    def write(filename, snpdata, hdf5_dtype=None, sid_major=True):
        """Writes a :class:`SnpData` to SnpHdf5 format.

        :param filename: the name of the file to create
        :type filename: string
        :param snpdata: The in-memory data that should be written to disk.
        :type snpdata: :class:`SnpData`
        :param hdf5_dtype: None (use the .val's dtype) or a Hdf5 dtype, e.g. 'f8','f4',etc.
        :type hdf5_dtype: string
        :param col_major: Tells if vals should be stored on disk in sid_major (default) or iid_major format.
        :type col_major: bool

        >>> from pysnptools.snpreader import SnpHdf5, Bed
        >>> import pysnptools.util as pstutil
        >>> snpdata = Bed('../examples/toydata.bed',count_A1=False)[:,:10].read()     # Read first 10 snps from Bed format
        >>> pstutil.create_directory_if_necessary("tempdir/toydata10.snp.hdf5")
        >>> SnpHdf5.write("tempdir/toydata10.snp.hdf5",snpdata)        # Write data in SnpHdf5 format
        """
        PstHdf5.write(filename,snpdata,hdf5_dtype=hdf5_dtype,col_major=sid_major)

class Hdf5(SnpHdf5):
    #!! warnings.warn("class 'Hdf5' is deprecated. Use the standard class 'SnpHdf5' instead", DeprecationWarning)
    def __init__(self, *args, **kwargs):
        super(Hdf5, self).__init__(*args, **kwargs)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
