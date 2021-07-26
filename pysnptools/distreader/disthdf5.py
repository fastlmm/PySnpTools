import logging
import numpy as np
from pysnptools.distreader import DistReader
from pysnptools.pstreader import PstHdf5
import warnings

class DistHdf5(PstHdf5,DistReader):
    '''
    A :class:`.DistReader` for reading \*.dist.hdf5 files from disk.

    See :class:`.DistReader` for general examples of using DistReaders.

    The general HDF5 format is described `here <http://www.hdfgroup.org/HDF5/>`__. The DistHdf5 format stores
    val, iid, sid, and pos information in Hdf5 format.
   
    **Constructor:**
        :Parameters: * **filename** (*string*) -- The DistHdf5 file to read.

        :Example:

        >>> from pysnptools.distreader import DistHdf5
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> hdf5_file = example_file("pysnptools/examples/toydata.snpmajor.dist.hdf5")
        >>> data_on_disk = DistHdf5(hdf5_file)
        >>> print((data_on_disk.iid_count, data_on_disk.sid_count))
        (25, 10000)

    **Methods beyond** :class:`.DistReader`

    '''
    def __init__(self, *args, **kwargs):
        super(DistHdf5, self).__init__(*args, **kwargs)


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
    def write(filename, distdata, hdf5_dtype=None, sid_major=True):
        """Writes a :class:`DistData` to DistHdf5 format and return a the :class:`.DistHdf5`.

        :param filename: the name of the file to create
        :type filename: string
        :param distdata: The in-memory data that should be written to disk.
        :type distdata: :class:`DistData`
        :param hdf5_dtype: None (use the .val's dtype) or a Hdf5 dtype, e.g. 'f8','f4',etc.
        :type hdf5_dtype: string
        :param sid_major: Tells if vals should be stored on disk in sid_major (default) or iid_major format.
        :type col_major: bool
        :rtype: :class:`.DistHdf5`

        >>> from pysnptools.distreader import DistHdf5, Bgen
        >>> import pysnptools.util as pstutil
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bgen_file = example_file("pysnptools/examples/2500x100.bgen")
        >>> distdata = Bgen(bgen_file)[:,:10].read()     # Read first 10 snps from BGEN format
        >>> pstutil.create_directory_if_necessary("tempdir/toydata10.dist.hdf5")
        >>> DistHdf5.write("tempdir/toydata10.dist.hdf5",distdata)        # Write data in DistHdf5 format
        DistHdf5('tempdir/toydata10.dist.hdf5')
        """
        PstHdf5.write(filename,distdata,hdf5_dtype=hdf5_dtype,col_major=sid_major)
        return DistHdf5(filename)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
