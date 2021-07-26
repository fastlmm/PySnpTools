from pysnptools.pstreader import PstNpz
from pysnptools.distreader import DistReader
import logging
import scipy as np
import warnings

class DistNpz(PstNpz,DistReader):
    '''
    A :class:`.DistReader` for reading \*.dist.npz files from disk.

    See :class:`.DistReader` for general examples of using DistReaders.

    The general NPZ format is described `here <http://docs.scipy.org/doc/numpy/reference/generated/numpy.savez.html>`__. The DistNpz format stores
    val, iid, sid, and pos information in NPZ format.
   
    **Constructor:**
        :Parameters: * **filename** (*string*) -- The DistNpz file to read.

        :Example:

        >>> from pysnptools.distreader import DistNpz
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> npz_file = example_file("pysnptools/examples/toydata10.dist.npz")
        >>> data_on_disk = DistNpz(npz_file)
        >>> print((data_on_disk.iid_count, data_on_disk.sid_count))
        (25, 10)

    **Methods beyond** :class:`.DistReader`

    '''

    def __init__(self, *args, **kwargs):
        super(DistNpz, self).__init__(*args, **kwargs)

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
    def write(filename, distdata):
        """Writes a :class:`DistData` to DistNpz format and returns the :class:`.DistNpz`

        :param filename: the name of the file to create
        :type filename: string
        :param distdata: The in-memory data that should be written to disk.
        :type distdata: :class:`DistData`
        :rtype: :class:`.DistNpz`

        >>> from pysnptools.distreader import DistNpz, DistHdf5
        >>> import pysnptools.util as pstutil
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> hdf5_file = example_file("pysnptools/examples/toydata.iidmajor.dist.hdf5")
        >>> distdata = DistHdf5(hdf5_file)[:,:10].read()     # Read first 10 snps from DistHdf5 format
        >>> pstutil.create_directory_if_necessary("tempdir/toydata10.dist.npz")
        >>> DistNpz.write("tempdir/toydata10.dist.npz",distdata)          # Write data in DistNpz format
        DistNpz('tempdir/toydata10.dist.npz')
        """
        row_ascii = np.array(distdata.row,dtype='S') #!!! would be nice to avoid this copy when not needed.
        col_ascii = np.array(distdata.col,dtype='S') #!!! would be nice to avoid this copy when not needed.
        np.savez(filename, row=row_ascii, col=col_ascii, row_property=distdata.row_property, col_property=distdata.col_property,val=distdata.val)
        logging.debug("Done writing " + filename)
        return DistNpz(filename)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
