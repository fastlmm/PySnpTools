import numpy as np
from pysnptools.pstreader import PstNpz
from pysnptools.eigenreader import EigenReader
import logging
import warnings

class EigenNpz(PstNpz,EigenReader):
    '''
    A :class:`.EigenReader` for reading \*.eigen.npz files from disk.

    See :class:`.EigenReader` for general examples of using EigenReaders.

    The general NPZ format is described `here <http://docs.scipy.org/doc/numpy/reference/generated/numpy.savez.html>`__. The EigenNpz format stores
    vectors, iid, eid, and values information in NPZ format.
   
    **Constructor:**
        :Parameters: * **filename** (*string*) -- The EigenNpz file to read.

        :Example:

        #!!!cmk
        #>>> from pysnptools.eigenreader import EigenNpz
        #>>> from pysnptools.util import example_file # Download and return local file name
        #>>> npz_file = example_file("pysnptools/examples/toydata10.eigen.npz")
        #>>> data_on_disk = EigenNpz(npz_file)
        #>>> print((data_on_disk.iid_count, data_on_disk.sid_count))
        #(25, 10)

    **Methods beyond** :class:`.EigenReader`

    '''

    def __init__(self, *args, **kwargs):
        super(EigenNpz, self).__init__(*args, **kwargs)

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
    def write(filename, eigendata):
        """Writes a :class:`EigenData` to EigenNpz format and returns the :class:`.EigenNpz`

        :param filename: the name of the file to create
        :type filename: string
        :param eigendata: The in-memory data that should be written to disk.
        :type eigendata: :class:`EigenData`
        :rtype: :class:`.EigenNpz`

        >>> from pysnptools.eigenreader import EigenNpz, EigenData
        >>> from pysnptools.kernelreader import KernelNpz
        >>> import pysnptools.util as pstutil
        >>> from pysnptools.util import example_file #  Download and return local file name
        >>> kernel_npz_file = example_file('pysnptools/examples/toydata.kernel.npz')
        >>> kernel_data = KernelNpz(kernel_npz_file).read()
        >>> values, vectors = np.linalg.eigh(kernel_data.read().val)
        >>> eigendata = EigenData(values=values,vectors=vectors,iid=kernel_data.iid)
        >>> pstutil.create_directory_if_necessary("tempdir/toydata.eigen.npz")
        >>> EigenNpz.write("tempdir/toydata.eigen.npz",eigendata)          #  Write data in EigenNpz format
        EigenNpz('tempdir/toydata.eigen.npz')

        """
        row_ascii = np.array(eigendata.row,dtype='S') #!!! would be nice to avoid this copy when not needed.
        col_ascii = np.array(eigendata.col,dtype='S') #!!! would be nice to avoid this copy when not needed.
        np.savez(filename, row=row_ascii, col=col_ascii, row_property=eigendata.row_property, col_property=eigendata.col_property,val=eigendata.vectors)
        logging.debug("Done writing " + filename)
        return EigenNpz(filename)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    if False: # !!!cmk
        from pysnptools.eigenreader import EigenNpz, EigenData
        from pysnptools.kernelreader import KernelNpz
        import pysnptools.util as pstutil
        from pysnptools.util import example_file # Download and return local file name
        kernel_npz_file = example_file('pysnptools/examples/toydata.kernel.npz')
        kernel_data = KernelNpz(kernel_npz_file).read()
        values, vectors = np.linalg.eigh(kernel_data.read().val)
        eigendata = EigenData(values=values,vectors=vectors,iid=kernel_data.iid)
        pstutil.create_directory_if_necessary("tempdir/toydata.eigen.npz")
        EigenNpz.write("tempdir/toydata.eigen.npz",eigendata)          # Write data in EigenNpz format

    import doctest
    doctest.testmod()
