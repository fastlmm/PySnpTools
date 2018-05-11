from pysnptools.pstreader import PstNpz
from pysnptools.snpreader.snpreader import SnpReader
import logging
import scipy as np
import warnings

class SnpNpz(PstNpz,SnpReader):
    '''
    A :class:`.SnpReader` for reading \*.snp.npz files from disk.

    See :class:`.SnpReader` for general examples of using SnpReaders.

    The general NPZ format is described in http://docs.scipy.org/doc/numpy/reference/generated/numpy.savez.html. The SnpNpz format stores
    val, iid, sid, and pos information in NPZ format.
   
    **Constructor:**
        :Parameters: * **filename** (*string*) -- The SnpNpz file to read.

        :Example:

        >>> from pysnptools.snpreader import SnpNpz
        >>> data_on_disk = SnpNpz('../examples/toydata10.snp.npz')
        >>> print data_on_disk.iid_count, data_on_disk.sid_count
        500 10

    **Methods beyond** :class:`.SnpReader`

    '''

    def __init__(self, *args, **kwargs):
        super(SnpNpz, self).__init__(*args, **kwargs)

    @staticmethod
    def write(filename, snpdata):
        """Writes a :class:`SnpData` to SnpNpz format.

        :param filename: the name of the file to create
        :type filename: string
        :param snpdata: The in-memory data that should be written to disk.
        :type snpdata: :class:`SnpData`

        >>> from pysnptools.snpreader import SnpNpz, Bed
        >>> import pysnptools.util as pstutil
        >>> snpdata = Bed('../examples/toydata.bed',count_A1=False)[:,:10].read()     # Read first 10 snps from Bed format
        >>> pstutil.create_directory_if_necessary("tempdir/toydata10.snp.npz")
        >>> SnpNpz.write("tempdir/toydata10.snp.npz",snpdata)          # Write data in SnpNpz format
        """
        PstNpz.write(filename,snpdata)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
