from pysnptools.kernelreader import KernelReader
from pysnptools.pstreader._subset import _PstSubset

class _KernelSubset(KernelReader,_PstSubset):
    def __init__(self, *args, **kwargs):
        super(_KernelSubset, self).__init__(*args, **kwargs)


    def _read_array(self, row_index, col_index, order, dtype,
                    force_python_only, view_ok, num_threads, xp):
        self._run_once()
        row_index = self.compose_indexer_with_indexer(
            self._internal.row_count, self._row_indexer, self.row_count, row_index)
        col_index = self.compose_indexer_with_indexer(
            self._internal.col_count, self._col_indexer, self.col_count, col_index)
        return self._internal._read_array(row_index, col_index, order, dtype,
                                          force_python_only, view_ok, num_threads, xp)
