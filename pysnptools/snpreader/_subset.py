from pysnptools.snpreader import SnpReader
from pysnptools.pstreader._subset import _PstSubset

class _SnpSubset(_PstSubset,SnpReader):
    def __init__(self, *args, **kwargs):
        super(_SnpSubset, self).__init__(*args, **kwargs)

    def _read_array(self, row_index, col_index, order, dtype,
                    force_python_only, view_ok, num_threads, xp):
        self._run_once()
        row_index = self.compose_indexer_with_indexer(
            self._internal.row_count, self._row_indexer, self.row_count, row_index)
        col_index = self.compose_indexer_with_indexer(
            self._internal.col_count, self._col_indexer, self.col_count, col_index)
        return self._internal._read_array(row_index, col_index, order, dtype,
                                          force_python_only, view_ok, num_threads, xp)

    def _read_with_standardizing(self, standardizer, order="F", dtype=None,
                                 xp=None, return_trained=False,
                                 force_python_only=False, num_threads=None,
                                 _row_index=None, _col_index=None):
        self._run_once()
        row_index = self.compose_indexer_with_indexer(
            self._internal.row_count, self._row_indexer, self.row_count, _row_index)
        col_index = self.compose_indexer_with_indexer(
            self._internal.col_count, self._col_indexer, self.col_count, _col_index)
        result = self._internal._read_with_standardizing(
            standardizer, order=order, dtype=dtype, xp=xp,
            return_trained=return_trained, force_python_only=force_python_only,
            num_threads=num_threads, _row_index=row_index, _col_index=col_index)
        return result
