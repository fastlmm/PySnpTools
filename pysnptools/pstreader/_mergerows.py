#!!! Has not been tested

import logging
import os
import numpy as np
from pysnptools.pstreader import PstReader

class _MergeRows(PstReader): #!!!why does this start with _
    def __init__(self,reader_list,cache_file=None,skip_check=False):
        super(_MergeRows, self).__init__()
        assert len(reader_list) > 0, "Expect at least one reader"
        self.skip_check = skip_check

        self.reader_list = reader_list
        self._repr_string = "_MergeRows({0})".format(",".join([str(s) for s in reader_list]))

        if cache_file is not None:
            if not os.path.exists(cache_file):
                self._run_once()
                np.savez(cache_file, _row=self._row, _row_property=self._row_property, _row_count_list=self._row_count_list)
            else:
                with np.load(cache_file,allow_pickle=True) as data:
                    self._row = data['_row']
                    self._row_property = data['_row_property']
                    self._row_count_list = data['_row_count_list']
                    self._has_run_once = True


    def _run_once(self):
        if hasattr(self,'_has_run_once'):
            return
        self._has_run_once = True
        #Check that all iids are distinct and that all sids and pos are the same and in the same order

        row_list = []
        row_property_list = []
        row_set = set()
        col = self.reader_list[0].col
        col_property = self.reader_list[0].col_property
        for reader_index,reader in enumerate(self.reader_list):
            if reader_index % 10 == 0: logging.info("_MergeRows looking at reader #{0}: {1}".format(reader_index,reader))
            if not self.skip_check:
                assert np.array_equal(col,reader.col), "Expect columns to be the same across all files"
                np.testing.assert_equal(col_property,reader.col_property) #"Expect column_property to be the same across all files"
                size_before = len(row_set)
                row_set.update((tuple(item) for item in reader.row))
                assert len(row_set) == size_before + reader.row_count, "Expect rows to be distinct in all files"
            row_list.append(reader.row)
            row_property_list.append(reader.row_property)
        self._row = np.concatenate(row_list)
        self._row_property = np.concatenate(row_property_list)
        self._row_count_list = [len(row) for row in row_list]

    @property
    def row(self):
        self._run_once()
        return self._row

    @property
    def col(self):
        self._run_once()
        return self.reader_list[0].col

    @property
    def row_property(self):
        self._run_once()
        return self._row_property

    @property
    def col_property(self):
        self._run_once()
        return self.reader_list[0].col_property

    def __repr__(self): 
        #Don't need _run_once because based only on initial info
        return self._repr_string

    def copyinputs(self, copier):
        self._run_once()
        for reader in self.reader_list:
            copier.input(reader)

    def _create_reader_and_row_index_list(self,row_index):
        result = []
        start = 0
        for reader_index in range(len(self.reader_list)): #!!!this needs test coverage
            stop = start + self._row_count_list[reader_index]
            is_here = (row_index >= start) * (row_index < stop)
            if any(is_here):
                row_index_rel = row_index[is_here]-start
                result.append((reader_index,is_here,row_index_rel))
            start = stop
        return result

    def _read(self, row_index_or_none, col_index_or_none, order, dtype, force_python_only, view_ok, num_threads):
        #!!!tests to do: no row's
        #!!!tests to do: no col's
        #!!!test to do: from file 1, file2, and then file1 again
        dtype = np.dtype(dtype)
        row_index = row_index_or_none if row_index_or_none is not None else np.arange(self.row_count) #!!!might want to special case reading all
        col_index_or_none_count = self.col_count if col_index_or_none is None else len(col_index_or_none)
        reader_and_row_index_list = self._create_reader_and_row_index_list(row_index)

        if len(reader_and_row_index_list) == 0:
            return self.reader_list[0]._read(row_index,col_index_or_none,order,dtype,force_python_only, view_ok, num_threads)
        elif len(reader_and_row_index_list) == 1:
            reader_index,row_index_in,row_index_rel = reader_and_row_index_list[0]
            reader = self.reader_list[reader_index]
            return reader._read(row_index_rel,col_index_or_none,order,dtype,force_python_only, view_ok, num_threads)
        else:
            logging.info("Starting read from {0} subreaders".format(len(reader_and_row_index_list)))
            if order == 'A' or order is None:
                order = 'F'
            val = np.empty((len(row_index),col_index_or_none_count),dtype=dtype,order=order)
            for reader_index,is_here,row_index_rel in reader_and_row_index_list:
                reader = self.reader_list[reader_index]
                if reader_index % 1 == 0: logging.info("Reading from #{0}: {1}".format(reader_index,reader))
                val[is_here,:] = reader._read(row_index_rel,col_index_or_none,order,dtype,force_python_only, view_ok=True, num_threads=num_threads)
            logging.info("Ended read from {0} subreaders".format(len(reader_and_row_index_list)))
            return val

