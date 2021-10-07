import numpy as np
import logging
from pysnptools.eigenreader import EigenReader
from pysnptools.snpreader import SnpData
from pysnptools.pstreader import PstData

# cmk


class EigenData(PstData, EigenReader):
    # !!!cmk should it default to eid_0 not eid0?
    # !!!cmk tell what eid defaults to
    """cmk iid to row
    A :class:`.EigenReader` for holding SNP distributions (or similar values) in-memory, along with related *iid*, *eid*, and *values* information.
    It is usually created by calling the :meth:`.EigenReader.read` method on another :class:`.EigenReader`, for example, :class:`.Bgen`.
    It can also be constructed directly.

    See :class:`.EigenReader` for details and examples.

    **Constructor:**
        :Parameters: * **iid** (an array of string pair) -- The :attr:`EigenReader.iid` information.
                     * **values** (optional, an array of strings) -- The :attr:`EigenReader.values` information
                     * **vectors** (a 2-D array of floats) -- The column eigenvectors
                     * **eid** (optional, an array of string pairs) -- The :attr:`EigenReader.eid` information.
                     * **name** (optional, string) -- Information to be display about the origin of this data
                     * **copyinputs_function** (optional, function) -- *Used internally by optional clustering code*

        :Example:

        >>> from pysnptools.eigenreader import EigenData
        >>> eigendata = EigenData(iid=[['fam0','iid0'],['fam0','iid1']],
        ...                     values=[0.88196601, 3.11803399],
        ...                     vectors=[[ 0.22975292, -0.97324899],
        ...                              [-0.97324899, -0.22975292]])
        >>> print((eigendata.values[0],eigendata.vectors[:,0], eigendata.iid_count, eigendata.eid_count))
        (0.88196601, array([ 0.22975292, -0.97324899]), 2, 2)

    **Equality:**

        Two EigenData objects are equal if their four arrays (:attr:`EigenReader.values`, :attr:`EigenData.vectors`, :attr:`EigenReader.iid`, and :attr:`EigenReader.eid`)
        are 'array_equal'. (Their 'name' does not need to be the same).
        If either :attr:`EigenData.vectors` contains NaN, the objects will not be equal. However, :meth:`.EigenData.allclose` can be used to treat NaN as
        regular values.

        :Example:

        >>> import numpy as np
        >>> from pysnptools.eigenreader import EigenData
        >>> eigendata1 = EigenData(iid=[['fam0','iid0'],['fam0','iid1']],
        ...                     values=[0.88196601, 3.11803399],
        ...                     vectors=[[ 0.22975292, -0.97324899],
        ...                              [-0.97324899, -0.22975292]])
        >>> eigendata2 = EigenData(iid=[['fam0','iid0'],['fam0','iid1']],
        ...                     values=[0.88196601, 3.11803399],
        ...                     vectors=[[ 0.22975292, -0.97324899],
        ...                              [-0.97324899, -0.22975292]])
        >>> print(eigendata1 == eigendata2) #True, because all the arrays have the same values.
        True
        >>> print(eigendata1.val is eigendata2.val) #False, because the two arrays have different memory.
        False
        >>> eigendata3 = EigenData(iid=[['a','b'],['c','d']],
        ...                     values=[0.88196601, 3.11803399],
        ...                     vectors=[[ 0.22975292, -0.97324899],
        ...                              [-0.97324899, -0.22975292]])
        >>> eigendata4 = EigenData(iid=[['fam0','iid0'],['fam0','iid1']],
        ...                     values=[0.88196601, 3.11803399],
        ...                     vectors=[[ 0.22975292, -0.97324899],
        ...                              [-0.97324899, -0.22975292]])
        >>> print(eigendata3 == eigendata4) #False, because the iids are different.
        False

    **Methods beyond** :class:`.EigenReader`

    """

    def __init__(
        self, row, values, vectors, eid=None, name=None, copyinputs_function=None
    ):

        # We don't have a 'super(EigenData, self).__init__()' here because EigenData takes full responsibility for initializing both its superclasses

        ##self._val = None

        self._row = PstData._fixup_input(
            row, empty_creator=lambda ignore: np.empty([0, 2], dtype="str"), dtype="str"
        )
        self._row_property = PstData._fixup_input(
            None,
            count=len(self._row),
            empty_creator=lambda count: np.empty([count, 0], dtype="str"),
            dtype="str",
        )

        self._col_property = PstData._fixup_input(
            values, dtype=np.float64
        )  # !!!cmk need to raise error of values is None
        self._col = PstData._fixup_input(
            eid,
            count=len(self._col_property),
            dtype="str",
            empty_creator=lambda count: np.array(
                [["", f"eid{eid_index}"] for eid_index in range(count)]
            ),
        )

        self._val = PstData._fixup_input_val(
            vectors,
            row_count=len(self._row),
            col_count=len(self._col),
            empty_creator=lambda row_count, col_count: np.empty(
                [row_count, col_count], dtype=np.float64
            ),
        )  #!!!Replace empty with my FillNA method?

        self._assert_row_eid_values(check_vectors=True)
        self._name = name or ""
        self._std_string_list = []

    @staticmethod
    def from_aka(aka, keep_above=np.NINF):
        # !!!cmk check that square aKa not just aKb???
        val = aka.val
        if len(val.shape) == 3:
            assert val.shape[2] == 1, "Expect to run on just one phenotype"
            val = np.squeeze(val, -1)
        w, v = np.linalg.eigh(val)  # !!! cmk do SVD sometimes?
        eigen = EigenData(values=w, vectors=v, row=aka.row)
        if keep_above > np.NINF:
            eigen = eigen[:, eigen.values > keep_above].read(view_ok=True)
        return eigen

    @property
    def vectors(self):
        """The 2D NumPy array of floats that represents the column eigen vectors. You can get or set this property.

        >>> from pysnptools.eigenreader import EigenData # !!!cmk later change this example to read from disk
        >>> eigendata = EigenData(row=[['fam0','iid0'],['fam0','iid1']],
        ...                     values=[0.88196601, 3.11803399],
        ...                     vectors=[[ 0.22975292, -0.97324899],
        ...                              [-0.97324899, -0.22975292]])
        >>> print(eigendata.vectors[:,0]) # print the first eigenvector
        [ 0.22975292 -0.97324899]
        """
        return self._val

    @vectors.setter
    def vectors(self, new_value):
        self._val = PstData._fixup_input_val(
            new_value,
            row_count=len(self._row),
            col_count=len(self._col),
            empty_creator=lambda row_count, col_count: np.empty(
                [row_count, col_count, 2], dtype=np.float64
            ),
        )  #!!!Replace empty with my FillNA method?
        self._assert_row_eid_values(check_val=True)

    def allclose(self, value, equal_nan=True):
        """
        :param value: Other object with which to compare.
        :type value: :class:`EigenData`
        :param equal_nan: (Default: True) Tells if NaN in :attr:`EigenData.vectors` should be treated as regular values when testing equality.
        :type equal_nan: bool

        >>> import numpy as np
        >>> eigendata5 = EigenData(row=[['fam0','iid0'],['fam0','iid1']],
        ...                     values=[0.88196601, 3.11803399],
        ...                     vectors=[[ np.nan, np.nan],
        ...                              [np.nan, np.nan]])
        >>> eigendata6 = EigenData(row=[['fam0','iid0'],['fam0','iid1']],
        ...                     values=[0.88196601, 3.11803399],
        ...                     vectors=[[ np.nan, np.nan],
        ...                              [np.nan, np.nan]])
        >>> print(eigendata5.allclose(eigendata6)) #True, if we consider the NaN as regular values, all the arrays have the same values.
        True
        >>> print(eigendata5.allclose(eigendata6,equal_nan=False)) #False, if we consider the NaN as special values, all the arrays are not equal.
        False

        """
        return PstData.allclose(self, value, equal_nan=equal_nan)

    def logdet(self, delta=None):
        # !!!cmk could have path for delta=0
        # "reshape" lets it broadcast
        if delta is None:
            Sd = self.values
        else:
            Sd = self.values + delta
        logdet = np.log(Sd).sum()
        if self.is_low_rank:  # !!!cmk test this
            logdet += (self.row_count - self.eid_count) * np.log(delta)
        return logdet, Sd

    #!!!cmk document
    #!!!cmk how to understand the low rank bit?
    def rotate_and_scale(self, pstdata, ignore_low_rank=False):
        rotation = self.rotate(pstdata, ignore_low_rank=ignore_low_rank)
        rotation.val[:, :] = rotation.val / self.values.reshape(-1, 1)
        return rotation

    def rotate(self, pstdata, ignore_low_rank=False):
        val = pstdata.val
        if len(val.shape) == 3:
            val = np.squeeze(val, -1)
        rotated_val = np.einsum("ae,ab->eb", self.vectors, val)
        rotated_pstdata = PstData(
            row=self.col, col=pstdata.col, val=rotated_val, name=f"rotated({pstdata})"
        )

        rotation = Rotation(rotated_pstdata, double=None)

        if self.is_low_rank and not ignore_low_rank:
            rotated_back_pstdata = self.rotate_back(rotation, check_low_rank=False)
            double_pstdata = rotated_back_pstdata.clone(
                val=pstdata.val - rotated_back_pstdata.val, name=f"double({pstdata})"
            )
            rotation = Rotation(rotated_pstdata, double=double_pstdata)

        ## !!!cmk make a test of this kludge
        # if not np.allclose(val, self.rotate_back(rotation).val, rtol=0, atol=1e-9):
        #    pstdata2 = self.rotate_back(rotation)
        #    assert False

        return rotation

    #!!!cmk how to understand the low rank bit?
    def rotate_back(self, rotation, check_low_rank=True):
        if check_low_rank:
            assert (
                rotation.double is not None
            ) == self.is_low_rank, "low rank eigens expect a non-empty rotation.double"

        val = rotation.val
        if len(val.shape) == 3:
            val = np.squeeze(val, -1)
        rotated_back_val = np.einsum("ae,eb->ab", self.vectors, val)

        if rotation.double is not None:
            rotated_back_val += rotation.double.val

        rotated_back_pstdata = PstData(
            row=self.row,
            col=rotation.col,
            val=rotated_back_val,
            name=f"rotated_back({rotation})",
        )
        return rotated_back_pstdata

    def __repr__(self):
        if self._name == "":
            if len(self._std_string_list) > 0:
                s = "{0}({1})".format(
                    self.__class__.__name__, ",".join(self._std_string_list)
                )
            else:
                s = "{0}()".format(self.__class__.__name__)
        else:
            if len(self._std_string_list) > 0:
                s = "{0}({1},{2})".format(
                    self.__class__.__name__, self._name, ",".join(self._std_string_list)
                )
            else:
                s = "{0}({1})".format(self.__class__.__name__, self._name)
        return s

    #!!!cmk document


class Rotation:
    #!!!cmk kludge diagonal_name = np.array(["diagonal"])  #!!!cmk similar code

    def __init__(self, rotated, double):
        self.rotated = rotated
        self.double = double

    def __getitem__(self, index):
        rotated = self.rotated[:, index : index + 1].read(view_ok=True)

        if self.double is not None:
            double = self.double[:, index : index + 1].read(view_ok=True)
        else:
            double = None
        return Rotation(rotated, double)

    @property
    def row_count(self):
        return self.rotated.row_count

    @property
    def col_count(self):
        return self.rotated.col_count

    @property
    def row(self):
        return self.rotated.row

    @property
    def col(self):
        return self.rotated.col

    def ein(self, s):
        assert len(s) == 1, "length of string must be 1"
        assert s != "d", "string can not be 'd'"  #!!!cmk remove kludge
        return s, np.s_[:]

    # @staticmethod#!!!cmk kludge
    # def ein_cat(*args):
    #    result = ""
    #    for i in range(len(args) - 1, -1, -1):
    #        arg = args[i]
    #        assert len(arg) == 1, "Expect inputs to be one letter"
    #        if arg in result:
    #            assert (
    #                arg == "d"
    #            ), "if a letter appears twice it should be 'd' for diagonal"
    #        else:
    #            result = arg + result
    #    return result

    @staticmethod
    def ein_d(a, b):
        result = "d"
        if b != "d":
            result = b + result
        if a != "d":
            result = a + result
        return result

    @property
    def val(self):
        return self.rotated.val

    @property
    def shape(self):
        return self.rotated.shape


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest

    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
