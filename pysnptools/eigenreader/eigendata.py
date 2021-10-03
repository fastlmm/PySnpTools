import numpy as np
import logging
from pysnptools.eigenreader import EigenReader
from pysnptools.snpreader import SnpData
from pysnptools.pstreader import PstData
#cmk

class EigenData(PstData,EigenReader):
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

    def __init__(self, row, values, vectors, eid=None, name=None, copyinputs_function=None):

        #We don't have a 'super(EigenData, self).__init__()' here because EigenData takes full responsibility for initializing both its superclasses

        ##self._val = None

        self._row = PstData._fixup_input(row,empty_creator=lambda ignore:np.empty([0,2],dtype='str'),dtype='str')
        self._row_property = PstData._fixup_input(None,count=len(self._row),empty_creator=lambda count:np.empty([count,0],dtype='str'),dtype='str')

        self._col_property = PstData._fixup_input(values,dtype=np.float64) # !!!cmk need to raise error of values is None
        self._col = PstData._fixup_input(eid,count=len(self._col_property),dtype='str',empty_creator=lambda count:np.array([["",f"eid{eid_index}"] for eid_index in range(count)]))

        self._val = PstData._fixup_input_val(vectors,row_count=len(self._row),col_count=len(self._col),
                                             empty_creator=lambda row_count,col_count:np.empty([row_count,col_count],
                                             dtype=np.float64))#!!!Replace empty with my FillNA method?

        self._assert_row_eid_values(check_vectors=True)
        self._name = name or ""
        self._std_string_list = []

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
        self._val = PstData._fixup_input_val(new_value,row_count=len(self._row),col_count=len(self._col),empty_creator=lambda row_count,col_count:np.empty([row_count,col_count,2],dtype=np.float64))#!!!Replace empty with my FillNA method?
        self._assert_row_eid_values(check_val=True)

    def allclose(self,value,equal_nan=True):
        '''
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

        '''
        return PstData.allclose(self,value,equal_nan=equal_nan)

    def logdet(self, delta=None):
        # !!!cmk could have path for delta=0
        # "reshape" lets it broadcast
        if delta is None:
            Sd = self.values.reshape(-1,1)
        else:
            Sd = self.values.reshape(-1, 1) + delta
        logdet = np.log(Sd).sum(axis=0).reshape(1,-1)
        if self.is_low_rank:  # !!!cmk test this
            logdet += (self.row_count - self.eid_count) * np.log(delta)
        return logdet, Sd # !!!cmk shape 1,3 50,3


    #!!!cmk document
    #!!!cmk should this take snpdata instead?
    #!!!cmk how to understand the low rank bit?
    def rotate(self, pstdata, is_diagonal=False):
        val = pstdata.val
        if len(val.shape)==3: #!!!cmk ugly
            val = np.squeeze(val,-1)
        rotated_val = self.vectors.T.dot(val)
        #!!!cmk make row calc faster
        rotated_pstdata = PstData(row=self.col, col=pstdata.col, val=rotated_val, name=f"rotated({pstdata})")

        if self.is_low_rank:
            double_val = pstdata.val - self.vectors.dot(rotated_val)
            double_pstdata = PstData(row=pstdata.row, col=pstdata.col, val=double_val, name=f"double({pstdata})")
        else:
            double_pstdata = None
        return Rotation(rotated_pstdata, double_pstdata, is_diagonal=is_diagonal)

    #!!!cmk should this take snpdata instead?
    #!!!cmk how to understand the low rank bit?
    def t_rotate(self, pstdata, is_diagonal=False): #!!!cmk not pstdata
        t_rotated_val = self.vectors.dot(pstdata.val)
        #!!!cmk make row calc faster
        rotated_pstdata = PstData(row=self.row, col=pstdata.col, val=t_rotated_val, name=f"t_rotated({pstdata})")

        assert(not self.is_low_rank) # !!!cmk need code or better error message
        double_pstdata = None
        return Rotation(rotated_pstdata, double_pstdata, is_diagonal=is_diagonal)


    def __repr__(self):
        if self._name == "":
            if len(self._std_string_list) > 0:
                s = "{0}({1})".format(self.__class__.__name__,",".join(self._std_string_list))
            else:
                s = "{0}()".format(self.__class__.__name__)
        else:
            if len(self._std_string_list) > 0:
                s = "{0}({1},{2})".format(self.__class__.__name__,self._name,",".join(self._std_string_list))
            else:
                s = "{0}({1})".format(self.__class__.__name__,self._name)
        return s

    #!!!cmk document
class Rotation:
    diagonal_name = np.array(["diagonal"])  #!!!cmk similar code

    def __init__(self, rotated, double, is_diagonal=False):
        self.rotated = rotated
        self.double = double
        self.is_diagonal = is_diagonal


    def __getitem__(self, index):
        rotated = self.rotated[:,index:index+1].read(view_ok=True)
        if self.is_diagonal:
            rotated = rotated.clone(col=self.diagonal_name)

        if self.double is not None:
            double = self.double[:,index:index+1].read(view_ok=True)
        else:
            double = None
        return Rotation(rotated,double,is_diagonal=self.is_diagonal)

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
