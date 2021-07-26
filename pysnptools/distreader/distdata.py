import numpy as np
from itertools import *
import pandas as pd
import logging
from pysnptools.distreader import DistReader
from pysnptools.standardizer import Unit
from pysnptools.standardizer import Identity
from pysnptools.pstreader import PstData
import warnings
import time

class DistData(PstData,DistReader):
    """
    A :class:`.DistReader` for holding SNP distributions (or similar values) in-memory, along with related *iid*, *sid*, and *pos* information.
    It is usually created by calling the :meth:`.DistReader.read` method on another :class:`.DistReader`, for example, :class:`.Bgen`.
    It can also be constructed directly.

    See :class:`.DistReader` for details and examples.

    **Constructor:**
        :Parameters: * **iid** (an array of string pair) -- The :attr:`DistReader.iid` information.
                     * **sid** (an array of strings) -- The :attr:`DistReader.sid` information.
                     * **val** (a 3-D array of floats) -- The SNP value distribution
                     * **pos** (optional, an array of strings) -- The :attr:`DistReader.pos` information
                     * **name** (optional, string) -- Information to be display about the origin of this data
                     * **copyinputs_function** (optional, function) -- *Used internally by optional clustering code*

        :Example:

        >>> from pysnptools.distreader import DistData
        >>> distdata = DistData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],
        ...                     val=[[[.5,.5,0],[0,0,1],[.5,.5,0]],
        ...                          [[0,1.,0],[0,.75,.25],[.5,.5,0]]])
        >>> print((distdata.val[1,1], distdata.iid_count, distdata.sid_count))
        (array([0.  , 0.75, 0.25]), 2, 3)

    **Equality:**

        Two DistData objects are equal if their four arrays (:attr:`DistData.val`, :attr:`DistReader.iid`, :attr:`DistReader.sid`, and :attr:`DistReader.pos`)
        are 'array_equal'. (Their 'name' does not need to be the same).
        If either :attr:`DistData.val` contains NaN, the objects will not be equal. However, :meth:`.DistData.allclose` can be used to treat NaN as
        regular values.

        :Example:

        >>> import numpy as np
        >>> from pysnptools.distreader import DistData
        >>> snpdata1 = DistData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'], 
        ...                     val=[[[.5,.5,0],[0,0,1],[.5,.5,0]],
        ...                          [[0,1.,0],[0,.75,.25],[.5,.5,0]]],
        ...                     pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> snpdata2 = DistData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],
        ...                     val=[[[.5,.5,0],[0,0,1],[.5,.5,0]],
        ...                          [[0,1.,0],[0,.75,.25],[.5,.5,0]]],
        ...                     pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> print(snpdata1 == snpdata2) #True, because all the arrays have the same values.
        True
        >>> print(snpdata1.val is snpdata2.val) #False, because the two arrays have different memory.
        False
        >>> snpdata3 = DistData(iid=[['a','0'],['b','0']], sid=['snp334','snp349','snp921'],
        ...                     val=[[[.5,.5,0],[0,0,1],[.5,.5,0]],
        ...                          [[0,1.,0],[0,.75,.25],[.5,.5,0]]],
        ...                     pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> snpdata4 = DistData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],
        ...                     val=[[[.5,.5,0],[0,0,1],[.5,.5,0]],
        ...                          [[0,1.,0],[0,.75,.25],[.5,.5,0]]],
        ...                     pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> print(snpdata3 == snpdata4) #False, because the iids are different.
        False
        >>> snpdata5 = DistData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],
        ...                     val=[[[.5,.5,0],[0,0,1],[.5,.5,0]],
        ...                          [[0,1.,0],[0,.75,.25],[np.nan,np.nan,np.nan]]],
        ...                     pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> snpdata6 = DistData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],
        ...                     val=[[[.5,.5,0],[0,0,1],[.5,.5,0]],
        ...                          [[0,1.,0],[0,.75,.25],[np.nan,np.nan,np.nan]]],
        ...                     pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> print(snpdata5 == snpdata6) #False, because the val's contain NaN
        False
        >>> print(snpdata5.allclose(snpdata6)) #True, if we consider the NaN as regular values, all the arrays have the same values.
        True

    **Methods beyond** :class:`.DistReader`

    """

    def __init__(self, iid, sid, val, pos=None, name=None, copyinputs_function=None):

        #We don't have a 'super(DistData, self).__init__()' here because DistData takes full responsibility for initializing both its superclasses

        self._val = None

        self._row = PstData._fixup_input(iid,empty_creator=lambda ignore:np.empty([0,2],dtype='str'),dtype='str')
        self._col = PstData._fixup_input(sid,empty_creator=lambda ignore:np.empty([0],dtype='str'),dtype='str')
        self._row_property = PstData._fixup_input(None,count=len(self._row),empty_creator=lambda count:np.empty([count,0],dtype='str'),dtype='str')
        self._col_property = PstData._fixup_input(pos,count=len(self._col),empty_creator=lambda count:np.full([count, 3], np.nan))
        self._val = PstData._fixup_input_val(val,row_count=len(self._row),col_count=len(self._col),empty_creator=lambda row_count,col_count:np.empty([row_count,col_count,3],dtype=np.float64))#!!!Replace empty with my FillNA method?
        self._assert_iid_sid_pos(check_val=True)
        self._name = name or ""
        self._std_string_list = []

    @property
    def val(self):
        """The 3D NumPy array of floats that represents the distribution of SNP values. You can get or set this property.

        >>> from pysnptools.distreader import Bgen
        >>> from pysnptools.util import example_file # Download and return local file name
        >>> bgen_file = example_file("pysnptools/examples/2500x100.bgen")
        >>> distdata = Bgen(bgen_file)[:5,:].read() #read data for first 5 iids
        >>> print(distdata.val[4,55]) #print one of the SNP values
        [0.23137255 0.65342184 0.11520562]
        """
        return self._val

    @val.setter
    def val(self, new_value):
        self._val = PstData._fixup_input_val(new_value,row_count=len(self._row),col_count=len(self._col),empty_creator=lambda row_count,col_count:np.empty([row_count,col_count,3],dtype=np.float64))#!!!Replace empty with my FillNA method?
        self._assert_iid_sid_pos(check_val=True)

    def allclose(self,value,equal_nan=True):
        '''
        :param value: Other object with which to compare.
        :type value: :class:`DistData`
        :param equal_nan: (Default: True) Tells if NaN in :attr:`DistData.val` should be treated as regular values when testing equality.
        :type equal_nan: bool

        >>> import numpy as np
        >>> snpdata5 = DistData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],
        ...                     val=[[[.5,.5,0],[0,0,1],[.5,.5,0]],
        ...                          [[0,1.,0],[0,.75,.25],[np.nan,np.nan,np.nan]]],
        ...                     pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> snpdata6 = DistData(iid=[['fam0','iid0'],['fam0','iid1']], sid=['snp334','snp349','snp921'],
        ...                     val=[[[.5,.5,0],[0,0,1],[.5,.5,0]],
        ...                          [[0,1.,0],[0,.75,.25],[np.nan,np.nan,np.nan]]],
        ...                     pos=[[0,0,0],[0,0,0],[0,0,0]])
        >>> print(snpdata5.allclose(snpdata6)) #True, if we consider the NaN as regular values, all the arrays have the same values.
        True
        >>> print(snpdata5.allclose(snpdata6,equal_nan=False)) #False, if we consider the NaN as special values, all the arrays are not equal.
        False

        '''
        return PstData.allclose(self,value,equal_nan=equal_nan)


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


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
    # There is also a unit test case in 'pysnptools\test.py' that calls this doc test
