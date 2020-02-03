import numpy as np
from itertools import chain,islice
from  more_itertools import unique_everseen
import unittest
import logging
import six #!!!cmk test Python2
import pysnptools.util as pstutil

class Pairs(object):
    """!!!cmkdescription of class"""

    def __init__(self, list0, list1, include_singles, duplicates_ok=True): #!!!cmk rename list to iterable?
        #!!!cmk add _ to start of some of these and some of the methods
        self.include_singles = include_singles
        self.duplicates_ok = duplicates_ok

        self.list0 = list(unique_everseen(list0))
        self.list1 = list(unique_everseen(list1))
        if not duplicates_ok:
            assert len(list0)==len(self.list0) and len(list1)==len(self.list1), "Expect lists to have no duplicates"
        self.dict0 = {k: v for v, k in enumerate(self.list0)} #!!!cmk do these two need to be remembered?
        self.dict1 = {k: v for v, k in enumerate(self.list1)}

        self.common = sorted(six.viewkeys(self.dict0) & self.list1, key=lambda k:self.dict0[k])
        self.only0 = sorted(six.viewkeys(self.dict0)-self.common,key=lambda k:self.dict0[k])
        self.only1 = sorted(six.viewkeys(self.dict1)-self.common,key=lambda k:self.dict1[k])
    
        #How many pairs?
        self.count0 = len(self.only0)*len(self.list1)
        self.count1 = self.count1_fun(len(self.common),len(self.common),len(self.only1),include_singles)
        self.count = self.count0+self.count1

    def pair_sequence(self, start=0,stop=None):
        stop = stop or count
        return islice(self._pair_sequence_inner(start=start),stop-start)

    def _pair_sequence_inner(self,start=0):
        only0_start = start // len(self.list1)
        start -= min(only0_start,len(self.only0))*len(self.list1)
        assert 0 <= start, "real assert"
        for v0 in self.only0[only0_start:]:
            if start > len(self.list1):
                assert False, "real assert"
                start -= len(self.list1)
            else:
                for v1 in islice(chain(self.only1,self.common),start,None):
                    yield v0,v1
                start = 0
        
        common_start = self.index_common_fun(start,len(self.common),len(self.only1),include_singles)
        start -= self.count1_fun(common_start,len(self.common),len(self.only1),include_singles)
        assert 0 <= start and start <= len(self.list1), "real assert"
        for index in range(common_start,len(self.common)):
            v0 = self.common[index]
            startx = index if include_singles else index+1
            if start > len(self.list1)-startx:
                assert False, "real assert"
                start -= (len(self.list1)-startx)
            else:
                for v1 in islice(chain(self.only1,self.common[startx:]),start,None):
                    yield v0,v1
                start = 0

    @staticmethod
    def count1_fun(index_common,len_common,len_only1,include_singles):
        a2 = -1
        b2 = len_only1*2+len_common*2+(1 if include_singles else 3)
        count1 = (a2*index_common*index_common + b2*index_common)//2
        return count1

    @staticmethod
    def isqrt(x): #https://code.activestate.com/recipes/577821-integer-square-root-function/
        if x < 0:
            raise ValueError('square root not defined for negative numbers')
        n = int(x)
        if n == 0:
            return 0
        a, b = divmod(n.bit_length(), 2)
        x = 2**(a+b)
        while True:
            y = (x + n//x)//2
            if y >= x:
                return x
            x = y

    @staticmethod
    def index_common_fun(count1,len_common,len_only1,include_singles):
        a2 = -1
        b2 = len_only1*2+len_common*2+(1 if include_singles else 3)
        c = -count1
        index_common = (b2-Pairs.isqrt(b2*b2-8*a2*c))//2
        return index_common

class TestPairs(unittest.TestCase):

    def test_index_and_count1_functions(self):
        for len_commonq in range(0,5):
            for len_only1q in range(0,5):
                for include_singlesq in [True,False]:
                    for index_commonq in range(len_commonq+1):
                        count1q = count1_fun(index_commonq,len_commonq,len_only1q,include_singlesq)
                        index_common2 = index_common_fun(count1q,len_commonq,len_only1q,include_singlesq)
                        assert index_commonq==index_common2



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod() #!!!cmk put this in main testing


    include_singles = True
    duplicates_ok = True

    if True: #!!!cmk make a test like this
        list0 = ['a','z','a','b','y']
        list1 = ['z','y','x','v']
    else: #!!!cmk make a test like this
        size = 500*1000
        seed = 0
        np.random.seed(seed)
        list0 = np.random.randint(size*10,size=size)
        list1 = np.random.randint(size*10,size=size)
        #list1 = list0 #!!!cmk make a test like this

    pairs = Pairs(list0, list1, include_singles, duplicates_ok)

    for goal_see in [1,2,4,8,15,16]:
        for start in [0,pairs.count//5,pairs.count//2,pairs.count-1,pairs.count]:
            print(start,goal_see,list(pairs.pair_sequence(start,start+goal_see))) #!!!cmk instead of "pair_sequence" how about an indexer? (and "pair" is redundant)


    print("done")
