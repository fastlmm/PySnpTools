import numpy as np
from itertools import chain,islice
from  more_itertools import unique_everseen
import doctest
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

    def pair_sequence(self, start=0,stop=None,speed='fast'):
        if speed=='fast':
            stop = min(stop,self.count) if stop is not None else self.count
            return islice(self._pair_sequence_inner(start=start),stop-start)
        if speed=='medium':
            stop = stop if stop is not None else self.count
            return islice(self._pair_sequence_inner_medium(start=start),stop-start)
        else:
            assert speed=='slow',"Don't know speed '{0}'".format(speed)
            return islice(self._pair_sequence_inner_slow(),start,stop)

    def _pair_sequence_inner(self,start=0):
        if start >= self.count0:
            start -= self.count0
        else:
            only0_start = start // len(self.list1)
            start -= only0_start*len(self.list1)
            assert 0 <= start and start <= len(self.list1), "real assert"
            for v0 in self.only0[only0_start:]:
                for v1 in islice(chain(self.only1,self.common),start,None):
                    yield v0,v1
                start = 0
        
        common_start = self.row_index_fun(start,len(self.common),len(self.only1),self.include_singles)
        start -= self.count1_fun(common_start,len(self.common),len(self.only1),self.include_singles)
        assert 0 <= start and start <= len(self.list1), "real assert"
        for index in range(common_start,len(self.common)):
            v0 = self.common[index]
            startx = index if self.include_singles else index+1
            assert start < len(self.list1)-startx, "real assert"
            for v1 in islice(chain(self.only1,self.common[startx:]),start,None):
                yield v0,v1
            start = 0

    def _pair_sequence_inner_medium(self,start=0):
        assert 0 <= start, "real assert"
        for v0 in self.only0:
            if start > len(self.list1):
                start -= len(self.list1)
            else:
                for v1 in islice(chain(self.only1,self.common),start,None):
                    yield v0,v1
                start = 0

        assert 0 <= start, "real assert"
        for index in range(len(self.common)):
            v0 = self.common[index]
            startx = index if self.include_singles else index+1
            if start > len(self.list1)-startx:
                start -= (len(self.list1)-startx)
            else:
                for v1 in islice(chain(self.only1,self.common[startx:]),start,None):
                    yield v0,v1
                start = 0

    def _pair_sequence_inner_slow(self):

        for v0 in self.only0:
            for v1 in chain(self.only1,self.common):
                yield v0,v1

        for index in range(len(self.common)):
            v0 = self.common[index]
            startx = index if self.include_singles else index+1
            for v1 in chain(self.only1,self.common[startx:]):
                yield v0,v1


    @staticmethod
    def count1_fun(row_index,len_common,len_only1,include_singles):
        a2 = -1
        b2 = len_only1*2+len_common*2+(1 if include_singles else -1)
        count1 = (a2*row_index*row_index + b2*row_index)//2
        return count1

    #@staticmethod #!!!cmk delete
    #def isqrt(x): #https://code.activestate.com/recipes/577821-integer-square-root-function/
    #    if x < 0:
    #        raise ValueError('square root not defined for negative numbers')
    #    n = int(x)
    #    if n == 0:
    #        return 0
    #    a, b = divmod(n.bit_length(), 2)
    #    x = 2**(a+b)
    #    while True:
    #        y = (x + n//x)//2
    #        if y >= x:
    #            return x
    #        x = y

    @staticmethod
    def row_index_fun(count1,len_common,len_only1,include_singles):
        a2 = -1
        b2 = len_only1*2+len_common*2+(1 if include_singles else -1)
        c = -count1
        row_index = int((b2-(b2*b2-8*a2*c)**.5)/2.0) #!!!cmk will this ever be off by one because of numerical problems?
        if row_index < 0:
            row_index = int((b2+(b2*b2-8*a2*c)**.5)/2.0)
        return row_index

class TestPairs(unittest.TestCase):

    def test_count1_fun(self):
        for len_common in range(5):
            for len_only1 in range(5):
                for include_singles in [True,False]:
                    list0 = ['common{0}'.format(i) for i in range(len_common)]
                    list1 = list0 + ['only1_{0}'.format(i) for i in range(len_only1)]
                    pairs = Pairs(list0,list1,include_singles)
                    pair_list = list(pairs.pair_sequence(0,speed='slow'))
                    for row_index in range(len_common+1):
                        count = Pairs.count1_fun(row_index,len_common,len_only1,include_singles)
                        count_ref = len([pair for pair in pair_list if pair[0] in list0[:row_index]])
                        assert count==count_ref, "count1_fun isn't giving the right answer"

    def test_row_index_fun(self):
        for len_common in range(5):
            for len_only1 in range(5):
                for include_singles in [True,False]:
                    list0 = ['common{0}'.format(i) for i in range(len_common)]
                    list1 = list0 + ['only1_{0}'.format(i) for i in range(len_only1)]
                    pairs = Pairs(list0,list1,include_singles)
                    pair_list = list(pairs.pair_sequence(0,speed='slow'))
                    for start in range(pairs.count+1):
                        row_index = pairs.row_index_fun(start,len_common,len_only1,include_singles)
                        assert start==len(pair_list) or pair_list[start][0]==list0[row_index], "row_index_fun isn't giving the right answer"


    def test_index_and_count1_functions_spot_test(self):
        index = Pairs.row_index_fun(3,2,2,True)
        count1 = Pairs.count1_fun(index,2,2,True)
        assert count1 <= 3


    def test_index_and_count1_functions(self):
        for len_commonq in range(0,5):
            for len_only1q in range(0,5):
                for include_singlesq in [True,False]:
                    for row_indexq in range(len_commonq+1):
                        count1q =Pairs.count1_fun(row_indexq,len_commonq,len_only1q,include_singlesq)
                        row_index2 = Pairs.row_index_fun(count1q,len_commonq,len_only1q,include_singlesq)
                        count2 = Pairs.count1_fun(row_index2,len_commonq,len_only1q,include_singlesq)
                        assert count1q==count2

    def test_index_and_count1_functions_1_0_False_1(self):
        len_commonq,len_only1q,include_singlesq,row_indexq = 1,0,False,1
        count1q = Pairs.count1_fun(row_indexq,len_commonq,len_only1q,include_singlesq)
        row_index2 = Pairs.row_index_fun(count1q,len_commonq,len_only1q,include_singlesq)
        count2 = Pairs.count1_fun(row_index2,len_commonq,len_only1q,include_singlesq)
        assert count1q==count2


    def test_index_and_count1_functions_0_0_False_0(self):
        len_commonq,len_only1q,include_singlesq,row_indexq = 0,0,False,0
        count1q =Pairs.count1_fun(row_indexq,len_commonq,len_only1q,include_singlesq)
        row_index2 = Pairs.row_index_fun(count1q,len_commonq,len_only1q,include_singlesq)
        assert row_indexq==row_index2

    def test_medium_fast(self):
        for len_common in range(2):
            common = ['common{0}'.format(i) for i in range(len_common)]
            for len_only0 in range(2):
                list0 = common + ['only0_{0}'.format(i) for i in range(len_only0)]
                for len_only1 in range(2):
                    list1 = common + ['only1_{0}'.format(i) for i in range(len_only1)]
                    for include_singles in [True,False]:
                        pairs = Pairs(list0, list1, include_singles, duplicates_ok=True)
                        for goal_see in range(pairs.count+1):
                            for start in range(pairs.count+1):
                                logging.info((len_common,len_only0,len_only1,include_singles,start,start+goal_see))
                                slow = np.array(list(pairs.pair_sequence(start,start+goal_see,speed='slow')))
                                medium = np.array(list(pairs.pair_sequence(start,start+goal_see,speed='medium')))
                                assert np.array_equal(slow,medium), 'Expect slow and medium to give the same answers'
                                fast = np.array(list(pairs.pair_sequence(start,start+goal_see,speed='fast')))
                                assert np.array_equal(medium,fast), 'Expect medium and fast to give the same answers'

    def test_11_12_True(self):
        self.little(11,12,include_singles=True)

    def test_12_13_False(self):
        self.little(12,13,include_singles=False)

    def little(self,start, goal_see,include_singles):
        list0 = ['a','z','a','b','y']
        list1 = ['z','y','x','v']
        
        pairs = Pairs(list0, list1, include_singles, duplicates_ok=True)
        fast = np.array(list(pairs.pair_sequence(start,start+goal_see,speed='fast')))
        medium = np.array(list(pairs.pair_sequence(start,start+goal_see,speed='medium')))
        assert np.array_equal(medium,fast), 'Expect medium and fast to give the same answers'

def getTestSuite():
    """
    set up composite test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPairs))
    return test_suite


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    TestPairs().test_index_and_count1_functions_1_0_False_1()
    TestPairs().test_index_and_count1_functions()
    TestPairs().test_row_index_fun()
    TestPairs().test_count1_fun()
    TestPairs().test_12_13_False()
    #TestPairs().test_index_and_count1_functions_spot_test()
    #TestPairs().test_11_12_True()
    TestPairs().test_medium_fast()

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=True)
    r.run(suites)

    result = doctest.testmod()
    assert result.failed == 0, "failed doc test: " + __file__

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

    speed = 'slow'
    for goal_see in [1,2,4,8,15,16]:
        for start in [0,pairs.count//5,pairs.count//2,pairs.count-1,pairs.count]:
            print(start,goal_see,list(pairs.pair_sequence(start,start+goal_see,speed=speed))) #!!!cmk instead of "pair_sequence" how about an indexer? (and "pair" is redundant)


    print("done")
