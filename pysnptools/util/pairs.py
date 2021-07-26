import numpy as np
from itertools import chain,islice
from  more_itertools import unique_everseen
import doctest
import unittest
import logging
import numbers
import pysnptools.util as pstutil

class _Pairs(object):
    """!!!description of class"""

    def __init__(self, list0, list1, include_singles, duplicates_ok=True): #!!! rename list to iterable?
        self.include_singles = include_singles
        self.duplicates_ok = duplicates_ok
        self.list0 = list(unique_everseen(list0))
        self.list1 = list(unique_everseen(list1))
        if not duplicates_ok:
            assert len(list0)==len(self.list0) and len(list1)==len(self.list1), "Expect lists to have no duplicates"

        dict0 = {k: v for v, k in enumerate(self.list0)}
        dict1 = {k: v for v, k in enumerate(self.list1)}
        self._common = sorted(six.viewkeys(dict0) & self.list1, key=lambda k:dict0[k])
        self._only0 = sorted(six.viewkeys(dict0)-self._common,key=lambda k:dict0[k])
        self._only1 = sorted(six.viewkeys(dict1)-self._common,key=lambda k:dict1[k])
    
        #How many pairs?
        self._count0 = len(self._only0)*len(self.list1)
        self._count1 = self._count1_fun(len(self._common),len(self._common),len(self._only1),include_singles)
        self._count = self._count0+self._count1

    def __len__(self):
        return self._count

    def __getitem__(self, indexer,speed='fast'):
        if isinstance(indexer,slice):
            assert indexer.step is None, "step is not supported"
            start,stop = indexer.start,indexer.stop
            if start is None:
                start = 0
            if stop is None:
                stop = self._count
            if start < 0:
                start = max(0,self._count+start)
            if stop < 0:
                stop = self._count+stop
            for pair in self._pair_sequence(start,stop,speed=speed):
                yield pair
        elif isinstance(indexer, numbers.Integral):
            start = indexer
            if start < 0:
                start = self._count+start
            assert 0<=start and start<self._count, 'Expect index between 0 (inclusive) and count of pairs (exclusive)'
            stop = start+1
            for pair in self._pair_sequence(start,stop,speed=speed):
                yield pair
        elif isinstance(indexer,np.ndarray): #!!! what about other arrays? iterator? negative elements?
            #!!! need test for this
            iterator = None
            previous = None
            for i in indexer:
                if iterator is None or i != previous+1:
                    iterator = self._pair_sequence(i,None,speed=speed)
                yield next(iterator)
                previous = i
        else:
            raise Exception("Don't understand indexer '{0}'".format(indexer))


    def _pair_sequence(self, start=0,stop=None,speed='fast'):
        start = int(start)
        if speed=='fast':
            stop = min(stop,self._count) if stop is not None else self._count
            return islice(self._pair_sequence_inner(start=start),max(stop-start,0))
        if speed=='medium':
            stop = stop if stop is not None else self._count
            return islice(self._pair_sequence_inner_medium(start=start),max(stop-start,0))
        else:
            assert speed=='slow',"Don't know speed '{0}'".format(speed)
            return islice(self._pair_sequence_inner_slow(),start,stop)

    def _pair_sequence_inner(self,start=0):
        if start >= self._count0:
            start -= self._count0
        else:
            row_start0 = start // len(self.list1)
            start -= row_start0*len(self.list1)
            assert 0 <= start and start <= len(self.list1), "real assert"
            for v0 in self._only0[row_start0:]:
                for v1 in islice(chain(self._only1,self._common),start,None):
                    yield v0,v1
                start = 0
        
        row_start1 = self._row_start1_fun(start,len(self._common),len(self._only1),self.include_singles)
        start -= self._count1_fun(row_start1,len(self._common),len(self._only1),self.include_singles)
        assert 0 <= start and start <= len(self.list1), "real assert"
        for common_index in range(row_start1,len(self._common)):
            v0 = self._common[common_index]
            common_start = common_index if self.include_singles else common_index+1
            assert start < len(self.list1)-common_start, "real assert"
            for v1 in islice(chain(self._only1,self._common[common_start:]),start,None):
                yield v0,v1
            start = 0

    def _pair_sequence_inner_medium(self,start=0):
        assert 0 <= start, "real assert"
        for v0 in self._only0:
            if start > len(self.list1):
                start -= len(self.list1)
            else:
                for v1 in islice(chain(self._only1,self._common),start,None):
                    yield v0,v1
                start = 0

        assert 0 <= start, "real assert"
        for index in range(len(self._common)):
            v0 = self._common[index]
            startx = index if self.include_singles else index+1
            if start > len(self.list1)-startx:
                start -= (len(self.list1)-startx)
            else:
                for v1 in islice(chain(self._only1,self._common[startx:]),start,None):
                    yield v0,v1
                start = 0

    def _pair_sequence_inner_slow(self):

        for v0 in self._only0:
            for v1 in chain(self._only1,self._common):
                yield v0,v1

        for index in range(len(self._common)):
            v0 = self._common[index]
            startx = index if self.include_singles else index+1
            for v1 in chain(self._only1,self._common[startx:]):
                yield v0,v1


    @staticmethod
    def _count1_fun(row_start,len_common,len_only1,include_singles):
        a2 = -1
        b2 = len_only1*2+len_common*2+(1 if include_singles else -1)
        count1 = (a2*row_start*row_start + b2*row_start)//2
        return count1

    @staticmethod
    def _row_start1_fun(count1,len_common,len_only1,include_singles):
        a2 = -1
        b2 = len_only1*2+len_common*2+(1 if include_singles else -1)
        c = -count1
        row_start = int((b2-(b2*b2-8*a2*c)**.5)/2.0) #!!! will this ever be off by one because of numerical problems?
        if row_start < 0:
            row_start = int((b2+(b2*b2-8*a2*c)**.5)/2.0)
        return row_start

class TestPairs(unittest.TestCase):

    def test_count1_fun(self):
        for len_common in range(5):
            for len_only1 in range(5):
                for include_singles in [True,False]:
                    list0 = ['common{0}'.format(i) for i in range(len_common)]
                    list1 = list0 + ['only1_{0}'.format(i) for i in range(len_only1)]
                    pairs = _Pairs(list0,list1,include_singles)
                    pair_list = list(pairs._pair_sequence(0,speed='slow'))
                    for row_start in range(len_common+1):
                        count = _Pairs._count1_fun(row_start,len_common,len_only1,include_singles)
                        count_ref = len([pair for pair in pair_list if pair[0] in list0[:row_start]])
                        assert count==count_ref, "_count1_fun isn't giving the right answer"

    def test_row_start_fun(self):
        for len_common in range(5):
            for len_only1 in range(5):
                for include_singles in [True,False]:
                    list0 = ['common{0}'.format(i) for i in range(len_common)]
                    list1 = list0 + ['only1_{0}'.format(i) for i in range(len_only1)]
                    pairs = _Pairs(list0,list1,include_singles)
                    pair_list = list(pairs._pair_sequence(0,speed='slow'))
                    for start in range(len(pairs)+1):
                        row_start = pairs._row_start1_fun(start,len_common,len_only1,include_singles)
                        assert start==len(pair_list) or pair_list[start][0]==list0[row_start], "_row_start1_fun isn't giving the right answer"

    def test_index_and_count1_functions(self):
        for len_commonq in range(0,5):  
            for len_only1q in range(0,5):
                for include_singlesq in [True,False]:
                    for row_startq in range(len_commonq+1):
                        count1q =_Pairs._count1_fun(row_startq,len_commonq,len_only1q,include_singlesq)
                        row_start2 = _Pairs._row_start1_fun(count1q,len_commonq,len_only1q,include_singlesq)
                        count2 = _Pairs._count1_fun(row_start2,len_commonq,len_only1q,include_singlesq)
                        assert count1q==count2

    def test_medium_fast(self):
        for len_common in range(2):
            common = ['common{0}'.format(i) for i in range(len_common)]
            for len_only0 in range(2):
                list0 = common + ['only0_{0}'.format(i) for i in range(len_only0)]
                for len_only1 in range(2):
                    list1 = common + ['only1_{0}'.format(i) for i in range(len_only1)]
                    for include_singles in [True,False]:
                        pairs = _Pairs(list0, list1, include_singles, duplicates_ok=True)
                        for goal_see in range(len(pairs)+1):
                            for start in range(len(pairs)+1):
                                logging.info((len_common,len_only0,len_only1,include_singles,start,start+goal_see))
                                slow = np.array(list(pairs._pair_sequence(start,start+goal_see,speed='slow')))
                                medium = np.array(list(pairs._pair_sequence(start,start+goal_see,speed='medium')))
                                assert np.array_equal(slow,medium), 'Expect slow and medium to give the same answers'
                                fast = np.array(list(pairs[start:start+goal_see]))
                                assert np.array_equal(medium,fast), 'Expect medium and fast to give the same answers'
    def test_big(self):
        size = 500*1000
        seed = 0
        np.random.seed(seed)
        list0 = np.random.randint(size*10,size=size)
        list1 = np.random.randint(size*10,size=size)
        for include_singles in [True,False]:
            pairs = _Pairs(list0, list1, include_singles, duplicates_ok=True)
            for start in [0,len(pairs)//5,len(pairs)-1,len(pairs)]:
                logging.info(("test_big",include_singles,start))
                medium = np.array(list(pairs._pair_sequence(start,start+10,speed='medium')))
                fast = np.array(list(pairs[start:start+10]))
                assert np.array_equal(medium,fast), 'Expect medium and fast to give the same answers'


    #!!! make doc string example with first and last names
    def test_slice(self):
        for len_common in range(3):
            common = ['common{0}'.format(i) for i in range(len_common)]
            for len_only0 in range(3):
                list0 = common + ['only0_{0}'.format(i) for i in range(len_only0)]
                for len_only1 in range(3):
                    list1 = common + ['only1_{0}'.format(i) for i in range(len_only1)]
                    for include_singles in [True,False]:
                        pairs = _Pairs(list0, list1, include_singles, duplicates_ok=True)
                        pair_list = np.array(list(pairs._pair_sequence(0,None,speed='slow')),dtype='U')
                        for slicer, start, stop in [
                            (-1,len(pairs)-1,len(pairs)),#[-1]
                            (slice(None,None),0,len(pairs)),#[:]
                            (slice(None,2),0,2),#[:2]
                            (2,2,3),#[2]
                            (slice(2,None),2,None),#[2:]
                            (slice(1,3),1,3),#[1:3]
                            (slice(3,3),3,3),#[3:3]
                            (slice(-3,-2),len(pairs)-3,len(pairs)-2),#[-3:-2]
                            (slice(-3,None),len(pairs)-3,None),#[-3:],
                            (slice(None,-3),None,len(pairs)-3),#[:-3],
                            (slice(1000,None),1000,None),#[1000:]
                            (slice(-3,None),len(pairs)-100000,None),#[-100000:]
                            (slice(4,1),4,1),#[4:1],
                            (slice(10,1000),10,1000)#[10,1000]
                            ]:
                            logging.info(("test_slice",slicer, start, stop,len_common,len_only0,len_only1,include_singles))
                            try:
                                slow = pair_list[slicer]
                                slow_fail = False
                            except:
                                slow_fail = True
                            try:
                                fast = np.array(list(pairs[slicer]),dtype='U')
                                fast_fail = False
                            except:
                                fast_fail = True
                            if slow_fail:
                                assert fast_fail, 'expect slow and fast to fail on the same inputs'
                            else:
                                if fast_fail or not np.array_equal(slow.reshape(-1,2),fast.reshape(-1,2)): #!!!:
                                    fast = np.array(list(pairs[slicer]))
                                assert not fast_fail, 'expect slow and fast to fail on the same inputs'
                                assert np.array_equal(slow.reshape(-1,2),fast.reshape(-1,2)), 'Expect slow and fast to give the same answers'

def getTestSuite():
    """
    set up composite test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPairs))
    return test_suite


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    #!!! this is work on an example for docstring, but a bit weird because names are ordered and these arent
    pairs = _Pairs(['Madison','Jill','Kelsey'],['Smith','Madison','Kelsey'],include_singles=True)#!!! is 'include_singles' the best name?
    for first,last in pairs[:]:
        print(first,last)

    pairs = _Pairs(range(0,2000),range(1000,3000),include_singles=True)
    print(len(pairs))
    for pair in pairs[1999999:2000004]:
        print(pair)
    print(list(pairs[-5:]))


    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=True)
    r.run(suites)

    result = doctest.testmod()
    assert result.failed == 0, "failed doc test: " + __file__


    print("done")