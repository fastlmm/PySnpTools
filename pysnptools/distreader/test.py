import numpy as np
import scipy as sp
import logging
import doctest
import shutil
import unittest
import os.path
import time

from pysnptools.distreader.distmemmap import TestDistMemMap
from pysnptools.distreader.bgen import TestBgen
from pysnptools.distreader.distgen import TestDistGen
from pysnptools.distreader import DistNpz, DistHdf5, DistMemMap, DistData, _DistMergeSIDs
from pysnptools.util import create_directory_if_necessary
from pysnptools.snpreader import SnpNpz, Bed
from pysnptools.kernelreader.test import _fortesting_JustCheckExists
from pysnptools.pstreader import PstData
from pysnptools.snpreader import SnpData

class TestDistReaders(unittest.TestCase):     


    @classmethod
    def setUpClass(self):
        self.currentFolder = os.path.dirname(os.path.realpath(__file__))
        #TODO: get data set with NANs!
        distreader = DistNpz(self.currentFolder + "/../examples/toydata.dist.npz")
        self.distdata = distreader.read(order='F',force_python_only=True)
        self.dist_values = self.distdata.val

    def test_3d(self):
        from pysnptools.distreader import DistData
        np.random.seed(0)
        row_count = 4
        col_count = 5
        val_shape = 3
        val = np.random.random((row_count,col_count,val_shape))
        val /= val.sum(axis=2,keepdims=True)  #make probabilities sum to 1
        distdata = DistData(val=val,iid=[['iid{0}'.format(i)]*2 for i in range(row_count)],sid=['sid{0}'.format(s) for s in range(col_count)]
                            )

    def test_scalar_index(self):
        distreader = DistNpz(self.currentFolder + "/../examples/toydata.dist.npz")
        arr=np.int64(1)
        distreader[arr,arr]

    def test_c_reader_hdf5(self):
        distreader = DistHdf5(self.currentFolder + "/../examples/toydata.snpmajor.dist.hdf5")
        self.c_reader(distreader)

    def test_hdf5_case3(self):
        distreader1 = DistHdf5(self.currentFolder + "/../examples/toydata.snpmajor.dist.hdf5")[::2,:]
        distreader2 = DistNpz(self.currentFolder + "/../examples/toydata.dist.npz")[::2,:]
        self.assertTrue(np.allclose(distreader1.read().val, distreader2.read().val, rtol=1e-05, atol=1e-05))


    def test_c_reader_npz(self):
        distreader = DistNpz(self.currentFolder + "/../examples/toydata10.dist.npz")
        distdata = distreader.read(order='F',force_python_only=False)
        snp_c = distdata.val
        
        self.assertEqual(np.float64, snp_c.dtype)
        self.assertTrue(np.allclose(self.dist_values[:,:10], snp_c, rtol=1e-05, atol=1e-05))

        distreader1 = DistNpz(self.currentFolder + "/../examples/toydata10.dist.npz")
        distreader2 = DistNpz(self.currentFolder + "/../examples/toydata.dist.npz")[:,:10]
        self.assertTrue(np.allclose(distreader1.read().val, distreader2.read().val, rtol=1e-05, atol=1e-05))


        distdata.val[1,2] = np.NaN # Inject a missing value to test writing and reading missing values
        output = "tempdir/distreader/toydata10.dist.npz"
        create_directory_if_necessary(output)
        DistNpz.write(output,distdata)
        snpdata2 = DistNpz(output).read()
        np.testing.assert_array_almost_equal(distdata.val, snpdata2.val, decimal=10)

        snpdata3 = distdata[:,0:0].read() #create distdata with no sids
        output = "tempdir/distreader/toydata0.dist.npz"
        DistNpz.write(output,snpdata3)
        snpdata4 = DistNpz(output).read()
        assert snpdata3 == snpdata4

               
    def c_reader(self,distreader):
        """
        make sure c-reader yields same result
        """
        distdata = distreader.read(order='F',force_python_only=False)
        snp_c = distdata.val
        
        self.assertEqual(np.float64, snp_c.dtype)
        self.assertTrue(np.allclose(self.dist_values, snp_c, rtol=1e-05, atol=1e-05))
        return distdata

    def test_write_distnpz_f64cpp_0(self):
        distreader = DistNpz(self.currentFolder + "/../examples/toydata.dist.npz")
        iid_index = 0
        logging.info("iid={0}".format(iid_index))
        #if distreader.iid_count % 4 == 0: # divisible by 4 isn't a good test
        #    distreader = distreader[0:-1,:]
        #assert distreader.iid_count % 4 != 0
        distdata = distreader[0:iid_index,:].read(order='F',dtype=np.float64)
        if distdata.iid_count > 0:
            distdata.val[-1,0] = float("NAN")
        output = "tempdir/toydata.F64cpp.{0}.dist.npz".format(iid_index)
        create_directory_if_necessary(output)
        DistNpz.write(output, distdata )
        snpdata2 = DistNpz(output).read()
        np.testing.assert_array_almost_equal(distdata.val, snpdata2.val, decimal=10)

    def test_write_distnpz_f64cpp_1(self):
        distreader = DistNpz(self.currentFolder + "/../examples/toydata.dist.npz")
        iid_index = 1
        logging.info("iid={0}".format(iid_index))
        #if distreader.iid_count % 4 == 0: # divisible by 4 isn't a good test
        #    distreader = distreader[0:-1,:]
        #assert distreader.iid_count % 4 != 0
        distdata = distreader[0:iid_index,:].read(order='F',dtype=np.float64)
        if distdata.iid_count > 0:
            distdata.val[-1,0] = float("NAN")
        output = "tempdir/toydata.F64cpp.{0}.dist.npz".format(iid_index)
        create_directory_if_necessary(output)
        DistNpz.write(output, distdata )
        snpdata2 = DistNpz(output).read()
        np.testing.assert_array_almost_equal(distdata.val, snpdata2.val, decimal=10)

    def test_write_distnpz_f64cpp_5(self):
        distreader = DistNpz(self.currentFolder + "/../examples/toydata.dist.npz")

        _fortesting_JustCheckExists().input(distreader)

        iid_index = 5
        logging.info("iid={0}".format(iid_index))
        #if distreader.iid_count % 4 == 0: # divisible by 4 isn't a good test
        #    distreader = distreader[0:-1,:]
        #assert distreader.iid_count % 4 != 0
        distdata = distreader[0:iid_index,:].read(order='F',dtype=np.float64)
        if distdata.iid_count > 0:
            distdata.val[-1,0] = float("NAN")
        output = "tempdir/toydata.F64cpp.{0}.dist.npz".format(iid_index)
        create_directory_if_necessary(output)
        DistNpz.write(output, distdata ) #,force_python_only=True)
        snpdata2 = DistNpz(output).read()
        np.testing.assert_array_almost_equal(distdata.val, snpdata2.val, decimal=10)



    def test_write_x_x_cpp(self):
        distreader = DistNpz(self.currentFolder + "/../examples/toydata.dist.npz")
        for order in ['C','F']:
            for dtype in [np.float32,np.float64]:
                distdata = distreader.read(order=order,dtype=dtype)
                distdata.val[-1,0] = float("NAN")
                output = "tempdir/toydata.{0}{1}.cpp.dist.npz".format(order,"32" if dtype==np.float32 else "64")
                create_directory_if_necessary(output)
                DistNpz.write(output, distdata)
                snpdata2 = DistNpz(output).read()
                np.testing.assert_array_almost_equal(distdata.val, snpdata2.val, decimal=10)

    def test_subset_view(self):
        distreader2 = DistNpz(self.currentFolder + "/../examples/toydata.dist.npz")[:,:]
        result = distreader2.read(view_ok=True)
        self.assertFalse(distreader2 is result)
        result2 = result[:,:].read()
        self.assertFalse(np.may_share_memory(result2.val,result.val))
        result3 = result[:,:].read(view_ok=True)
        self.assertTrue(np.may_share_memory(result3.val,result.val))
        result4 = result3.read()
        self.assertFalse(np.may_share_memory(result4.val,result3.val))
        result5 = result4.read(view_ok=True)
        self.assertTrue(np.may_share_memory(result4.val,result5.val))




    def test_writes(self):
        from pysnptools.distreader import DistData, DistHdf5, DistNpz, DistMemMap, Bgen
        from pysnptools.kernelreader.test import _fortesting_JustCheckExists

        the_class_and_suffix_list = [(DistNpz,"npz",None,None),
                                     (Bgen,"bgen",None,lambda filename,distdata: Bgen.write(filename,distdata,bits=32)),
                                     (DistHdf5,"hdf5",None,None),
                                     (DistMemMap,"memmap",None,None)]
        cant_do_col_prop_none_set = {'bgen'}
        cant_do_col_len_0_set = {'bgen'}
        cant_do_row_count_zero_set = {'bgen'}
        can_swap_0_2_set = {}
        can_change_col_names_set = {}
        ignore_fam_id_set = {}
        ignore_pos1_set = {'bgen'}
        ignore_pos_set = {}
        erase_any_write_dir = {}

        
        #===================================
        #    Starting main function
        #===================================
        logging.info("starting 'test_writes'")
        np.random.seed(0)
        output_template = "tempdir/distreader/writes.{0}.{1}"
        create_directory_if_necessary(output_template.format(0,"npz"))
        i = 0
        for row_count in [0,5,2,1]:
            for col_count in [4,2,1,0]:
                val=np.random.random(size=[row_count,col_count,3])
                val /= val.sum(axis=2,keepdims=True)  #make probabilities sum to 1

                val[val==3]=np.NaN
                row = [('0','0'),('1','1'),('2','2'),('3','3'),('4','4')][:row_count]
                col = ['s0','s1','s2','s3','s4'][:col_count]
                for is_none in [True,False]:
                    row_prop = None
                    col_prop = None if is_none else [(x,x,x) for x in range(5)][:col_count]
                    distdata = DistData(iid=row,sid=col,val=val,pos=col_prop,name=str(i))
                    for the_class,suffix,constructor,writer in the_class_and_suffix_list:
                        constructor = constructor or (lambda filename: the_class(filename))
                        writer = writer or (lambda filename,distdata: the_class.write(filename,distdata))
                        if col_count == 0 and suffix in cant_do_col_len_0_set:
                            continue
                        if col_prop is None and suffix in cant_do_col_prop_none_set:
                            continue
                        if row_count==0 and suffix in cant_do_row_count_zero_set:
                            continue
                        filename = output_template.format(i,suffix)
                        logging.info(filename)
                        i += 1
                        if suffix in erase_any_write_dir and os.path.exists(filename):
                            shutil.rmtree(filename)
                        ret = writer(filename,distdata)
                        assert ret is not None
                        for subsetter in [None, np.s_[::2,::3]]:
                            reader = constructor(filename)
                            _fortesting_JustCheckExists().input(reader)
                            subreader = reader if subsetter is None else reader[subsetter[0],subsetter[1]]
                            readdata = subreader.read(order='C')
                            expected = distdata if subsetter is None else distdata[subsetter[0],subsetter[1]].read()
                            if not suffix in can_swap_0_2_set:
                                assert np.allclose(readdata.val,expected.val,equal_nan=True)
                            else:
                                for col_index in range(readdata.col_count):
                                    assert (np.allclose(readdata.val[:,col_index],expected.val[:,col_index],equal_nan=True) or
                                            np.allclose(readdata.val[:,col_index]*-1+2,expected.val[:,col_index],equal_nan=True))
                            if not suffix in ignore_fam_id_set:
                                assert np.array_equal(readdata.row,expected.row)
                            else:
                                assert np.array_equal(readdata.row[:,1],expected.row[:,1])
                            if not suffix in can_change_col_names_set:
                                assert np.array_equal(readdata.col,expected.col)
                            else:
                                assert readdata.col_count==expected.col_count
                            assert np.array_equal(readdata.row_property,expected.row_property) or (readdata.row_property.shape[1]==0 and expected.row_property.shape[1]==0)

                            if suffix in ignore_pos1_set:
                                assert np.allclose(readdata.col_property[:,[0,2]],expected.col_property[:,[0,2]],equal_nan=True) or (readdata.col_property.shape[1]==0 and expected.col_property.shape[1]==0)
                            elif not suffix in ignore_pos_set:
                                assert np.allclose(readdata.col_property,expected.col_property,equal_nan=True) or (readdata.col_property.shape[1]==0 and expected.col_property.shape[1]==0)
                            else:
                                assert len(readdata.col_property)==len(expected.col_property)
                        try:
                            os.remove(filename)
                        except:
                            pass
        logging.info("done with 'test_writes'")

    def test_block_size_DistData(self):
        np.random.seed(0)
        sid_count = 20
        val=np.array(np.random.random(size=[3,sid_count,3]),dtype=np.float64,order='F')
        val /= val.sum(axis=2,keepdims=True)  #make probabilities sum to 1
        distreader = DistData(iid=[["0","0"],["1","1"],["2","2"]],sid=[str(i) for i in range(sid_count)],val=val)
        snpdata0 = distreader.as_snp(max_weight=100,block_size=1).read()
        snpdata1 = distreader.as_snp(max_weight=100,block_size=None).read()
        np.testing.assert_array_almost_equal(snpdata0.val,snpdata1.val, decimal=10)

    def test_block_size_Snp2Dist(self):
        from pysnptools.snpreader import SnpData
        from pysnptools.distreader._snp2dist import _Snp2Dist
        np.random.seed(0)
        sid_count = 20
        val=np.array(np.random.randint(0,3,size=[3,sid_count]),dtype=np.float64,order='F')
        snpreader = SnpData(iid=[["0","0"],["1","1"],["2","2"]],sid=[str(i) for i in range(sid_count)],val=val)
        distdata0 = snpreader.as_dist(max_weight=2,block_size=1).read()
        distdata1 = snpreader.as_dist(max_weight=2,block_size=None).read()
        np.testing.assert_array_almost_equal(distdata0.val,distdata1.val, decimal=10)

    def test_intersection_Dist2Snp(self):
        from pysnptools.snpreader._dist2snp import _Dist2Snp
        from pysnptools.snpreader import Pheno
        from pysnptools.distreader._subset import _DistSubset
        from pysnptools.snpreader._subset import _SnpSubset
        from pysnptools.util import intersect_apply

        dist_all = DistNpz(self.currentFolder + "/../examples/toydata.dist.npz")
        k = dist_all.as_snp(max_weight=25)

        pheno = Pheno(self.currentFolder + "/../examples/toydata.phe")
        pheno = pheno[1:,:] # To test intersection we remove a iid from pheno

        k1,pheno = intersect_apply([k,pheno]) 
        assert isinstance(k1.distreader,_DistSubset) and not isinstance(k1,_SnpSubset)

        #What happens with fancy selection?
        k2 = k[::2,:]
        assert isinstance(k2,_Dist2Snp)

        logging.info("Done with test_intersection")

    def test_intersection_Snp2Dist(self):
        from pysnptools.distreader._snp2dist import _Snp2Dist
        from pysnptools.snpreader import Pheno, Bed
        from pysnptools.distreader._subset import _DistSubset
        from pysnptools.snpreader._subset import _SnpSubset
        from pysnptools.util import intersect_apply

        snp_all = Bed(self.currentFolder + "/../examples/toydata.5chrom.bed",count_A1=True)
        k = snp_all.as_dist(max_weight=2)

        pheno = Pheno(self.currentFolder + "/../examples/toydata.phe")
        pheno = pheno[1:,:] # To test intersection we remove a iid from pheno

        k1,pheno = intersect_apply([k,pheno]) 
        assert isinstance(k1.snpreader,_SnpSubset) and not isinstance(k1,_DistSubset)

        #What happens with fancy selection?
        k2 = k[::2,:]
        assert isinstance(k2,_Snp2Dist)

        logging.info("Done with test_intersection")

    def test_roundtrip(self):
        max_weight = 2
        snpreader1 = Bed(self.currentFolder + "/../examples/toydata.5chrom.bed",count_A1=True)
        snpdata1 = snpreader1.read()
        distreader1 = snpreader1.as_dist(max_weight)
        snpreader2 = distreader1.as_snp(max_weight)
        assert snpdata1.allclose(snpreader2.read(),equal_nan=True)
        snpdata1.val[0,0] = np.nan
        assert snpdata1.allclose(snpdata1.as_dist(max_weight).as_snp(max_weight).read(),equal_nan=True)

    def test_dist_snp2(self):
        logging.info("in test_dist_snp2")
        distreader = DistNpz(self.currentFolder + "/../examples/toydata.dist.npz")
        dist2snp = distreader.as_snp(max_weight=33)
        s = str(dist2snp)
        _fortesting_JustCheckExists().input(dist2snp)

    def test_snp_dist2(self):
        logging.info("in test_snp_dist2")
        snpreader = Bed(self.currentFolder + "/../examples/toydata.5chrom.bed",count_A1=False)
        snp2dist = snpreader.as_dist(max_weight=2)
        s = str(snp2dist)
        _fortesting_JustCheckExists().input(snp2dist)

    def test_subset_Dist2Snp(self):
        logging.info("in test_subset")
        distreader = DistNpz(self.currentFolder + "/../examples/toydata.dist.npz")
        dist2snp = distreader.as_snp(max_weight=10)
        dssub = dist2snp[::2,::2]
        snpdata1 = dssub.read()
        expected = distreader.as_snp(max_weight=10)[::2,::2].read()
        np.testing.assert_array_almost_equal(snpdata1.val, expected.val, decimal=10)

        logging.info("done with test")

    def test_subset_Snp2Dist(self): #!!!move these to another test class
        logging.info("in test_subset")
        snpreader = Bed(self.currentFolder + "/../examples/toydata.5chrom.bed",count_A1=False)
        snp2dist = snpreader.as_dist(max_weight=2)
        sub = snp2dist[::2,::2]
        distdata1 = sub.read()
        expected = snpreader.as_dist(max_weight=2)[::2,::2].read()
        np.testing.assert_array_almost_equal(distdata1.val, expected.val, decimal=10)

        logging.info("done with test")

    def test_respect_read_inputs(self):
        from pysnptools.distreader import Bgen,DistGen,DistHdf5,DistMemMap,DistNpz
        from pysnptools.snpreader import Bed

        previous_wd = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        for distreader in [
                           _DistMergeSIDs([Bgen('../examples/example.bgen')[:,:5].read(),Bgen('../examples/example.bgen')[:,5:].read()]),
                           Bed('../examples/toydata.5chrom.bed',count_A1=True).as_dist(block_size=2000),
                           Bed('../examples/toydata.5chrom.bed',count_A1=True).as_dist(),
                           Bgen('../examples/example.bgen').read(),
                           Bgen('../examples/bits1.bgen'),                          
                           DistGen(seed=0,iid_count=500,sid_count=50),
                           DistGen(seed=0,iid_count=500,sid_count=50)[::2,::2],
                           DistHdf5('../examples/toydata.snpmajor.dist.hdf5'),
                           DistMemMap('../examples/tiny.dist.memmap'),
                           DistNpz('../examples/toydata10.dist.npz')
                          ]:
            logging.info(str(distreader))
            for order in ['F','C','A']:
                for dtype in [np.float32,np.float64]:
                    for force_python_only in [True,False]:
                        for view_ok in [True,False]:
                            val = distreader.read(order=order,dtype=dtype,force_python_only=force_python_only,view_ok=view_ok).val
                            has_right_order = order=="A" or (order=="C" and val.flags["C_CONTIGUOUS"]) or (order=="F" and val.flags["F_CONTIGUOUS"])
                            if hasattr(distreader,'val') and not view_ok:
                                assert distreader.val is not val
                            if (hasattr(distreader,'val') and view_ok and distreader.val is not val and
                                (order == 'A' or (order == 'F' and distreader.val.flags['F_CONTIGUOUS']) or (order == 'C' and distreader.val.flags['C_CONTIGUOUS'])) and
                                (dtype is None or  distreader.val.dtype == dtype)):
                                logging.info("{0} could have read a view, but didn't".format(distreader))
                            assert val.dtype == dtype and has_right_order

        os.chdir(previous_wd)


    def test_respect_inputs_DistData(self):
        np.random.seed(0)
        for dtype_start,decimal_start in [(np.float32,5),(np.float64,10)]:
            for order_start in ['F','C','A']:
                for sid_count in [20,2]:
                    val=np.array(np.random.random(size=[3,sid_count,3]),dtype=dtype_start,order=order_start)
                    val /= val.sum(axis=2,keepdims=True)  #make probabilities sum to 1
                    distdataX = DistData(iid=[["0","0"],["1","1"],["2","2"]],sid=[str(i) for i in range(sid_count)],val=val)
                    for max_weight in [1.0,2.0]:
                        weights = np.array([0,.5,1])*max_weight
                        for distreader0 in [distdataX,distdataX[:,1:]]:
                            distreader1 = distreader0[1:,:]
                            refdata0 = distreader0.read()
                            refval0 = (refdata0.val * weights).sum(axis=-1)
                            for dtype_goal,decimal_goal in [(np.float32,5),(np.float64,10)]:
                                for order_goal in ['F','C','A']:
                                    k = distreader0.as_snp(max_weight=max_weight,block_size=1).read(order=order_goal,dtype=dtype_goal)
                                    DistData._array_properties_are_ok(k.val,order_goal,dtype_goal)
                                    np.testing.assert_array_almost_equal(refval0,k.val, decimal=min(decimal_start,decimal_goal))


    def test_npz(self):
        logging.info("in test_npz")
        distreader = DistNpz(self.currentFolder + "/../examples/toydata.dist.npz")
        snpdata1 = distreader.as_snp(max_weight=1.0).read()
        s = str(snpdata1)
        output = "tempdir/distreader/toydata.snp.npz"
        create_directory_if_necessary(output)
        SnpNpz.write(output,snpdata1)
        snpreader2 = SnpNpz(output)
        snpdata2 = snpreader2.read()
        np.testing.assert_array_almost_equal(snpdata1.val, snpdata2.val, decimal=10)
        logging.info("done with test")

class TestDistNaNCNC(unittest.TestCase):
    def __init__(self, iid_index_list, snp_index_list, distreader, dtype, order, force_python_only, reference_snps, reference_dtype):
        self.iid_index_list = iid_index_list
        self.snp_index_list = snp_index_list
        self.distreader = distreader
        self.dtype = np.dtype(dtype)
        self.order = order
        self.force_python_only = force_python_only
        self.reference_snps = reference_snps
        self.reference_dtype = reference_dtype

    _testMethodName = "runTest"
    _testMethodDoc = None

    @staticmethod
    def factory_iterator():

        snp_reader_factory_distnpz = lambda : DistNpz("../examples/toydata.dist.npz")
        snp_reader_factory_snpmajor_hdf5 = lambda : DistHdf5("../examples/toydata.snpmajor.dist.hdf5")
        snp_reader_factory_iidmajor_hdf5 = lambda : DistHdf5("../examples/toydata.iidmajor.dist.hdf5")

        previous_wd = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        distreader0 = snp_reader_factory_distnpz()
        S_original = distreader0.sid_count
        N_original = distreader0.iid_count

        snps_to_read_count = min(S_original, 100)

        for iid_index_list in [range(N_original), range(N_original//2), range(N_original - 1,0,-2)]:
            for snp_index_list in [range(snps_to_read_count), range(snps_to_read_count//2), range(snps_to_read_count - 1,0,-2)]:
                reference_snps, reference_dtype = TestDistNaNCNC(iid_index_list, snp_index_list, snp_reader_factory_distnpz(), np.float64, "C", "False", None, None).read_and_standardize()
                for distreader_factory in [snp_reader_factory_distnpz, 
                                            snp_reader_factory_snpmajor_hdf5, snp_reader_factory_iidmajor_hdf5
                                            ]:
                    for dtype in [np.float64,np.float32]:
                        for order in ["C", "F"]:
                            for force_python_only in [False, True]:
                                distreader = distreader_factory()
                                test_case = TestDistNaNCNC(iid_index_list, snp_index_list, distreader, dtype, order, force_python_only, reference_snps, reference_dtype)
                                yield test_case
        os.chdir(previous_wd)

    def __str__(self):
        iid_index_list = self.iid_index_list
        snp_index_list = self.snp_index_list
        distreader = self.distreader
        dtype = self.dtype
        order = self.order
        force_python_only = self.force_python_only
        return "{0}(iid_index_list=[{1}], snp_index_list=[{2}], distreader={3}, dtype={4}, order='{5}', force_python_only=={6})".format(
            self.__class__.__name__,
            ",".join([str(i) for i in iid_index_list]) if len(iid_index_list) < 10 else ",".join([str(i) for i in iid_index_list[0:10]])+",...",
            ",".join([str(i) for i in snp_index_list]) if len(snp_index_list) < 10 else ",".join([str(i) for i in snp_index_list[0:10]])+",...",
            distreader, dtype, order, force_python_only)

    def read_and_standardize(self):
        previous_wd = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        iid_index_list = self.iid_index_list
        snp_index_list = self.snp_index_list
        distreader = self.distreader
        dtype = self.dtype
        order = self.order
        force_python_only = self.force_python_only
        
        snps = distreader[iid_index_list,snp_index_list].read(order=order, dtype=dtype, force_python_only=force_python_only).val
        snps[0,0] = [np.nan]*3 # add a NaN
        snps[:,1] = [.1,.2,.7] # make a SNC
        os.chdir(previous_wd)
        return snps, dtype

    def doCleanups(self):
        pass
        #return super(NaNCNCTestCases, self).doCleanups()

    def runTest(self, result = None):
        snps, dtype = self.read_and_standardize()
        assert not np.array_equal(snps[0,0],snps[0,0]) #without SnpReader's standardization NaN's stay NaN's
        assert np.allclose(snps[0,1],[.1,.2,.7])
        if self.reference_snps is not None:
            self.assertTrue(np.allclose(self.reference_snps, snps, rtol=1e-04 if dtype == np.float32 or self.reference_dtype == np.float32 else 1e-12,equal_nan=True))


# We do it this way instead of using doctest.DocTestSuite because doctest.DocTestSuite requires modules to be pickled, which python doesn't allow.
# We need tests to be pickleable so that they can be run on a cluster.
class TestDistReaderDocStrings(unittest.TestCase):
    def test_distreader(self):
        import pysnptools.distreader.distreader
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(pysnptools.distreader.distreader,optionflags=doctest.ELLIPSIS)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_distdata(self):
        import pysnptools.distreader.distdata
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(pysnptools.distreader.distdata,optionflags=doctest.ELLIPSIS)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__


    def test_disthdf5(self):
        import pysnptools.distreader.disthdf5
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(pysnptools.distreader.disthdf5,optionflags=doctest.ELLIPSIS)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

    def test_distnpz(self):
        import pysnptools.distreader.distnpz
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        result = doctest.testmod(pysnptools.distreader.distnpz,optionflags=doctest.ELLIPSIS)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__


def getTestSuite():
    """
    set up composite test suite
    """

    test_suite = unittest.TestSuite([])

    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestBgen))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestDistReaderDocStrings))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestDistGen))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestDistMemMap))
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestDistReaders))
    test_suite.addTests(TestDistNaNCNC.factory_iterator())
    

    return test_suite

if __name__ == '__main__':
    logging.basicConfig(level=logging.WARN)

    if False:
        from pysnptools.snpreader import Bed
        from pysnptools.distreader import DistData, DistNpz
        # Create toydata.dist.npz
        currentFolder = os.path.dirname(os.path.realpath(__file__))
        if True:
            snpreader = Bed(currentFolder + "/../examples/toydata.5chrom.bed",count_A1=True)[:25,:]
            np.random.seed(392)
            val = np.random.random((snpreader.iid_count,snpreader.sid_count,3))
            val /= val.sum(axis=2,keepdims=True)  #make probabilities sum to 1
            distdata = DistData(iid=snpreader.iid,sid=snpreader.sid,pos=snpreader.pos,val=val)
            DistNpz.write(currentFolder + "/../examples/toydata.dist.npz",distdata)
        if True:
            distdata = DistNpz(currentFolder + "/../examples/toydata.dist.npz").read()
            for sid_major,name_bit in [(False,'iidmajor'),(True,'snpmajor')]:
                DistHdf5.write(currentFolder + "/../examples/toydata.{0}.dist.hdf5".format(name_bit),distdata,sid_major=sid_major)
        if True:
            distdata = DistNpz(currentFolder + "/../examples/toydata.dist.npz")[:,:10].read()
            DistNpz.write(currentFolder + "/../examples/toydata10.dist.npz",distdata)
        if True:
            distdata = DistNpz(currentFolder + "/../examples/toydata.dist.npz")[:,:10].read()
            DistMemMap.write(currentFolder + "/../examples/tiny.dist.memmap",distdata)
        print('done')

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=False)
    ret = r.run(suites)
    assert ret.wasSuccessful()
