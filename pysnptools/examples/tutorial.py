# Need to upgrade to Python 3
# import os
# import numpy as np

# def randlong(l): #Because of round off error will not get some at end
#     sr = int(np.sqrt(l))
#     result = np.random.randint(sr) + sr * np.random.randint(sr)
#     assert 0 <= result < l
#     return result

# def region_gen(scale=1,seed=0):
#     range_count = int(150000*scale)
#     position_count = int(3000000000*scale) #3 billion
#     region_max_length = int(2000000*scale) #2 million

#     np.random.seed(seed)
#     for range_index in range(range_count):
#         length = int(np.exp(np.random.random()*np.log(region_max_length)))
#         start = randlong(position_count-length) #does randint really go up to 3 billon?
#         stop = start+length
#         yield start,stop

# from pysnptools.util import IntRangeSet

# geneset = IntRangeSet()
# for start,stop in region_gen(scale=.1):
#     geneset |= (start,stop)
# print geneset
# print geneset.ranges_len

# print("done")



# os.chdir(r"C:\Source\carlk\fastlmm2\tests\datasets\synth")

# from pysnptools.snpreader import Bed

# # Use "Bed" to access file "all.bed"
# snpreader = Bed("all.bed")

# # What is snpreader?
# print snpreader
# # Bed("all.bed")


# # Find out about iids and sids
# print snpreader.iid_count
# print snpreader.sid_count
# print snpreader.iid[:3]
# print snpreader.sid[:3]
# #500
# #5000
# #[['cid0P0' 'cid0P0']
# # ['cid1P0' 'cid1P0']
# # ['cid2P0' 'cid2P0']]
# #['snp625_m0_.03m1_.07' 'snp1750_m0_.02m1_.04' 'snp0_m0_.37m1_.24']

# #Read all the SNP data in to memory
# snpdata = snpreader.read()
# #What is snpdata?
# # SnpData(Bed("all.bed"))

# #What do the iids and sid of snpdata look like?
# print snpdata.iid_count, snpdata.sid_count
# print snpdata.iid[:3]
# print snpdata.sid[:3]
# # The same.

# # print the SNP data
# print snpdata.val
# #[[ 2.  2.  1. ...,  2.  1.  2.]
# # [ 2.  2.  1. ...,  2.  0.  2.]
# # [ 2.  2.  1. ...,  1.  1.  1.]
# # ...,
# # [ 2.  2.  2. ...,  1.  2.  2.]
# # [ 2.  2.  2. ...,  1.  2.  2.]
# # [ 2.  2.  2. ...,  2.  0.  2.]]

# # snpdata.val is a NumPy array. Can apply any np functions
# print np.mean(snpdata.val)
# #1.478588

# #If all you want is to read data in a Numpy array, here it is one line:
# print np.mean(Bed("all.bed").read().val)

# #You can also create a SnpData object from scratch (without reading from a SnpReader)
# from pysnptools.snpreader import SnpData
# snpdata1 = SnpData(iid=[['f1','c1'],['f1','c2'],['f2','c1']],sid=['snp1','snp2'],val=[[0,1],[2,.5],[.5,np.nan]])
# print np.nanmean(snpdata1.val)
# # 0.8



# #Review SnpReader and Bed and SnpData, and common attributes including val

# #Topics: Reading subsets of data, reading with re-ordering iids & sids (rows & cols), stacking

# #Reading just one snp
# snpreader = Bed("all.bed")
# snp0reader = snpreader[:,0]
# print snp0reader, snp0reader.iid_count, snp0reader.sid_count, snp0reader.sid
# # Bed('all.bed')[:,0] 500 1 ['snp625_m0_.03m1_.07']
# print snpreader
# # Bed("all.bed")
# snp0data = snp0reader.read()
# print snp0data, snp0data.iid_count, snp0data.sid_count, snp0data.sid
# #SnpData(Bed('all.bed')[:,0]) 500 1 ['snp625_m0_.03m1_.07']
# print snp0data.val
# #[[ 2.]
# # [ 2.]
# # [ 2.]
# # [ 2.]
# # [ 2.]
# # [ 2.]
# # [ 2.]
# # [ 2.]
# #....]]

# # print the data for iid #9
# print Bed("all.bed")[9,:].read().val
# # [[ 2.  2.  1. ...,  1.  2.  1.]]

# # Read the data for the first 5 iids AND the first 5 sids:
# snp55data = Bed("all.bed")[:5,:5].read()
# print snp55data, snp55data.iid_count, snp55data.iid, snp55data.sid_count, snp55data.sid
# #SnpData(Bed('all.bed')[:5,:5]) 5 [['cid0P0' 'cid0P0']
# # ['cid1P0' 'cid1P0']
# # ['cid2P0' 'cid2P0']
# # ['cid3P0' 'cid3P0']
# # ['cid4P0' 'cid4P0']] 5 ['snp625_m0_.03m1_.07' 'snp1750_m0_.02m1_.04' 'snp0_m0_.37m1_.24'
# # 'snp375_m0_.52m1_.68' 'snp1125_m0_.26m1_.27']
# print snp55data.val
# #[[ 2.  2.  1.  0.  2.]
# # [ 2.  2.  1.  1.  2.]
# # [ 2.  2.  1.  0.  1.]
# # [ 2.  2.  2.  0.  2.]
# # [ 2.  2.  2.  2.  2.]]

# #Stacked is OK and efficient
# snpreaderA = Bed("all.bed") # read all
# snpreaderB = snpreaderA[:,:250] #read first sids
# snpreaderC = snpreaderB[:10,:] # reader first 10 iids
# # recall NumPy slice notation:   start:stop:step, so ::2 is every other
# snpreaderD = snpreaderC[::2,::2]
# print snpreaderD, snpreaderD.iid_count, snpreaderD.iid, snpreaderD.sid_count
# #Bed('all.bed')[:,:250][:10,:][::2,::2] 5 [['cid0P0' 'cid0P0']
# # ['cid2P0' 'cid2P0']
# # ['cid4P0' 'cid4P0']
# # ['cid6P0' 'cid6P0']
# # ['cid8P0' 'cid8P0']] 125
# print snpreaderD.read().val
# #[[ 2.  1.  2.  0.  1....
# #   2.  2.  2.  1.  1.
# #   2.  1.  1.  0.  1.
# #   2.  2.  2.  2.  2.
# #   2.  1.  2.  2.  2.
# #   2.  2.  2.  1.  2.
# #   2.  2.  1.  0.  1.
# # [ 2.  1.  1.  2.  0.
# #   1.  2.  2.  2.  2.
# #   2.  1.  0.  0.  1.
# #   1.  1.  1.  2.  2.
# #....

# #Fancy indexing - list of indexes, slices, list of booleans, negatives
# # on iid or sid or both

# # List of indexes (can use to reorder)
# snpdata43210 = Bed("all.bed")[[4,3,2,1,0],:].read()
# print snpdata43210.iid
# #[['cid4P0' 'cid4P0']
# # ['cid3P0' 'cid3P0']
# # ['cid2P0' 'cid2P0']
# # ['cid1P0' 'cid1P0']
# # ['cid0P0' 'cid0P0']]
# # List of booleans to select
# snpdata43210B = snpdata43210[[False,True,True,False,False],:]
# print  snpdata43210B, snpdata43210B.iid
# #SnpData(Bed('all.bed')[[4,3,2,1,0],:])[[1,2],:] [['cid3P0' 'cid3P0']
# # ['cid2P0' 'cid2P0']]

# #Question: Does snpdata43210B have a val property?
# #Answer: No. It's a subset of a SnpData, so it reads from a SnpData, but it is not a SnpData.
# #Use .read() to get its values
# print snpdata43210B.read(view_ok=True)

# #Negatives -- missing some support
# # 'start','stop': negative means counting from the end
# # 'step': negative means count backwards
# print range(10)
# #[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
# print range(10)[-4:]
# # [6, 7, 8, 9]
# print range(10)[::-1]
# #[9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
# print range(10)[-2:-7:-2]
# # [8, 6, 4]

# print Bed("all.bed")[::-10,:].iid #Slices work
# #[['cid499P1' 'cid499P1']
# # ['cid489P1' 'cid489P1']
# # ['cid479P1' 'cid479P1']
# # ...
# # ['cid19P0' 'cid19P0']
# # ['cid9P0' 'cid9P0']]

# # Indexing at negative or list of negative, not yet supported
# # print Bed("all.bed")[-1,-1].read().val



# #Summary:
# # You can use indexing (fancy!) on a reader to specify a reader for a subset of the data
# # Subset reader can also reorder and sids
# # No SNP data is actually read until you say ".read()"
# # Every read gives fresh data (possible exception: view_ok=True)
# # The result of reading is a SnpData which is a SnpReader with a .val attribute


# #Topic more properties and attributes of SnpReaders

# #read() supports both memory layouts and 8 byte or 4 byte floats
# print Bed("all.bed").read().val.flags
#   #C_CONTIGUOUS : False
#   #F_CONTIGUOUS : True
#   #OWNDATA : True
#   #WRITEABLE : True
#   #ALIGNED : True
#   #UPDATEIFCOPY : False
# snpdata32c = Bed("all.bed").read(order='C',dtype=np.float32)
# print snpdata32c.val.dtype
# #float32
# print snpdata32c.val.flags
#   #C_CONTIGUOUS : True
#   #F_CONTIGUOUS : False
#   #OWNDATA : True
#   #WRITEABLE : True
#   #ALIGNED : True
#   #UPDATEIFCOPY : False

# # Every reader includes an array of SNP properties called ".pos"
# print Bed("all.bed").pos
# #[[   1    0    0]
# # [   1    1    1]
# # [   1    2    2]
# # ...,
# # [   5 4997 4997]
# # [   5 4998 4998]
# # [   5 4999 4999]]

# # chromosome, genetic distance, basepair distance
# # Accessible without a SNP data read.

# # So, using Python fancy indexing, how to we read all SNPs at Chrom 5?
# snpreader = Bed("all.bed")
# chr5_bools = (snpreader.pos[:,0] == 5)
# print chr5_bools
# # array([False, False, False, ...,  True,  True,  True], dtype=bool)
# chr5reader = snpreader[:,chr5_bools]
# print chr5reader
# #Bed('all.bed')[:,[4000,4001,4002,4003,4004,4005,4006,4007,4008,4009,...]]
# chr5data = chr5reader.read()
# print chr5data.pos
# #[[   5 4000 4000]
# # [   5 4001 4001]
# # [   5 4002 4002]
# # ...,
# # [   5 4997 4997]
# # [   5 4998 4998]
# # [   5 4999 4999]]

# #In one-line:
# chr5data = Bed("all.bed")[:,snpreader.pos[:,0] == 5].read()

# # You can turn iid or sid names into indexes
# snpreader = Bed("all.bed")
# iid0 =[['cid499P1','cid499P1'],
#       ['cid489P1','cid489P1'],
#       ['cid479P1','cid479P1']]
# indexes0 = snpreader.iid_to_index(iid0)
# print indexes0
# #array([499, 489, 479])
# snpreader0 = snpreader[indexes0,:]
# print snpreader0.iid
# #[['cid499P1' 'cid499P1']
# # ['cid489P1' 'cid489P1']
# # ['cid479P1' 'cid479P1']]

# # more condensed
# snpreader0 = snpreader[snpreader.iid_to_index(iid0),:]

# #both a once
# snpdata0chr5 = snpreader[snpreader.iid_to_index(iid0),snpreader.pos[:,0] == 5].read()
# print np.mean(snpdata0chr5.val)
# # 1.493

# #Summary:
# # Every SnpReader has a .pos attribute
# #      a .snpcount x 3 array of genetic info
# #      available before any reads
# # The iid_to_index and sid_to_index methods turn iid's and sid's into indexes
# # Via NumPy-style indexing, these allow reading by name and genetic property

# #Topic: Other SnpReaders and how to write

# #Read from the PLINK phenotype file (text) instead of a Bed file
# # Looks like:
# #cid0P0 cid0P0 0.4853395139922632
# #cid1P0 cid1P0 -0.2076984565752155
# #cid2P0 cid2P0 1.4909084058931985
# #cid3P0 cid3P0 -1.2128996652683697
# #cid4P0 cid4P0 0.4293203431508744
# #...

# from pysnptools.snpreader import Pheno
# phenoreader = Pheno("pheno_10_causals.txt")
# print phenoreader, phenoreader.iid_count, phenoreader.sid_count, phenoreader.sid, phenoreader.pos
# #Pheno('pheno_10_causals.txt') 500 1 ['pheno0'] [[ nan  nan  nan]]
# phenodata = phenoreader.read()
# print phenodata.val
# #[[  4.85339514e-01]
# # [ -2.07698457e-01]
# # [  1.49090841e+00]
# # [ -1.21289967e+00]
# # ...

# # Write 1st 10 iids and sids of Bed data into Pheno format
# snpdata1010 = Bed("all.bed")[:10,:10].read()
# Pheno.write("deleteme1010.txt",snpdata1010)

# #Write it to Bed format
# Bed.write("deleteme1010.bed",snpdata1010)

# # Create a snpdata on the fly and write to Bed
# snpdata1 = SnpData(iid=[['f1','c1'],['f1','c2'],['f2','c1']],sid=['snp1','snp2'],val=[[0,1],[2,1],[1,np.nan]])
# Bed.write("deleteme1.bed",snpdata1)


# #Pheno is slow because its txt. Bed format can only hold 0,1,2,missing.
# # Use SnpNpz for fastest read/write times, smallest file size
# from pysnptools.snpreader import SnpNpz
# SnpNpz.write("deleteme1010.snp.npz", snpdata1010)

# # Use SnpHdf5 for random-access reads, good speed and size, and compatibility outside Python
# from pysnptools.snpreader import SnpHdf5
# SnpHdf5.write("deleteme1010.snp.hdf5", snpdata1010)

# #Summary: Every format has its own SnpReader class
# #       Table: Pheno, SnpNpz, SnpHdf5
# #   That SnpReader has a static write method for SnpData


# #Topics: Intersecting iids
# #What if we have two data sources with slightly different iids in different order?
# snpreader = Bed("all.bed")
# phenoreader = Pheno("pheno_10_causals.txt")[::-2,:]
# print snpreader.iid_count, phenoreader.iid_count, snpreader.iid, phenoreader.iid
# #Create an intersecting and reordering reader with
# import pysnptools.util as pstutil
# snpreader_i,phenoreader_i  = pstutil.intersect_apply([snpreader,phenoreader])
# assert np.array_equal(snpreader_i.iid,phenoreader_i.iid)
# snpdata_i = snpreader_i.read()
# phenodata_i = phenoreader_i.read()

# bs = np.linalg.lstsq(snpdata_i.val, phenodata_i.val,rcond=-1)[0] #usually would add a 1's column
# predicted = snpdata_i.val.dot(bs)
# import matplotlib.pyplot as plt
# plt.plot(phenodata_i.val, predicted, '.', markersize=10)
# #plt.show() #Easy to 'predict' seen 250 cases with 5000 variables.
# # How does it predict unseen cases?
# phenoreader_unseen = Pheno("pheno_10_causals.txt")[-2::-2,:]
# snpreader_u,phenoreader_u  = pstutil.intersect_apply([snpreader,phenoreader_unseen])
# snpdata_u = snpreader_u.read()
# phenodata_u = phenoreader_u.read()
# predicted_u = snpdata_u.val.dot(bs)
# plt.plot(phenodata_u.val, predicted_u, '.', markersize=10)
# #plt.show() #Hard to predict unseen 250 cases with 5000 variables.

# #Summary:
# # pstutil.intersect_apply can insect and reorder readers by iid
# # Because it works before reading SNP data from disk could be more efficient
# # Also works with None, and tuples of (iid list,numpy array)
# # Can stack more subsetting (e.g. SNP selection) before the reads


# #Topics: Standardization, Kernels

# # To Unit standardize: read data, ".standardize()"

# snpreader = Bed("all.bed")
# snpdata = snpreader.read()
# snpdata = snpdata.standardize() #In place AND returns self
# print snpdata.val
# #[[ 0.30156099  0.2481353  -0.50673344 ...,  0.92208184 -0.1266665   0.55601103]
# # [ 0.30156099  0.2481353  -0.50673344 ...,  0.92208184 -1.5034763   0.55601103]
# #...

# # In one-line:
# snpdata = Bed("all.bed").read().standardize()

# # Beta standardization
# from pysnptools.standardizer import Beta
# snpdataB = Bed("all.bed").read().standardize(Beta(1,25))
# print snpdataB.val
# #[[  7.40112054e-01   7.15532756e-01  -5.02003205e-04 ...,   4.40649336e-03   -1.13331663e-06   1.87525732e-01]
# # [  7.40112054e-01   7.15532756e-01  -5.02003205e-04 ...,   4.40649336e-03   -1.34519756e-05   1.87525732e-01]
# # ...

# # To create an kernel (the relatedness of each iid pair as the dot product of their standardized SNP values)
# from pysnptools.standardizer import Unit
# kerneldata = Bed("all.bed").read_kernel(standardizer=Unit())
# print kerneldata.val
# #array([[ 5081.6121922 ,   253.32922313,   165.9842232 , ...,  -130.76998392,  -298.66392286,  -287.66887036],
# #       [  253.32922313,  5061.87849635,   384.04149913, ...,  -334.33599388,  -127.02308706,  -291.41483161]
# #
# #...

# # Low memory:
# kerneldata = Bed("all.bed").read_kernel(standardizer=Unit(),block_size=500)



# # Summary
# #  Standardization
# #     default Unit - mean 0, stdev 1, THEN fill with 0
# #     In place, and returns self
# #     Other standardizers: Beta, Unit, DiagKToN
# #  Kernels
# #     Returns a kerneldata something new with a .val property
# #     Does the read itself, allowing memory to be saved
# #     Requires a standardizer. Use Identity() for none


# print("done!")
