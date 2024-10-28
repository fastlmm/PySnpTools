# Download a sample file
from pysnptools.util import example_file
bgen_file = example_file("pysnptools/examples/example.bgen")

# Read from the file
from pysnptools.distreader import Bgen
bgen = Bgen(bgen_file)          # Create a reader
probs0 = bgen[:,0].read().val   # Read 1st SNP
print(probs0.shape)             # Shape of the NumPy array
assert probs0.shape==(500,1,3)

probs_all = bgen.read().val     # Read all variants
print(probs_all.shape)          # Shape of the NumPy array
assert probs_all.shape==(500,199,3)
