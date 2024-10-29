"""Tools for reading and manipulating SNP data.
"""

def _snps_fixup(snp_input, iid_if_none=None, count_A1=None):
    """Helper function for SNP data input."""
    import numpy as np
    from pysnptools.snpreader.bed import Bed
    from pysnptools.snpreader.snpdata import SnpData

    if isinstance(snp_input, str):
        return Bed(snp_input, count_A1=count_A1)

    if isinstance(snp_input, dict):
        return SnpData(iid=snp_input['iid'], sid=snp_input['header'], val=snp_input['vals'])

    if snp_input is None:
        assert iid_if_none is not None, "snp_input cannot be None here"
        return SnpData(
            iid_if_none, sid=np.empty((0), dtype="str"),
            val=np.empty((len(iid_if_none), 0)), pos=np.empty((0, 3)), name=""
        )

    return snp_input

# Import other necessary components, with noqa comments to ignore unused import warnings
from pysnptools.snpreader.snpreader import SnpReader  # noqa: E402, F401
from pysnptools.snpreader.snpdata import SnpData  # noqa: E402, F401
from pysnptools.snpreader.bed import Bed, plink_chrom_map, reverse_plink_chrom_map  # noqa: E402, F401
from pysnptools.snpreader.ped import Ped  # noqa: E402, F401
from pysnptools.snpreader.dat import Dat  # noqa: E402, F401
from pysnptools.snpreader.snphdf5 import SnpHdf5, Hdf5  # noqa: E402, F401
from pysnptools.snpreader.snpnpz import SnpNpz  # noqa: E402, F401
from pysnptools.snpreader.dense import Dense  # noqa: E402, F401
from pysnptools.snpreader.pheno import Pheno  # noqa: E402, F401
from pysnptools.snpreader.snpmemmap import SnpMemMap  # noqa: E402, F401
from pysnptools.snpreader._mergesids import _MergeSIDs  # noqa: E402, F401
from pysnptools.snpreader._mergeiids import _MergeIIDs  # noqa: E402, F401
from pysnptools.snpreader.snpgen import SnpGen  # noqa: E402, F401
from pysnptools.snpreader.distributedbed import DistributedBed, _Distributed1Bed  # noqa: E402, F401

__all__ = [
    "SnpReader",
    "SnpData",
    "Bed",
    "plink_chrom_map",
    "reverse_plink_chrom_map",
    "Ped",
    "Dat",
    "SnpHdf5",
    "Hdf5",
    "SnpNpz",
    "Dense",
    "Pheno",
    "SnpMemMap",
    "_MergeSIDs",
    "_MergeIIDs",
    "SnpGen",
    "DistributedBed",
    "_Distributed1Bed",
]
