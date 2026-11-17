"""Resolve precision consistently across readers, kernels, and constructors."""
import numpy as np
import pytest
import bed_reader
from pysnptools.kernelreader import KernelData, SnpKernel
from pysnptools.pstreader import PstData
from pysnptools.snpreader import SnpData, Bed, Pheno
from pysnptools.standardizer import Identity, Unit
from pysnptools.util import array_dtype, sub_matrix


def source(xp=np):
    return SnpData(iid=[["f", str(sample)] for sample in range(3)], sid=["a", "b"],
                   val=np.array([[0., 1.], [1., 0.], [2., 2.]], dtype=np.float64), xp=xp)


@pytest.mark.parametrize("configured", ["float32", "float64"])
@pytest.mark.parametrize("backend", ["numpy", "cupy"])
def test_read_and_kernel_defaults(monkeypatch, configured, backend):
    xp = np if backend == "numpy" else pytest.importorskip("cupy")
    if backend == "cupy" and not bed_reader.has_gpu_support():
        pytest.skip("requires the cubecl-cuda Cargo feature")
    monkeypatch.setenv("ARRAY_MODULE", backend)
    monkeypatch.setenv("ARRAY_DTYPE", configured)
    snps = source(xp)
    readers = [snps, Bed(bed_reader.sample_file("small.bed"), count_A1=True),
               KernelData(iid=snps.iid, val=np.eye(3), xp=xp), SnpKernel(snps, Identity())]
    for reader in readers:
        result = reader.read(xp=xp)
        assert result.val.dtype == np.dtype(configured)
        assert isinstance(result.val, xp.ndarray)
        assert reader.read(dtype="float64", xp=xp).val.dtype == np.float64
        assert reader.read(dtype="float32", xp=xp).val.dtype == np.float32
    assert snps.read_kernel(Identity(), xp=xp).val.dtype == np.dtype(configured)
    assert SnpKernel(snps, Unit()).read_snps(xp=xp).val.dtype == np.dtype(configured)
    assert sub_matrix(snps.val, [2, 0], [1, 0]).dtype == np.dtype(configured)
    # Existing arrays/standardizers keep their dtype instead of following a new default.
    assert snps.val.dtype == np.float64
    snps.standardize()
    assert snps.val.dtype == np.float64


def test_generic_reader_and_untyped_containers(monkeypatch):
    monkeypatch.setenv("ARRAY_MODULE", "numpy")
    monkeypatch.setenv("ARRAY_DTYPE", "float32")
    generic = PstData(row=["r"], col=["c"], val=np.array([[1.]], dtype=np.float64))
    assert generic.read().val.dtype == np.float32
    assert generic.val.dtype == np.float64
    snps = SnpData(iid=[["f", "0"]], sid=["s"], val=[[1.]])
    assert snps.val.dtype == np.float32
    assert SnpData(iid=[], sid=[], val=None).val.dtype == np.float32
    assert KernelData(iid=[], val=None).val.dtype == np.float32


def test_legacy_dependency_and_invalid_configuration(monkeypatch):
    monkeypatch.delenv("ARRAY_DTYPE", raising=False)
    assert array_dtype() == np.float64
    monkeypatch.delattr(bed_reader, "array_dtype")
    assert array_dtype() == np.float64
    assert array_dtype("float32") == np.float32
    monkeypatch.setenv("ARRAY_DTYPE", "float32")
    with pytest.raises(RuntimeError, match="bed-reader"):
        array_dtype()
    assert array_dtype("float64") == np.float64


def test_pheno_metadata_does_not_allocate_gpu_float64(monkeypatch):
    cupy = pytest.importorskip("cupy")
    monkeypatch.setenv("ARRAY_MODULE", "cupy")
    monkeypatch.setenv("ARRAY_DTYPE", "float32")
    with pytest.MonkeyPatch.context() as patch:
        def unexpected(*args, **kwargs):
            raise AssertionError("metadata inspection must not allocate GPU values")
        patch.setattr(bed_reader, "empty_gpu", unexpected)
        pheno = Pheno({"iid": np.array([["f", "0"]]), "header": ["y"], "vals": np.array([[1e-50]])})
        assert len(pheno.iid) == 1
        assert isinstance(pheno._data.val, np.ndarray)
        assert pheno._data.val[0, 0] == 1e-50
    result = pheno.read(xp=cupy)
    assert result.val.dtype == np.float32


@pytest.mark.parametrize("format_name", ["Dense", "Dat", "Ped"])
def test_text_parsing_keeps_raw_values_on_cpu(tmp_path, monkeypatch, format_name):
    import warnings
    from pysnptools import snpreader

    cupy = pytest.importorskip("cupy")
    reader_type = getattr(snpreader, format_name)
    snps = source(np)
    snps._col_property = np.array([[1, 0, 1], [1, 0, 2]], dtype=np.float64)
    path = str(tmp_path / ("source." + format_name.lower()))
    reader_type.write(path, snps)
    expected = reader_type(path).read(dtype=np.float32, xp=np).val
    monkeypatch.setenv("ARRAY_MODULE", "cupy")
    monkeypatch.setenv("ARRAY_DTYPE", "float32")
    reader = reader_type(path)
    with warnings.catch_warnings():
        warnings.simplefilter("error", bed_reader.GpuFloat64Warning)
        assert reader.iid_count == snps.iid_count
        assert isinstance(reader._data.val, np.ndarray)
        actual = reader.read(xp=cupy)
    assert actual.val.dtype == np.float32
    np.testing.assert_array_equal(cupy.asnumpy(actual.val), expected)


@pytest.mark.parametrize("configured", ["float32", "float64"])
def test_distribution_reader_precision(monkeypatch, configured):
    from pysnptools.distreader import DistData

    monkeypatch.setenv("ARRAY_MODULE", "numpy")
    monkeypatch.setenv("ARRAY_DTYPE", configured)
    values = np.array([[[.5, .5, 0]], [[0, 0, 1]]], dtype=np.float64)
    distributions = DistData(iid=[["f", "0"], ["f", "1"]], sid=["s"], val=values)
    assert distributions.val.dtype == np.float64
    actual = distributions.read()
    assert actual.val.dtype == np.dtype(configured)
    np.testing.assert_array_equal(actual.val, values)
    assert distributions.read(dtype="float64").val.dtype == np.float64
    assert distributions.as_snp().read().val.dtype == np.dtype(configured)
