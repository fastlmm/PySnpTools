"""GPU API regressions. Requires a working CuPy GPU."""
import numpy as np
import pytest
from unittest.mock import patch

from bed_reader import sample_file
from pysnptools.snpreader import Bed, SnpData
from pysnptools.kernelreader import SnpKernel
from pysnptools.standardizer import Unit, Beta


@pytest.fixture
def cupy():
    # Optional GPU qualification is explicitly reported on CPU-only installations.
    from pysnptools.util.gpu import has_gpu_support
    if not has_gpu_support():
        pytest.skip("requires CuPy and a usable CUDA device")
    return pytest.importorskip("cupy")


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("standardizer", [Unit(), Beta(1, 25)])
def test_subset_and_blocked(cupy, dtype, standardizer):
    reader = Bed(sample_file("toydata.5chrom.bed"), count_A1=True)[::-3, ::97][::-1, ::-1]
    actual, gpu_trained = reader._read_with_standardizing(
        standardizer, dtype=dtype, xp=cupy, return_trained=True)
    # CPU read and transfer are exact, before the CuPy standardizer is applied.
    np.testing.assert_array_equal(cupy.asnumpy(reader.read(dtype=dtype, xp=cupy).val),
                                  reader.read(dtype=dtype, xp=np).val)
    # The removed native-dispatch test required bitwise CPU/native parity.
    # Verify the replacement CPU-read/transfer/CuPy pipeline exactly; independent
    # CPU/reference accuracy bounds are covered in the general standardizer tests.
    explicit, explicit_trained = reader.read(dtype=dtype, xp=cupy).standardize(
        standardizer, return_trained=True)
    np.testing.assert_array_equal(cupy.asnumpy(actual.val), cupy.asnumpy(explicit.val))
    np.testing.assert_array_equal(gpu_trained.stats, explicit_trained.stats)
    kernel = reader.read_kernel(standardizer, block_size=7, dtype=dtype, xp=cupy)
    # CPU Python and CuPy both use SciPy Beta PDFs. Compare their independent
    # reductions/products at the original blocked-kernel tolerance; CPU Rust's
    # statrs PDF has a known different normalization in the last bits.
    expected = reader.read(dtype=dtype, xp=np).standardize(
        standardizer, force_python_only=True)
    tolerance = 32 * np.finfo(dtype).eps
    np.testing.assert_allclose(cupy.asnumpy(kernel.val), expected.val @ expected.val.T,
                               atol=tolerance * reader.sid_count, rtol=tolerance)
    for to_kerneldata in [False, True]:
        result = SnpKernel(reader, standardizer, block_size=7)._read_with_standardizing(
            to_kerneldata, dtype=dtype, xp=cupy)
        val = result.val if to_kerneldata else result.snpreader.val
        assert val.dtype == dtype and isinstance(val, cupy.ndarray)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("order", ["C", "F"])
@pytest.mark.parametrize("standardizer", [Beta(1, 25), Beta(0.5, 0.5)])
def test_blocked_beta_trained_bed(cupy, tmp_path, dtype, order, standardizer):
    from bed_reader import to_bed

    rng = np.random.default_rng(37)
    values = rng.integers(0, 3, size=(13, 9)).astype(dtype)
    values[:8, :2] = [0, 2]  # training constants can vary in the test cohort
    values[8:, 4] = np.nan  # applying moments must not refit an empty test SNP
    path = tmp_path / "trained.bed"
    to_bed(path, values, count_A1=True)
    reader = Bed(path, count_A1=True)
    _, trained = reader[:8, :].read(dtype=dtype, xp=np).standardize(standardizer, return_trained=True)
    test_reader = reader[8:, :][::-1, :]
    standardized = test_reader.read(dtype=dtype, xp=np).standardize(trained).val
    expected = standardized @ standardized.T
    # Nine well-scaled SNP contributions and three block additions round in dtype.
    tolerance = 64 * np.finfo(dtype).eps
    for xp in (np, cupy):
        kernel = test_reader.read_kernel(trained, block_size=4, dtype=dtype, order=order, xp=xp)
        assert isinstance(kernel.val, xp.ndarray) and kernel.val.dtype == dtype
        assert kernel.val.flags["C_CONTIGUOUS" if order == "C" else "F_CONTIGUOUS"]
        actual = kernel.val if xp is np else cupy.asnumpy(kernel.val)
        np.testing.assert_allclose(actual, expected, rtol=tolerance, atol=tolerance * reader.sid_count)


def test_explicit_module_and_existing_array(cupy):
    reader = Bed(sample_file("small.bed"), count_A1=True)
    with patch.dict("os.environ", {"ARRAY_MODULE": "cupy"}):
        assert isinstance(reader.read(xp=np).val, np.ndarray)
        assert isinstance(reader.read().val, cupy.ndarray)
        data = reader.read(xp=np)
        data.standardize()
        assert isinstance(data.val, np.ndarray)
        assert isinstance(reader.read_kernel(Unit(), xp=np).val, np.ndarray)
    with patch.dict("os.environ", {"ARRAY_MODULE": "numpy"}):
        data = reader.read(xp=cupy)
        data.standardize()
        assert isinstance(data.val, cupy.ndarray)
        assert isinstance(data[::2, :].read(xp=cupy).val, cupy.ndarray)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("order", ["F", "C"])
def test_foreign_dosages_and_trained(cupy, dtype, order):
    rng = np.random.default_rng(1)
    val = rng.uniform(0.0, 2.0, (41, 6)).astype(dtype)
    val[::5, 0] = np.nan
    val[:, 1] = 1.0
    iid = np.array([[str(index), str(index)] for index in range(len(val))])
    sid = np.array([str(index) for index in range(val.shape[1])])
    for standardizer in [Unit(), Beta(1, 25)]:
        # Beta is only defined for allele frequencies in [0, 1].
        source = val if isinstance(standardizer, Unit) else val[:, :2]
        data_np = SnpData(iid=iid, sid=sid[:source.shape[1]], val=np.array(source, order=order), xp=np)
        data_cp = SnpData(iid=iid, sid=data_np.sid, val=cupy.array(source, order=order), xp=cupy)
        expected, trained_np = data_np.standardize(standardizer, return_trained=True)
        actual, trained_cp = data_cp.standardize(standardizer, return_trained=True)
        tolerance = 64 * np.finfo(dtype).eps
        np.testing.assert_allclose(cupy.asnumpy(actual.val), expected.val, rtol=tolerance, atol=tolerance)
        fresh_np = SnpData(iid=iid, sid=data_np.sid, val=np.array(source, order=order), xp=np).standardize(trained_np)
        fresh_cp = SnpData(iid=iid, sid=data_np.sid, val=cupy.array(source, order=order), xp=cupy).standardize(trained_cp)
        np.testing.assert_allclose(cupy.asnumpy(fresh_cp.val), fresh_np.val, rtol=tolerance, atol=tolerance)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_large_offset_against_centered_reference(cupy, dtype):
    # The legacy Rust CPU standardizer uses E[x*x]-E[x]**2 and loses
    # precision here, so use an independent centered float64 reference.
    rng = np.random.default_rng(12)
    values = (1e4 + rng.random((100, 4))).astype(dtype)
    mean = values.astype(np.float64).mean(axis=0).astype(dtype)
    deviation = values.astype(np.float64).std(axis=0).astype(dtype)
    expected = (values - mean) * (1 / deviation)
    data = SnpData(iid=[[str(index), str(index)] for index in range(100)],
                   sid=["a", "b", "c", "d"], val=cupy.array(values, order="F"), xp=cupy)
    actual = data.standardize(Unit()).val
    # Parallel reductions can round the large mean by a few ulps. Subtracting
    # it and dividing by a small deviation amplifies that error by this ratio.
    absolute_tolerance = 8 * np.finfo(dtype).eps * np.max(np.abs(values)) / np.min(deviation)
    np.testing.assert_allclose(cupy.asnumpy(actual), expected,
                               rtol=8 * np.finfo(dtype).eps, atol=absolute_tolerance)


@pytest.mark.parametrize("order", ["C", "F"])
def test_standardize_workspace_is_not_matrix_sized(cupy, order):
    rng = np.random.default_rng(8)
    values = rng.integers(0, 3, size=(1000, 5000), dtype=np.int8)
    pool = cupy.cuda.MemoryPool()
    with cupy.cuda.using_allocator(pool.malloc):
        data = SnpData(
            iid=[[str(index), str(index)] for index in range(values.shape[0])],
            sid=[str(index) for index in range(values.shape[1])],
            val=cupy.array(values, dtype=np.float32, order=order), xp=cupy)
        cupy.cuda.get_current_stream().synchronize()
        baseline = pool.total_bytes()
        data.standardize(Unit())
        cupy.cuda.get_current_stream().synchronize()
        # The 19 MiB matrix must not need a matrix-sized mask/temporary.
        # Four MiB allows reduction scratch space and several per-SNP vectors.
        assert pool.total_bytes() - baseline < 4 * 2**20
        assert bool(cupy.isfinite(data.val).all())


def _non_bed_reader(kind, tmp_path):
    from pysnptools.distreader import DistNpz
    from pysnptools.snpreader import SnpMemMap, SnpNpz
    from pysnptools.util import example_file

    if kind == "dist_as_snp":
        return DistNpz(example_file("pysnptools/examples/toydata.dist.npz")).as_snp()[:, :40]
    source = Bed(sample_file("small.bed"), count_A1=True).read(dtype="float32", xp=np)
    if kind == "npz":
        return SnpNpz.write(str(tmp_path / "small.snp.npz"), source)
    if kind == "memmap":
        return SnpMemMap.write(str(tmp_path / "small.snp.memmap"), source)
    raise ValueError(f"unknown reader kind {kind!r}")


@pytest.mark.parametrize("kind", ["npz", "memmap", "dist_as_snp"])
@pytest.mark.parametrize("order", ["C", "F"])
def test_non_bed_readers_decode_on_cpu_then_use_cupy(cupy, tmp_path, kind, order):
    # Every reader shares SnpReader._read_array: CPU decoding, then transfer.
    reader = _non_bed_reader(kind, tmp_path)
    rows = list(range(reader.iid_count - 1, -1, -2))
    subset = reader[rows, 1:]
    gpu = subset.read(dtype="float32", order=order, xp=cupy)
    assert isinstance(gpu.val, cupy.ndarray) and gpu.val.dtype == np.float32
    assert gpu.val.flags["C_CONTIGUOUS" if order == "C" else "F_CONTIGUOUS"]
    np.testing.assert_array_equal(cupy.asnumpy(gpu.val), subset.read(dtype="float32", xp=np).val)

    cpu_kernel = reader.read_kernel(Unit(), block_size=2, dtype="float64", xp=np).val
    gpu_kernel = reader.read_kernel(Unit(), block_size=2, dtype="float64", xp=cupy).val
    assert isinstance(gpu_kernel, cupy.ndarray)
    # Float64 blocked accumulation; only the reduction order differs.
    np.testing.assert_allclose(cupy.asnumpy(gpu_kernel), cpu_kernel,
                               rtol=1e-12, atol=1e-12 * np.abs(cpu_kernel).max())


def _containers(cupy, xp):
    """Return a factory for each in-memory container, holding ``xp`` values."""
    from pysnptools.distreader import DistData
    from pysnptools.kernelreader import KernelData
    from pysnptools.pstreader import PstData

    iid = [["f", "0"], ["f", "1"]]
    snp_values = xp.asarray(np.array([[0.0, 1.0], [2.0, 1.0]], dtype=np.float32))
    kernel_values = xp.asarray(np.array([[1.0, 0.5], [0.5, 1.0]], dtype=np.float32))
    dist_values = xp.asarray(np.full((2, 2, 3), 1 / 3, dtype=np.float32))
    return {
        "SnpData": lambda **kw: SnpData(iid=iid, sid=["s0", "s1"], val=snp_values, **kw),
        "KernelData": lambda **kw: KernelData(iid=iid, val=kernel_values, **kw),
        "PstData": lambda **kw: PstData(row=["r0", "r1"], col=["c0", "c1"], val=snp_values, **kw),
        "DistData": lambda **kw: DistData(iid=iid, sid=["s0", "s1"], val=dist_values, **kw),
    }


@pytest.mark.parametrize("environment", [None, "numpy", "cupy"])
def test_implicit_device_transfer_raises(cupy, environment):
    # xp=None means ARRAY_MODULE, then NumPy. An existing NumPy or CuPy array
    # from the other module must not be copied between devices implicitly.
    import os

    with patch.dict("os.environ"):
        os.environ.pop("ARRAY_MODULE", None)
        if environment is not None:
            os.environ["ARRAY_MODULE"] = environment
        default = cupy if environment == "cupy" else np
        other = np if default is cupy else cupy
        for name, factory in _containers(cupy, other).items():
            with pytest.raises(ValueError, match=f"{name} val is a {other.__name__} array, "
                                                 f"but xp=None selects {default.__name__}"):
                factory()
            for xp in (np, cupy):  # an explicit module keeps or copies the values
                data = factory(xp=xp)
                assert isinstance(data.val, xp.ndarray), (name, xp)
            # Reads from in-memory values follow the same rule.
            source = factory(xp=other)
            with pytest.raises(ValueError, match="xp=None selects"):
                source.read()
            with pytest.raises(ValueError, match="xp=None selects"):
                source[:1, :].read()
            for xp in (np, cupy):
                assert isinstance(source.read(xp=xp).val, xp.ndarray), (name, xp)
                assert isinstance(source[:1, :].read(xp=xp).val, xp.ndarray), (name, xp)
            # Same-module values are untouched by the default.
            same = factory(xp=default)
            assert isinstance(same.read().val, default.ndarray), name


def test_implicit_kernel_reads_raise(cupy, monkeypatch):
    monkeypatch.delenv("ARRAY_MODULE", raising=False)  # these cases assume the NumPy default
    from pysnptools.kernelreader import SnpKernel

    iid = [["f", str(row)] for row in range(4)]
    values = np.array([[0, 1, 2], [1, 1, 0], [2, 0, 1], [0, 2, 2]], dtype=np.float32)
    gpu = SnpData(iid=iid, sid=["s0", "s1", "s2"], val=cupy.asarray(values), xp=cupy)
    for block_size in (None, 1):
        with pytest.raises(ValueError, match="SnpData val is a cupy array"):
            gpu.read_kernel(Unit(), block_size=block_size)
        with pytest.raises(ValueError, match="SnpData val is a cupy array"):
            SnpKernel(gpu, Unit(), block_size=block_size).read()
        # An explicit outer request covers the internal block reads it makes.
        for xp in (np, cupy):
            kernel = gpu.read_kernel(Unit(), block_size=block_size, xp=xp)
            assert isinstance(kernel.val, xp.ndarray)
            snpkernel = SnpKernel(gpu, Unit(), block_size=block_size).read(xp=xp)
            assert isinstance(snpkernel.val, xp.ndarray)
        cpu_kernel = gpu.read_kernel(Unit(), block_size=block_size, xp=np).val
        np.testing.assert_allclose(cupy.asnumpy(gpu.read_kernel(Unit(), block_size=block_size, xp=cupy).val),
                                   cpu_kernel, rtol=1e-5, atol=1e-5)


def test_val_assignment_keeps_module(cupy, monkeypatch):
    monkeypatch.delenv("ARRAY_MODULE", raising=False)  # these cases assume the NumPy default
    for name, factory in _containers(cupy, np).items():
        for xp, other in ((np, cupy), (cupy, np)):
            data = factory(xp=xp)
            replacement = other.asarray(np.asarray(factory(xp=np).val))
            with pytest.raises(ValueError, match=f"{name} val is a {xp.__name__} array; "
                                                 f"cannot assign a {other.__name__} array"):
                data.val = replacement
            data.val = xp.asarray(np.asarray(factory(xp=np).val))  # same module is fine
            assert isinstance(data.val, xp.ndarray)


def test_files_and_memory_maps_follow_the_default(cupy, monkeypatch, tmp_path):
    monkeypatch.delenv("ARRAY_MODULE", raising=False)  # these cases assume the NumPy default
    from pysnptools.snpreader import SnpMemMap

    # Values decoded from a file have no device, so ARRAY_MODULE decides.
    reader = Bed(sample_file("small.bed"), count_A1=True)
    memmap = SnpMemMap.write(str(tmp_path / "small.snp.memmap"), reader.read(xp=np))
    with patch.dict("os.environ", {"ARRAY_MODULE": "cupy"}):
        assert isinstance(reader.read().val, cupy.ndarray)
        reopened = SnpMemMap(str(tmp_path / "small.snp.memmap"))
        assert isinstance(reopened.val, np.ndarray)  # a memory map stays NumPy
        assert isinstance(reopened.read().val, cupy.ndarray)
    assert isinstance(memmap.read().val, np.ndarray)


@pytest.mark.parametrize("kind", ["Dat", "Dense", "Ped", "Pheno"])
def test_text_readers_follow_array_module(cupy, tmp_path, kind):
    # Text readers parse the whole file into a cached NumPy SnpData. That cache
    # holds file values, not a caller's array, so ARRAY_MODULE=cupy may copy it.
    from pysnptools import snpreader

    data = Bed(sample_file("small.bed"), count_A1=True).read(xp=np)
    if kind == "Pheno":
        data = data[:, :1].read(xp=np)
    reader_class = getattr(snpreader, kind)
    reader = reader_class.write(str(tmp_path / f"small.{kind.lower()}"), data)
    with patch.dict("os.environ", {"ARRAY_MODULE": "cupy"}):
        for result in (reader.read(), reader[:2, :].read()):
            assert isinstance(result.val, cupy.ndarray)
        np.testing.assert_array_equal(cupy.asnumpy(reader.read().val), reader.read(xp=np).val)


def test_writers_copy_gpu_values_to_numpy(cupy, tmp_path, monkeypatch):
    # Writing a file is an explicit request to leave the GPU: every writer
    # accepts CuPy-backed data (and GPU subsets) and writes the same values.
    # Bgen, which needs qctool, has its own test below.
    from pysnptools.distreader import DistData, DistHdf5, DistMemMap, DistNpz
    from pysnptools.kernelreader import KernelData, KernelHdf5, KernelNpz
    from pysnptools.pstreader import PstData, PstHdf5, PstMemMap, PstNpz
    from pysnptools.snpreader import (Dat, Dense, DistributedBed, Ped, Pheno, SnpHdf5,
                                      SnpMemMap, SnpNpz)

    monkeypatch.delenv("ARRAY_MODULE", raising=False)
    iid = [["f", str(row)] for row in range(4)]
    sid = ["s0", "s1", "s2"]
    pos = [[1, 0, 1], [1, 0, 2], [1, 0, 3]]
    snp_values = np.array([[0, 1, 2], [2, 1, 0], [1, 1, 1], [0, 2, np.nan]], dtype=np.float32)
    snp_gpu = SnpData(iid=iid, sid=sid, pos=pos, val=cupy.asarray(snp_values), xp=cupy)
    kernel_values = np.array([[2.0, 1.0], [1.0, 2.0]])
    kernel_gpu = KernelData(iid=iid[:2], val=cupy.asarray(kernel_values), xp=cupy)
    pst_gpu = PstData(row=["r0", "r1", "r2", "r3"], col=sid, val=cupy.asarray(snp_values), xp=cupy)
    dist_values = np.full((4, 3, 3), 1 / 3)
    dist_gpu = DistData(iid=iid, sid=sid, pos=pos, val=cupy.asarray(dist_values), xp=cupy)

    def path(name):
        return str(tmp_path / name)

    snp_writers = {
        "SnpNpz": lambda data, label: SnpNpz.write(path(label + ".snp.npz"), data),
        "Bed": lambda data, label: Bed.write(path(label + ".bed"), data, count_A1=True),
        "Dat": lambda data, label: Dat.write(path(label + ".dat"), data),
        "Dense": lambda data, label: Dense.write(path(label + ".dense.txt"), data),
        "Ped": lambda data, label: Ped.write(path(label + ".ped"), data),
        "SnpHdf5": lambda data, label: SnpHdf5.write(path(label + ".snp.hdf5"), data),
        "SnpMemMap": lambda data, label: SnpMemMap.write(path(label + ".snp.memmap"), data),
        "DistributedBed": lambda data, label: DistributedBed.write(path(label + "-distributed"), data),
    }
    reader_writers = {"SnpMemMap", "DistributedBed"}  # the others take in-memory SnpData only
    for name, write in snp_writers.items():
        sources = [("data", snp_gpu)] + ([("subset", snp_gpu[:, :])] if name in reader_writers else [])
        for label, source in sources:
            written = write(source, label)
            values = written.read(xp=np).val
            if name == "Ped":  # PED may count either allele, so a column can read back as 2 - value
                for column in range(snp_values.shape[1]):
                    expected = snp_values[:, column]
                    assert (np.array_equal(values[:, column], expected, equal_nan=True)
                            or np.array_equal(values[:, column], 2 - expected, equal_nan=True)), column
            else:
                np.testing.assert_array_equal(values, snp_values, err_msg=f"{name} {label}")
    pheno = Pheno.write(path("a.pheno.txt"), snp_gpu[:, :1].read(xp=cupy))
    np.testing.assert_array_equal(pheno.read(xp=np).val, snp_values[:, :1])

    for written in (KernelNpz.write(path("a.kernel.npz"), kernel_gpu),
                    KernelHdf5.write(path("a.kernel.hdf5"), kernel_gpu)):
        np.testing.assert_array_equal(written.read(xp=np).val, kernel_values)
    for written in (PstNpz.write(path("a.pst.npz"), pst_gpu), PstHdf5.write(path("a.pst.hdf5"), pst_gpu),
                    PstMemMap.write(path("a.pst.memmap"), pst_gpu)):
        np.testing.assert_array_equal(written.read(xp=np).val, snp_values)
    for written in (DistNpz.write(path("a.dist.npz"), dist_gpu), DistHdf5.write(path("a.dist.hdf5"), dist_gpu),
                    DistMemMap.write(path("a.dist.memmap"), dist_gpu),
                    DistMemMap.write(path("b.dist.memmap"), dist_gpu[:, :])):
        np.testing.assert_allclose(written.read(xp=np).val, dist_values, rtol=1e-15)


def test_bgen_writer_copies_gpu_values_to_numpy(cupy, tmp_path, monkeypatch):
    # Bgen.write runs the external qctool; report a skip rather than passing silently.
    import os
    from pysnptools.distreader import Bgen, DistData

    if not os.environ.get("QCTOOLPATH"):
        pytest.skip("Bgen.write requires external qctool; set QCTOOLPATH (see the Testing section in README.md)")
    monkeypatch.delenv("ARRAY_MODULE", raising=False)
    iid = [["f", str(row)] for row in range(4)]
    sid = ["s0", "s1", "s2"]
    dist_values = np.full((4, 3, 3), 1 / 3)
    dist_gpu = DistData(iid=iid, sid=sid, pos=[[1, 0, 1], [1, 0, 2], [1, 0, 3]],
                        val=cupy.asarray(dist_values), xp=cupy)
    for label, source in (("data", dist_gpu), ("subset", dist_gpu[:, :])):
        written = Bgen.write(str(tmp_path / f"{label}.bgen"), source, bits=32)
        np.testing.assert_allclose(written.read(xp=np).val, dist_values, atol=1e-6)
        del written


def test_writers_keep_layout_and_reversed_arguments(cupy, tmp_path, monkeypatch):
    # The NumPy copy keeps F/C layout, and the deprecated (data, filename)
    # order is normalized before the copy.
    from pysnptools.snpreader import Dat, Dense, Ped
    from pysnptools.util import _device, asnumpy

    monkeypatch.delenv("ARRAY_MODULE", raising=False)
    values = np.array([[0, 1, 2], [2, 1, 0], [1, 1, 1], [0, 2, 2]], dtype=np.float64)
    for order in ("C", "F"):
        gpu_values = cupy.array(values, order=order)
        host = asnumpy(gpu_values)
        assert host.flags["C_CONTIGUOUS" if order == "C" else "F_CONTIGUOUS"], order
        data = SnpData(iid=[["f", str(row)] for row in range(4)], sid=["s0", "s1", "s2"],
                       val=gpu_values, xp=cupy)
        copied = _device.numpy_for_write(data)
        assert isinstance(copied.val, np.ndarray)
        assert copied.val.flags["C_CONTIGUOUS" if order == "C" else "F_CONTIGUOUS"], order
    for name, reader_class in (("bed", Bed), ("dat", Dat), ("ped", Ped), ("dense", Dense)):
        filename = str(tmp_path / f"reversed.{name}")
        options = {"count_A1": True} if reader_class is Bed else {}
        with pytest.warns(DeprecationWarning, match="filename before data"):
            written = reader_class.write(data, filename, **options)
        result = written.read(xp=np).val
        for column in range(values.shape[1]):  # PED may count either allele
            assert (np.array_equal(result[:, column], values[:, column])
                    or (name == "ped" and np.array_equal(result[:, column], 2 - values[:, column]))), (name, column)
