"""Independent expectations for fractional Unit/Beta arrays and trained statistics."""
import numpy as np
import pytest
from scipy.stats import beta as beta_distribution

from pysnptools.snpreader import SnpData
from pysnptools.standardizer import Unit, Beta


def reference(values, standardizer):
    expected = np.zeros_like(values)
    stats = np.empty((values.shape[1], 2), dtype=values.dtype)
    for column in range(values.shape[1]):
        valid = ~np.isnan(values[:, column])
        observed = values[valid, column].astype(np.longdouble)
        mean = observed.mean()
        variance = ((observed - mean) ** 2).mean()
        stats[column] = [mean, np.sqrt(variance) if variance > 0 else np.inf]
        if variance == 0:
            continue
        mean, deviation = stats[column]
        if isinstance(standardizer, Unit):
            factor = values.dtype.type(1) / deviation
        else:
            maf = float(mean) / 2
            factor = values.dtype.type(beta_distribution.pdf(min(maf, 1 - maf), standardizer.a, standardizer.b))
        expected[valid, column] = (values[valid, column] - mean) * factor
    return expected, stats


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("order", ["C", "F"])
@pytest.mark.parametrize("standardizer", [Unit(), Beta(2, 3), Beta(0.5, 0.5)])
def test_fractional_standardization_and_training(dtype, order, standardizer):
    cp = pytest.importorskip("cupy")
    values = np.array([
        [0.25, 0.75, 0, 2, np.nan],
        [0.5, 0.75, 0, 2, np.nan],
        [np.nan, 0.75, 0, 2, 1.125],
        [1.75, 0.75, 0, 2, np.nan],
        [1.125, 0.75, 0, 2, np.nan],
        [1.5, 0.75, 0, 2, np.nan],
    ], dtype=dtype)
    if isinstance(standardizer, Unit):
        values[:, 0] -= 4
        values[:, 1] += 3
    elif standardizer.a == 2:
        values[:, 0] = [-1.25, 0.25, np.nan, 2.5, 1.125, 1.5]
    expected, expected_stats = reference(values, standardizer)
    iid = [[str(index), str(index)] for index in range(len(values))]
    sid = [str(index) for index in range(values.shape[1])]
    # Compare the CPU and GPU public API independently to the same reference.
    for xp in [np, cp]:
        data = SnpData(iid=iid, sid=sid, val=xp.array(values, order=order), xp=xp)
        actual, trained = data.standardize(standardizer, return_trained=True)
        tolerance = 16 * np.finfo(dtype).eps
        actual_values = actual.val if xp is np else cp.asnumpy(actual.val)
        np.testing.assert_allclose(actual_values, expected, rtol=tolerance, atol=tolerance)
        np.testing.assert_allclose(trained.stats, expected_stats, rtol=tolerance, atol=tolerance)
        assert actual.val.dtype == values.dtype
        assert actual.val.flags["C_CONTIGUOUS" if order == "C" else "F_CONTIGUOUS"]
        fresh_values = values.copy()
        fresh_values[:, 1:] = 0.5  # trained constant columns must still become zero
        fresh = SnpData(iid=iid, sid=sid, val=xp.array(fresh_values, order=order), xp=xp)
        fresh.standardize(trained)
        np.testing.assert_allclose(fresh.val if xp is np else cp.asnumpy(fresh.val),
                                   expected, rtol=tolerance, atol=tolerance)


def test_beta_errors_match_cpu():
    cp = pytest.importorskip("cupy")
    iid = [["a", "a"], ["b", "b"]]
    for xp in [np, cp]:
        for standardizer, values, message in [
            (Beta(1, 25), [[-2.5], [-1.5]], "Illegal SNP mean"),
            (Beta(-1, 25), [[0.25], [1.25]], "Beta"),
            (Unit(), [[np.nan], [np.nan]], "No individual observed"),
        ]:
            data = SnpData(iid=iid, sid=["s"], val=xp.array(values, dtype=np.float64), xp=xp)
            with pytest.raises(ValueError, match=message):
                data.standardize(standardizer)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("order", ["C", "F"])
@pytest.mark.parametrize("standardizer", [Unit(), Beta(0.5, 0.5)])
def test_native_statistics_through_bed_helper(tmp_path, dtype, order, standardizer):
    from unittest.mock import patch
    from bed_reader import has_gpu_support, to_bed
    from bed_reader import bed_reader as native
    from pysnptools.snpreader import Bed

    if not has_gpu_support():
        pytest.skip("requires bed-reader built with the cubecl-cuda Cargo feature")
    cp = pytest.importorskip("cupy")
    values = np.array([[0.25, 0], [0.5, 0], [1.75, 0], [np.nan, 0]], dtype=dtype)
    if isinstance(standardizer, Unit):
        values[:, 0] -= 4
    path = tmp_path / "owned.bed"
    to_bed(path, np.zeros(values.shape, dtype=np.int8))
    reader = Bed(path, count_A1=True)
    reader._run_once()
    owner = reader._open_bed._read_gpu(dtype=dtype, order=order)
    cp.asarray(owner)[...] = cp.asarray(values)
    # BED files themselves contain only genotypes. Supply an exported, modified
    # owner to test that PySnpTools's native handoff adds no input restriction.
    with patch.object(reader._open_bed, "_read_gpu", return_value=owner), \
         patch.object(native, "standardize_gpu", wraps=native.standardize_gpu) as standardize:
        data, trained = reader._read_with_standardizing(
            standardizer, dtype=dtype, order=order, xp=cp, return_trained=True)
        assert standardize.call_count == 1
    expected, stats = reference(values, standardizer)
    tolerance = 16 * np.finfo(dtype).eps
    np.testing.assert_allclose(cp.asnumpy(data.val), expected, rtol=tolerance, atol=tolerance)
    np.testing.assert_allclose(trained.stats, stats, rtol=tolerance, atol=tolerance)
    fresh_values = values.copy()
    fresh_values[:, 1] = 0.5
    fresh = SnpData(iid=data.iid, sid=data.sid, val=cp.array(fresh_values, order=order), xp=cp)
    fresh.standardize(trained)
    np.testing.assert_allclose(cp.asnumpy(fresh.val), expected, rtol=tolerance, atol=tolerance)
