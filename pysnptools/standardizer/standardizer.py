import warnings
import numpy as np
import logging
import pysnptools.util as pstutil
from bed_reader import get_num_threads, standardize_f64, standardize_f32


def _check_observed(imissX):
    # Match the Rust standardizer, which raises for a SNP with no observed values. (The
    # Python CPU path would otherwise silently set that SNP to 0.)
    observed_count = (~imissX).sum(axis=0)
    if bool((observed_count == 0).any()):
        raise ValueError("No individual observed for the SNP.")


def _is_constant(snps, axis=0):
    # A SNP whose observed values are all equal is constant. Test that directly rather
    # than relying on a zero variance: a fractional mean rounds (0.1 averages to
    # 0.09999999999999999), leaving a tiny nonzero variance that would standardize
    # every value to +/-1. Requires at least one observed value per SNP.
    xp = pstutil.get_array_module(snps)
    return xp.nanmin(snps, axis=axis) == xp.nanmax(snps, axis=axis)


class Standardizer(object):
    """
    A Standardizer is a class such as :class:`.Unit` and :class:`.Beta` to be used by the :meth:`.SnpData.standardize` and :meth:`.SnpReader.read_kernel` method to standardize SNP data.

    :Example:

    Read and standardize SNP data.

    >>> from pysnptools.standardizer import Unit
    >>> from pysnptools.snpreader import Bed
    >>> from pysnptools.util import example_file # Download and return local file name
    >>> bedfile = example_file("tests/datasets/all_chr.maf0.001.N300.*","*.bed")
    >>> snpdata1 = Bed(bedfile,count_A1=False).read().standardize(Unit())
    >>> print('{0:.6f}'.format(snpdata1.val[0,0]))
    0.229416

    Create a kernel from SNP data on disk.

    >>> bedfile2 = example_file("pysnptools/examples/toydata.5chrom.*","*.bed")
    >>> kerneldata = Bed(bedfile2,count_A1=False).read_kernel(Unit())
    >>> print('{0:.6f}'.format(kerneldata.val[0,0]))
    9923.069928


    Can also return a constant SNP standardizer that can be applied to other :class:`.SnpData`.

    >>> snp_whole = Bed(bedfile,count_A1=False)
    >>> train_idx, test_idx = range(10,snp_whole.iid_count), range(0,10) #test on the first 10, train on the rest
    >>> snp_train, trained_standardizer = Unit().standardize(snp_whole[train_idx,:].read(),return_trained=True)
    >>> print('{0:.6f}'.format(snp_train.val[0,0]))
    0.233550
    >>> print(trained_standardizer.stats[0,:]) #The mean and stddev of the 1st SNP on the training data # '...' for a possible space character
    [...1.94827586  0.22146953]
    >>> snp_test = snp_whole[test_idx,:].read().standardize(trained_standardizer)
    >>> float(snp_test.val[0,0])
    0.23354968324845735

    Standardize any Numpy array.

    >>> val = Bed(bedfile,count_A1=False).read().val
    >>> print(val[0,0])
    2.0
    >>> val = Unit().standardize(val)
    >>> print('{0:.6f}'.format(val[0,0]))
    0.229416

    Details of Methods & Properties:
    """

    def __init__(self):
        super(Standardizer, self).__init__()

    def standardize(
        self,
        snps,
        block_size=None,
        return_trained=False,
        force_python_only=False,
        num_threads=None,
    ):
        """
        Applies standardization in place on the existing NumPy/CuPy device.
        Returns the same input, or (input, trained standardizer) when requested.
        Raw arrays remain accepted but are deprecated; return_trained=True for
        Unit/Beta requires SnpData. Values must be contiguous float32/float64.
        See :doc:`/gpu` for complete examples and Unit/Beta exceptions.

        :param snps: SNP values to standardize, on their existing device
        :type snps: :class:`.SnpData` (or a NumPy/CuPy array)

        :param block_size: Deprecated and ignored; direct Unit/Beta calls warn.
        :type block_size: None

        :param return_trained: If true, returns a second value containing a constant standardizer trained on this data.
        :type return_trained: bool

        :param force_python_only: Request Python CPU standardization where supported.
            Ignored by CuPy, which stays on its current device. This is not a
            backend switch; see :doc:`/gpu` for CPU Beta fallback limitations.
        :type force_python_only: bool

        :param num_threads: optional -- The CPU standardization thread count (ignored for CuPy). Defaults to all available
            processors. Can also be set with these environment variables (listed in priority order):
            'PST_NUM_THREADS', 'NUM_THREADS', 'MKL_NUM_THREADS'.
        :type num_threads: None or int

        :rtype: SnpData or array, or a tuple with a trained Standardizer

        """
        if block_size is not None:
            warnings.warn(
                "block_size is deprecated (and not needed, since standardization is in-place",
                DeprecationWarning,
            )
        raise NotImplementedError(
            "subclass {0} needs to implement method '.standardize'".format(
                self.__class__.__name__
            )
        )

    @staticmethod
    # changes snps in place
    def _standardize_unit_and_beta(
        snps,
        is_beta,
        a,
        b,
        apply_in_place,
        use_stats,
        stats,
        num_threads,
        force_python_only=False,
    ):
        """
        Dispatch Unit/Beta training and application from the actual array module.
        """
        xp = pstutil.get_array_module(snps)

        assert (
            snps.flags["C_CONTIGUOUS"] or snps.flags["F_CONTIGUOUS"]
        ), "Expect snps to be order 'C' or order 'F'"

        if xp is not np:
            from pysnptools.util._gpu_warnings import warn_gpu_float64
            warn_gpu_float64(snps.dtype, xp)
            return Standardizer._standardize_cupy(
                snps, is_beta, a, b, apply_in_place, use_stats, stats, xp)

        # Make sure stats is the same type as snps. Because we might be creating a new array, we return it
        if stats is None:
            stats = xp.empty(
                [snps.shape[1], 2],
                dtype=snps.dtype,
                order="F" if snps.flags["F_CONTIGUOUS"] else "C",
            )
        elif not (
            stats.dtype == snps.dtype  # stats must have the same dtype as snps
            and (stats.flags["OWNDATA"])  # stats must own its data
            and (
                (snps.flags["C_CONTIGUOUS"] and stats.flags["C_CONTIGUOUS"])
                or (snps.flags["F_CONTIGUOUS"] and stats.flags["F_CONTIGUOUS"])
            )  # Both layouts also require the matching dtype and ownership above.
        ):
            stats = xp.array(
                stats,
                dtype=snps.dtype,
                order="F" if snps.flags["F_CONTIGUOUS"] else "C",
            )
        # Raise explicitly, not with assert, so the check survives python -O.
        if stats.shape != (snps.shape[1], 2):
            raise AssertionError("stats must have size [sid_count,2]")

        if not force_python_only and xp is np:
            num_threads = get_num_threads(num_threads)

            if snps.dtype == np.float64:
                if (snps.flags["F_CONTIGUOUS"] or snps.flags["C_CONTIGUOUS"]) and (
                    snps.flags["OWNDATA"] or snps.base.nbytes == snps.nbytes
                ):  # !!create a method called is_single_segment
                    standardize_f64(
                        snps,
                        is_beta,
                        a,
                        b,
                        apply_in_place,
                        use_stats,
                        stats,
                        num_threads,
                    )
                    return stats
                else:
                    logging.info(
                        "Array is not contiguous, so will standardize with python only instead of C++"
                    )
            elif snps.dtype == np.float32:
                if (snps.flags["F_CONTIGUOUS"] or snps.flags["C_CONTIGUOUS"]) and (
                    snps.flags["OWNDATA"] or snps.base.nbytes == snps.nbytes
                ):
                    standardize_f32(
                        snps,
                        is_beta,
                        a,
                        b,
                        apply_in_place,
                        use_stats,
                        stats,
                        num_threads,
                    )
                    return stats
                else:
                    logging.info(
                        "Array is not contiguous, so will standardize with python only instead of C++"
                    )
            else:
                logging.info(
                    "Array type is not float64 or float32, so will standardize with python only instead of C++"
                )


        if is_beta:
            Standardizer._standardize_beta_python(
                snps, a, b, apply_in_place, use_stats=use_stats, stats=stats
            )
            return stats
        else:
            Standardizer._standardize_unit_python(
                snps, apply_in_place, use_stats=use_stats, stats=stats
            )
            return stats

    @staticmethod
    def _standardize_cupy(snps, is_beta, beta_a, beta_b, apply_in_place,
                          use_stats, stats, xp):
        """Use reductions and an in-place update, with only per-SNP workspaces.

        Centered float64 reductions also support dosages and large offsets without
        the cancellation in E[x*x] - E[x]**2. Transpose F-order inputs so the
        elementwise update traverses contiguous memory.
        """
        from pysnptools.standardizer._cupy import kernels

        if snps.dtype not in (np.float32, np.float64):
            raise ValueError("Standardizing requires float32 or float64")
        count_kernel, sum_kernel, variance_kernel, apply_kernel = kernels()
        values = snps.T if snps.flags["F_CONTIGUOUS"] else snps
        axis = 1 if snps.flags["F_CONTIGUOUS"] else 0
        expand = (lambda vector: vector[:, None]) if axis == 1 else (lambda vector: vector[None, :])
        if use_stats:
            stats = xp.asarray(stats, dtype=snps.dtype)
        else:
            count = count_kernel(values, axis=axis)
            if bool((count == 0).any()):
                raise ValueError("No individual observed for the SNP.")
            mean = sum_kernel(values, axis=axis) / count
            variance = variance_kernel(values, expand(mean), axis=axis) / count
            deviation = xp.sqrt(variance).astype(snps.dtype)
            deviation[(deviation == 0) | _is_constant(values, axis)] = xp.inf
            stats = xp.empty((snps.shape[1], 2), dtype=snps.dtype)
            stats[:, 0] = mean
            stats[:, 1] = deviation
            if bool(xp.isnan(stats[:, 0]).any()) or (is_beta and bool(
                    ((stats[:, 0] < 0) | (stats[:, 0] > 2)).any())):
                raise ValueError("Illegal SNP mean.")
        if stats.shape != (snps.shape[1], 2):
            raise ValueError("stats must have size [sid_count,2]")
        if apply_in_place:
            if is_beta:
                import scipy.stats as st

                # Match statrs' Beta parameter validation used by native Rust.
                if (np.isnan(beta_a) or np.isnan(beta_b) or beta_a <= 0 or beta_b <= 0
                        or (np.isinf(beta_a) and np.isinf(beta_b))):
                    raise ValueError(f"Cannot create Beta distribution with a={beta_a}, b={beta_b}")
                mean = pstutil.asnumpy(stats[:, 0])
                maf = mean / 2.0
                maf = np.minimum(maf, 1.0 - maf)
                factor = xp.asarray(st.beta.pdf(maf, beta_a, beta_b), dtype=snps.dtype)
                # CPU Rust zeros constant columns even for singular Beta PDFs
                # at mean 0 or 2, and when applying trained statistics.
                factor[xp.isinf(stats[:, 1])] = 0
            else:
                factor = 1 / stats[:, 1]
            apply_kernel(expand(stats[:, 0]), expand(factor), values)
        return stats

    @staticmethod
    def _standardize_unit_python(snps, apply_in_place, use_stats, stats):
        """
        Standardize snps to zero-mean and unit variance.

        Will work with both numpy and cupy ndarray.
        """
        assert snps.dtype in [
            np.float64,
            np.float32,
        ], "snps must be a float in order to standardize in place."
        xp = pstutil.get_array_module(snps)

        imissX = xp.isnan(snps)

        if use_stats:
            snp_mean = stats[:, 0]
            snp_std = stats[:, 1]
        else:
            _check_observed(imissX)
            snp_std = xp.nanstd(snps, axis=0)
            snp_mean = xp.nanmean(snps, axis=0)
            # avoid div by 0 when standardizing
            # Don't need this warning because SNCs are still meaning full in QQ plots because they should be thought of as SNPs without enough data.
            # logging.warn("A least one snps has only one value, that is, its standard deviation is zero")
            snp_std[(snp_std == 0.0) | _is_constant(snps)] = (
                xp.inf
            )  # We make the stdev infinity so that applying as a trained_standardizer will turn any input to 0. Thus if a variable has no variation in the training data, then it will be set to 0 in test data, too.
            stats[:, 0] = snp_mean
            stats[:, 1] = snp_std

        if apply_in_place:
            snps -= snp_mean
            snps /= snp_std
            snps[imissX] = 0

    @property
    def is_constant(self):
        """
        Tell if a standardizer's statistics are constant. For example, :class:`.UnitTrained` contains a pre-specified mean and stddev for each SNP, so :attr:`UnitTrained.is_constant` is True.
        :class:`.Unit`, on the other hand, measures the mean and stddev for each SNP, so :attr:`Unit.is_constant` is False.

        :rtype: bool
        """
        return False

    @staticmethod
    def _standardize_beta_python(snps, betaA, betaB, apply_in_place, use_stats, stats):
        """
        standardize snps with Beta prior
        """
        assert snps.dtype in [
            np.float64,
            np.float32,
        ], "snps must be a float in order to standardize in place."

        xp = pstutil.get_array_module(snps)
        imissX = np.isnan(snps)
        snp_sum = np.nansum(snps, axis=0)
        n_obs_sum = (~imissX).sum(0)

        if use_stats:
            snp_mean = stats[:, 0]
            snp_std = stats[:, 1]
        else:
            _check_observed(imissX)
            snp_mean = (snp_sum * 1.0) / n_obs_sum
            # Match the native and CuPy paths, which reject a mean outside [0, 2]
            # before changing any values.
            if bool(((snp_mean < 0) | (snp_mean > 2)).any()):
                raise ValueError("Illegal SNP mean.")
            snp_std = np.sqrt(np.nansum((snps - snp_mean) ** 2, axis=0) / n_obs_sum)
            # Don't need this warning because SNCs are still meaning full in QQ plots because they should be thought of as SNPs without enough data.
            # logging.warn("A least one snps has only one value, that is, its standard deviation is zero")
            snp_std[(snp_std == 0) | _is_constant(snps)] = np.inf
            stats[:, 0] = snp_mean
            stats[:, 1] = snp_std

        if apply_in_place:
            maf = snp_mean / 2.0
            maf[maf > 0.5] = 1.0 - maf[maf > 0.5]

            # avoid div by 0 when standardizing
            import scipy.stats as st

            # scipy needs numpy; maf has one value per SNP, so the round trip is cheap.
            maf_beta = xp.asarray(st.beta.pdf(pstutil.asnumpy(maf), betaA, betaB))
            # Constant training SNPs need zero weight before multiplication:
            # endpoint Beta densities can be infinite, so even centered zeros
            # would otherwise produce NaNs (and varying test values infinities).
            maf_beta[snp_std == np.inf] = 0.0
            # print("BetaPdf[{0},{1},{2}]={3}".format(maf,betaA,betaB,maf_beta))
            snps -= snp_mean
            snps *= maf_beta
            snps[imissX] = 0.0
            if (
                use_stats
            ):  # If we're applying to test data, set any variables with to 0 if they have no variation in the training data.
                snps[:, snp_std == np.inf] = 0.0

    def _merge_trained(self, trained_list):
        raise Exception("Not defined")


class _CannotBeTrained(Standardizer):
    def __init__(self, name):
        super(_CannotBeTrained, self).__init__()
        self.name = name

    def __repr__(self):
        return "{0}({1})".format(self.__class__.__name__, self.name)

    def standardize(
        self,
        snps,
        block_size=None,
        return_trained=False,
        force_python_only=False,
        num_threads=None,
    ):
        if block_size is not None:
            warnings.warn(
                "block_size is deprecated (and not needed, since standardization is in-place",
                DeprecationWarning,
            )
        raise Exception("Standardizer '{0}' cannot be trained", self)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    import doctest

    doctest.testmod(optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE)
