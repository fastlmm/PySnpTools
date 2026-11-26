PySnpTools
====================

PySnpTools is a library for reading and manipulating genetic data.

Main Features:

* **New: [GPU standardization and kernels](https://github.com/fastlmm/PySnpTools/blob/main/doc/source/gpu.rst).** PySnpTools reads genotypes on the CPU, then transfers them to [CuPy](https://cupy.dev/), the NumPy-compatible GPU array library, for GPU standardization and kernel construction. A main use is FaST-LMM's single-kernel `single_snp`, which ran over 10 times as fast with float32 on an RTX 3060 as with the default CPU float64, in a synthetic benchmark with 14,925 individuals and 100,000 SNPs. See also the [GPU workflow notebook](https://github.com/fastlmm/PySnpTools/blob/main/doc/ipynb/gpu_workflow.ipynb).

* [SnpReader](http://fastlmm.github.io/PySnpTools): Efficiently read genetic PLINK formats including \*.bed/bim/fam files.
          Also, efficiently read parts of files, read kernel data, and standardize data.
          Features include multi-threaded BED reading, cluster-ready BED data, on-the-fly SNP generation,
          and larger in-memory data.

* [DistReader](https://fastlmm.github.io/PySnpTools/#module-pysnptools.distreader): Efficiently work with
         unphased BGEN format and other diploid, biallelic distribution data.
          Also, efficiently read parts of files. See [Distribution IPython Notebook](https://nbviewer.jupyter.org/github/fastlmm/PySnpTools/blob/main/doc/ipynb/Dist.ipynb).

* [util](https://fastlmm.github.io/PySnpTools/#module-pysnptools.util): In one line, intersect and re-order IIDs from snpreader and other sources.
          Also, efficiently extract a submatrix from an ndarray.

* [IntRangeSet](https://fastlmm.github.io/PySnpTools/#util-intrangeset): Efficiently manipulate ranges of integers - for example, genetic position - with set operators including union, intersection, and set difference.

* [mapreduce1](https://fastlmm.github.io/PySnpTools/#module-pysnptools.util.mapreduce1): Run loops locally, on multiple processors, or on any cluster.

* [filecache](https://fastlmm.github.io/PySnpTools/#module-pysnptools.util.filecache):  Read and write files locally or from/to any remote storage.

Install
-------

    pip install pysnptools

Requires Python 3.11.5 or later. Core Python 3.15 support excludes HDF5, BGEN,
and CuPy because their required wheels are not currently available together.

Install the `cupy` extra to use an NVIDIA GPU.

    pip install "pysnptools[cupy]==1.1.0b2"

GPU support requires NVIDIA CUDA on x86-64 Windows or Linux, including WSL2.
Python 3.15 is currently CPU-only. See the GPU guide for driver and CUDA setup.

HDF5 support is optional. To read or write HDF5 files, install:

    pip install pysnptools[hdf5]

Existing HDF5 users upgrading from an earlier release should install this extra.

*If you need support for BGEN files, instead do:*

    pip install pysnptools[bgen]

Documentation
-------------

* [Main Documentation](http://fastlmm.github.io/PySnpTools/) with examples. It includes links to tutorial slides, notebooks, and video.
* [Project Home and Full Annotated Bibliography](https://fastlmm.github.io/)

Code
----

* [PyPi](https://pypi.org/project/pysnptools/)
* [GitHub](https://github.com/fastlmm/PySnpTools)
* [Change Log](https://github.com/fastlmm/PySnpTools/blob/main/CHANGELOG.md)

Contacts
--------

* Email the developers at <fastlmm-dev@python.org>.
* [Join](mailto:fastlmm-user-join@python.org?subject=Subscribe) the user discussion and announcement list (or use [web sign up](https://mail.python.org/mailman3/lists/fastlmm-user.python.org)).
* [Open an issue](https://github.com/fastlmm/PySnpTools/issues) on GitHub.

## Testing

PySnpTools has two test runners. Run `pytest` from the repository root for
pytest tests and doctests in Python modules, reStructuredText, and Markdown.
Run `python -m pysnptools.test` for the legacy unittest suite, which assembles
reader, standardizer, utility, optional-format, and doctest checks. These runners
cover overlapping but different sets of tests.

Two doctest examples in `pysnptools.distreader.bgen` require external QCTool v2:
`Bgen.genwrite` and `Bgen.write`. Install QCTool v2 from the [official QCTool
site](https://www.well.ox.ac.uk/~gav/qctool/), then set `QCTOOLPATH` to the
command or executable path used to launch it. For example, on Linux or macOS:

```sh
export QCTOOLPATH=/absolute/path/to/qctool
```

On Windows, install QCTool in WSL and set `QCTOOLPATH` to the WSL launch command
and Linux executable path, as described in the `Bgen.write` API documentation.
The BGEN doctest fails with a message naming qctool when this configuration is
missing; it does not silently skip.
