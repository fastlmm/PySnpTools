PySnpTools
====================

PySnpTools is a library for reading and manipulating genetic data.

Main Features:

* [SnpReader](http://fastlmm.github.io/PySnpTools): Efficiently read genetic PLINK formats including \*.bed/bim/fam files.
          Also, efficiently read parts of files, read kernel data, and standardize data. 
          New features include on-the-fly SNP generation, larger in-memory data, and
          cluster-ready BED data.

* [DistReader](https://fastlmm.github.io/PySnpTools/#module-pysnptools.distreader): Efficiently work with
         unphased BGEN format and other diploid, biallelic distribution data.
          Also, efficiently read parts of files. See [Distribution IPython Notebook](https://nbviewer.jupyter.org/github/fastlmm/PySnpTools/blob/master/doc/ipynb/Dist.ipynb).

* [util](https://fastlmm.github.io/PySnpTools/#module-pysnptools.util): In one line, intersect and re-order IIDs from snpreader and other sources.
          Also, efficiently extract a submatrix from an ndarray. 

* [IntRangeSet](https://fastlmm.github.io/PySnpTools/#util-intrangeset): Efficiently manipulate ranges of integers - for example, genetic position - with set operators including union, intersection, and set difference. 

* [mapreduce1](https://fastlmm.github.io/PySnpTools/#module-pysnptools.util.mapreduce1): Run loops locally, on multiple processors, or on any cluster. 

* [filecache](https://fastlmm.github.io/PySnpTools/#module-pysnptools.util.filecache):  Read and write files locally or from/to any remote storage.

Documentation
=================================

* [API Documentation](http://fastlmm.github.io/PySnpTools/) with examples. It includes links to tutorial slides, notebooks, and video.
* [Project Home and Full Annotated Bibliography](https://fastlmm.github.io/)

Code
=================================
* [PyPi](https://pypi.org/project/pysnptools/)
* [GitHub](https://github.com/fastlmm/PySnpTools)

Contacts
=================================

* Email the developers at fastlmm-dev@python.org.
* [Join](mailto:fastlmm-user-join@python.org?subject=Subscribe) the user discussion and announcement list (or use [web sign up](https://mail.python.org/mailman3/lists/fastlmm-user.python.org)).
* [Open an issue](https://github.com/fastlmm/PySnpTools/issues) on GitHub.


Quick install:
====================

`pip install pysnptools`
