***********
Quick Start
***********

We first download the ``example.bgen``, ``haplotypes.bgen``, and ``complex.bgen`` files:

.. code-block:: python

    >>> from bgen_reader import example_filepath
    >>>
    >>> filepath = {}
    >>> filepath["example.bgen"] = example_filepath("example.bgen")
    >>> filepath["haplotypes.bgen"] = example_filepath("haplotypes.bgen")
    >>> filepath["complex.bgen"] = example_filepath("complex.bgen")

The above-mentioned files are read by this software in the same way but they present
different levels of genotype complexity.
The variant loci of the genotype stored in ``example.bgen`` have the same ploidy equal
to two and are all unphased.
The ``haplotypes.bgen`` file stores a phased, diploid genotype.
Both genotypes stored in ``example.bgen`` and ``haplotypes.bgen`` have the same number
of alleles for each loci.
The ``complex.bgen`` file, on the other hand, have different ploidy levels and number of
alleles across different variant loci.


Unphased genotype
=================

With unphased genotype we don't have the information about which chromosome holds a
given allele.
Let's read the ``example.bgen`` file and print out some information.

.. doctest::

   >>> from bgen_reader import read_bgen
   >>>
   >>> bgen = read_bgen(filepath["example.bgen"], verbose=False)
   >>>
   >>> # Variants metadata.
   >>> print(bgen["variants"].head())
           id    rsid chrom   pos  nalleles allele_ids  vaddr
   0  SNPID_2  RSID_2    01  2000         2        A,G   6069
   1  SNPID_3  RSID_3    01  3000         2        A,G   9635
   2  SNPID_4  RSID_4    01  4000         2        A,G  12956
   3  SNPID_5  RSID_5    01  5000         2        A,G  16034
   4  SNPID_6  RSID_6    01  6000         2        A,G  19377
   >>> # Samples read from the bgen file.
   >>> print(bgen["samples"].head())
   0    sample_001
   1    sample_002
   2    sample_003
   3    sample_004
   4    sample_005
   Name: id, dtype: object
   >>> # There are 199 variants in total.
   >>> print(len(bgen["genotype"]))
   199
   >>> # This library avoid as much as possible accessing the bgen file for performance
   >>> # and memory reasons. The `compute` function actually tells the library to
   >>> # access the file to retrieve some data.
   >>> geno = bgen["genotype"][0].compute()
   >>> print(geno.keys())
   dict_keys(['probs', 'phased', 'ploidy', 'missing'])
   >>> # Let's have a look at the probabilities regarding the first variant.
   >>> print(geno["probs"])
   [[       nan        nan        nan]
    [0.02780236 0.00863674 0.9635609 ]
    [0.01736504 0.04968414 0.93295083]
    ...
    [0.01419069 0.02810669 0.95770262]
    [0.91949463 0.05206298 0.02844239]
    [0.00244141 0.98410029 0.0134583 ]]
   >>> # The above matrix is of size samples-by-(combination-of-alleles).
   >>> print(geno["probs"].shape)
   (500, 3)

The columns dimension of the probability matrix of a given variant depend on the
number of alleles, the ploidy, and whether the locus is phased or unphased.
Please, refer to |bgen specification| for a detailed description.

Phased genotype
===============

.. doctest::

   >>> bgen = read_bgen(filepath["haplotypes.bgen"], verbose=False)
   >>>
   >>> print(bgen["variants"].head(4))
        id rsid chrom  pos  nalleles allele_ids  vaddr
   0  SNP1  RS1     1    1         2        A,G    102
   1  SNP2  RS2     1    2         2        A,G    159
   2  SNP3  RS3     1    3         2        A,G    216
   3  SNP4  RS4     1    4         2        A,G    273
   >>> print(bgen["samples"].head())
    0    sample_0
    1    sample_1
    2    sample_2
    3    sample_3
    Name: id, dtype: object
   >>> # Print the estimated probabilities for the first variant
   >>> # and second individual.
   >>> geno = bgen["genotype"][0].compute()
   >>> print(geno["probs"][1])
   [0. 1. 1. 0.]
   >>> # Is it a phased one?
   >>> print(geno["phased"])
   1
   >>> # How many haplotypes for each sample?
   >>> print(geno["ploidy"])
   [2 2 2 2]
   >>> # And how many alleles?
   >>> variant = bgen["variants"].compute()
   >>> print(variant.loc[0, "nalleles"])
   2
   >>> # Therefore, the first haplotype has probability 100%
   >>> # of having the allele
   >>> alleles = variant.loc[0, "allele_ids"].split(",")
   >>> print(alleles[1])
   G
   >>> # And the second haplotype has probability 100% of having
   >>> # the first allele
   >>> print(alleles[0])
   A

Please, refer to |bgen specification| for a detailed description.

Complex file
============

The bgen file format allows the storage of very heterogeneous genetic data.
In the ``complex.bgen`` file we have variants with different ploidy and number of
alleles, as well as phased\ *ness*.

.. doctest::

   >>> bgen = read_bgen(filepath["complex.bgen"], verbose=False)
   >>>
   >>> # Note how the number of alleles very widely across loci.
   >>> print(bgen["variants"].compute())
         id rsid chrom  pos  nalleles                            allele_ids  vaddr
   0          V1    01    1         2                                   A,G     98
   1   V2.1   V2    01    2         2                                   A,G    175
   2          V3    01    3         2                                   A,G    232
   3          M4    01    4         3                                 A,G,T    305
   ..   ...  ...   ...  ...       ...                                   ...    ...
   6          M7    01    7         6                 A,G,GT,GTT,GTTT,GTTTT    557
   7          M8    01    8         7          A,G,GT,GTT,GTTT,GTTTT,GTTTTT    663
   8          M9    01    9         8  A,G,GT,GTT,GTTT,GTTTT,GTTTTT,GTTTTTT    783
   9         M10    01   10         2                                   A,G    863
   <BLANKLINE>
   [10 rows x 7 columns]
   >>> print(bgen["samples"])
   0    sample_0
   1    sample_1
   2    sample_2
   3    sample_3
   Name: id, dtype: object
   >>> # Print the estimated probabilities for the first variant
   >>> # and second individual.
   >>> geno = bgen["genotype"][0].compute()
   >>> print(geno["probs"][1])
   [1. 0. 0.]
   >>> # The 9th variant for the 4th individual has ploidy
   >>> geno = bgen["genotype"][8].compute()
   >>> ploidy = geno["ploidy"][3]
   >>> print(ploidy)
   2
   >>> # and number of alleles equal to
   >>> nalleles = bgen["variants"].loc[8, "nalleles"].compute().values[0]
   >>> print(nalleles)
   8
   >>> # Its probability distribution is given by the array
   >>> p = geno["probs"][3]
   >>> print(p)
   [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.
    0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
   >>> # Since the 9th variant for the 4th individual is unphased,
   >>> print(geno["phased"])
   0
   >>> # we can pick an alternative allele and compute the dosage
   >>> # from allele expectation.
   >>> # If we select the third allele as being the alternative one, we have
   >>> from bgen_reader import allele_expectation, compute_dosage
   >>> e = allele_expectation(bgen, 8)
   >>> print(compute_dosage(e, 2))
   [0. 0. 0. 1.]

Please, refer to :ref:`Dosage` section for further details.

.. |bgen specification| raw:: html

   <a href="https://github.com/limix/bgen" target="_blank">bgen specification⧉</a>
