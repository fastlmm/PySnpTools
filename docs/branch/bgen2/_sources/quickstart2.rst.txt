.. _quick_start2:

***********
Quick Start
***********

We first download the ``example.bgen``, ``haplotypes.bgen``, and ``complex.bgen`` files:

.. code-block:: python

    >>> from bgen_reader import example_filepath
    >>>
    >>> example_bgen_path = example_filepath("example.bgen")
    >>> haplotypes_bgen_path = example_filepath("haplotypes.bgen")
    >>> complex_bgen_path = example_filepath("complex.bgen")

This software reads the files the same way, but they present
different levels of genotype complexity.
The variant loci of the genotype stored in ``example.bgen`` have the ploidy two and are all unphased.
The ``haplotypes.bgen`` file stores a phased, diploid genotype.
Both the genotypes stored in ``example.bgen`` and ``haplotypes.bgen`` are biallelic, that is with two alleles
at each loci.
The ``complex.bgen`` file, on the other hand, has varying ploidy levels and alleles numbers
across the variant loci.


.. _unphased_genotype2:

Unphased genotype
=================

With unphased genotype we don't have the information about which chromosome holds a
given allele.
Let's read the ``example.bgen`` file and print out some information.

.. doctest::

   >>> from bgen_reader import open_bgen
   >>>
   >>> bgen = open_bgen(example_bgen_path, verbose=False)
   >>>
   >>> # Samples
   >>> print(bgen.samples[:5]) #first 5
   ['sample_001' 'sample_002' 'sample_003' 'sample_004' 'sample_005']
   >>>
   >>> # Variants ids.
   >>> print(bgen.ids[:5]) #first 5
   ['SNPID_2' 'SNPID_3' 'SNPID_4' 'SNPID_5' 'SNPID_6']
   >>>
   >>> # There are 199 variants in total.
   >>> print(bgen.nvariants)
   199
   >>> # Read the probabilities of the first variant
   >>> probs = bgen.read(0)
   >>> print(probs)
   [[[       nan        nan        nan]]
   <BLANKLINE>
    [[0.02780236 0.00863674 0.9635609 ]]
   <BLANKLINE>
    [[0.01736504 0.04968414 0.93295083]]
   <BLANKLINE>
    ...
   <BLANKLINE>
    [[0.01419069 0.02810669 0.95770262]]
   <BLANKLINE>
    [[0.91949463 0.05206298 0.02844239]]
   <BLANKLINE>
    [[0.00244141 0.98410029 0.0134583 ]]]
   >>> # The above matrix is of size samples-by-1-variants-by-(combination-of-alleles).
   >>> print(probs.shape)
   (500, 1, 3)
   >>>
   >>> # Read the probabilities of all the samples and all variants
   >>> probs = bgen.read()
   >>> # The matrix is of size samples-by-variants-(combination-of-alleles).
   >>> print(probs.shape)
   (500, 199, 3)
   >>>
   >>> # Read the probabilities of the first 5 samples and
   >>> #variants 15 (inclusive) to 25 (exclusive)
   >>> probs = bgen.read((slice(5),slice(15,25)))
   >>> print(probs.shape)
   (5, 10, 3)
   >>>
   >>> del bgen  #close the file and delete the open_bgen object

The 3rd dimension of the probability matrix of a given variant depends on the
number of alleles, the ploidy, and whether the locus is phased or unphased.
See :attr:`bgen_reader.open_bgen.ncombinations`. For details, see the
|bgen specification|.

.. _phased_genotype2:

Phased genotype
===============

.. doctest::

   >>> bgen = open_bgen(haplotypes_bgen_path, verbose=False)
   >>>
   >>> # Samples
   >>> print(bgen.samples)
   ['sample_0' 'sample_1' 'sample_2' 'sample_3']
   >>>
   >>> # Variants ids.
   >>> print(bgen.ids)
   ['SNP1' 'SNP2' 'SNP3' 'SNP4']
   >>>
   >>> # Read the probabilities and ploidies for the second individual and first variant
   >>> probs,ploidy = bgen.read((1,0),return_ploidies=True)
   >>> print(probs)
   [[[0. 1. 1. 0.]]]
   >>> # How many haplotypes?
   >>> print(ploidy)
   [[2]]
   >>> # Is the first variant phased?
   >>> print(bgen.phased[0])
   True
   >>> # And how many alleles for the first variant?
   >>> print(bgen.nalleles[0])
   2
   >>> # Therefore, the first haplotype has probability 100%
   >>> # of having the allele
   >>> alleles = bgen.allele_ids[0].split(",")
   >>> print(alleles[1])
   G
   >>> # And the second haplotype has probability 100% of having
   >>> # the first allele
   >>> print(alleles[0])
   A
   >>>
   >>> del bgen  #close the file and delete the open_bgen object


Please refer to |bgen specification| for a detailed description.

.. _complex_file2:

Complex file
============

The bgen file format allows the storage of very heterogeneous genetic data.
In the ``complex.bgen`` file we have variants with different ploidy and numbers of
alleles, as well as phased\ *ness*.

.. doctest::

   >>> bgen = open_bgen(complex_bgen_path, verbose=False)
   >>>
   >>> # Note how the number of alleles very widely across loci.
   >>> print(bgen.allele_ids)
   ['A,G' 'A,G' 'A,G' 'A,G,T' 'A,G' 'A,G,GT,GTT' 'A,G,GT,GTT,GTTT,GTTTT'
    'A,G,GT,GTT,GTTT,GTTTT,GTTTTT' 'A,G,GT,GTT,GTTT,GTTTT,GTTTTT,GTTTTTT'
    'A,G']
   >>> print(bgen.samples)
   ['sample_0' 'sample_1' 'sample_2' 'sample_3']
   >>> # Print the probabilities for the second individual and first variant
   >>> probs = bgen.read((1,0),max_combinations=bgen.ncombinations[0])
   >>> print(probs)
   [[[1. 0. 0.]]]
   >>> # The 4th individual and 9th variant has ploidy ...
   >>> probs, ploidy = bgen.read((3,8),max_combinations=bgen.ncombinations[8],return_ploidies=True)
   >>> print(ploidy)
   [[2]]
   >>> # and number of alleles equal to ...
   >>> print(bgen.nalleles[8])
   8
   >>> # Its probability distribution is given by the array
   >>> print(probs)
   [[[0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0.
      0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]]]
   >>> # Since the 9th variant is unphased,
   >>> print(bgen.phased[8])
   False
   >>> # we can pick an alternative allele and compute the dosage
   >>> # from allele expectation.
   >>> # If we select the third allele as being the alternative one, we have
   >>> e = bgen.allele_expectation(8)
   >>> dosage = e[:,0,2]
   >>> print(dosage)
   [0. 0. 0. 1.]

Please, refer to :meth:`bgen_reader.open_bgen.allele_expectation` for further details.

.. |bgen specification| raw:: html

   <a href="https://github.com/limix/bgen" target="_blank">bgen specification⧉</a>
