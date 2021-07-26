import os
import logging
import unittest
from pysnptools.util.filecache import Hashdown
import doctest

pysnptools_hashdown = Hashdown.load_hashdown(
    os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "pysnptools.hashdown.json"
    ),
    directory=os.environ.get("PYSNPTOOLS_CACHE_HOME", None),
)

def example_file(pattern, endswith=None):
    """
    Returns the local location of a PySnpTools example file, downloading it
    if needed.

    :param pattern: The name of the example file of interest. A 
        `file name pattern <https://docs.python.org/3.7/library/fnmatch.html>`__
        may be given. All matching files will be downloaded (if needed) and
        the name of one will be returned.
    :type file_name: string

    :param endswith: The pattern of the file name to return. By default, when
        no `endswith` is given, the name of the first matching file will be
        returned.
    :type file_name: string

    :rtype: string

    By default, the local location will be under the system temp directory
    (typically controlled with the TEMP environment variable).
    Alternatively, the directory can be set with the PYSNPTOOLS_CACHE_HOME
    environment variable.

    This function knows the MD5 hash of all PySnpTools example files and uses
    that content-based hash to decide if a file needs to be downloaded.

    >>> from pysnptools.util import example_file # Download and return local file name
    >>> # Download the phenotype file if necessary. Return its local location.
    >>> pheno_fn = example_file("pysnptools/examples/toydata.phe")
    >>> print('The local file name is ', pheno_fn)
    The local file name is ...pysnptools/examples/toydata.phe
    >>>
    >>> # Download the bed,bim,&fam files if necessary. Return the location of bed file.
    >>> bedfile = example_file("tests/datasets/all_chr.maf0.001.N300.*","*.bed")
    >>> print('The local file name is ', bedfile)
    The local file name is ...tests/datasets/all_chr.maf0.001.N300.bed

    """
    return pysnptools_hashdown._example_file(pattern, endswith=endswith)


class TestExampleFile(unittest.TestCase):
    def test_doc_test(self):
        import pysnptools.util._example_file as example_mod

        result = doctest.testmod(
            example_mod, optionflags=doctest.ELLIPSIS
        )
        assert result.failed == 0, "failed doc test: " + __file__


def getTestSuite():
    """
    set up test suite
    """

    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestExampleFile))

    return test_suite


if __name__ == "__main__":
    # This creates a Json based on all files, but it may have the wrong line-endings and thus hash on Windows. It can be edited to files of interest
    if True:
        update = Hashdown(url=pysnptools_hashdown.url, file_to_hash=pysnptools_hashdown.file_to_hash, allow_unknown_files=True)
        for file0 in list(os.walk(r'D:\OneDrive\programs\pysnptools\pysnptools\examples'))[0][2:][0]:
            file = 'pysnptools/examples/'+file0
            print(file)
            update.file_exists(file)
        update.save_hashdown("deldir/updated.hashdown.json")
    if False:
        Hashdown.scan_local(
            r"D:\OneDrive\programs\pysnptools",
            url="https://github.com/fastlmm/PySnpTools/raw/cf248cbf762516540470d693532590a77c76fba2",
        ).save_hashdown("deldir/pysnptools.hashdown.json")
    # This fixes the hashes of the files
    if False:
        update = Hashdown(url=pysnptools_hashdown.url, allow_unknown_files=True)
        for file in pysnptools_hashdown.file_to_hash:
            print(file)
            update.file_exists(file)
        update.save_hashdown("deldir/updated.hashdown.json")

    logging.basicConfig(level=logging.INFO)

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=False)
    ret = r.run(suites)
    assert ret.wasSuccessful()

    doctest.testmod(optionflags=doctest.ELLIPSIS)
