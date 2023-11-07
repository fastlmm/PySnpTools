from io import StringIO
import sys
import unittest
import logging
import os
import tempfile
import shutil
import time
from pysnptools.util.filecache import LocalCache, PeerToPeer, Hashdown


class TestFileCache(unittest.TestCase):
    def test_local_file(self):
        logging.info("test_local_file")

        temp_dir = self._temp_dir()

        def storage_closure():
            return LocalCache(temp_dir)

        self._write_and_read(storage_closure())
        self._distribute(storage_closure)

    def test_hashdown(self):
        logging.info("test_hashdown")

        from pysnptools.util.filecache.test import TestFileCache as self

        logging.basicConfig(level=logging.INFO)
        logging.info("test_hashdown")

        url = "https://github.com/fastlmm/PySnpTools/raw/9de8e93a91b330b064b482c918a38104904b45c0/pysnptools"
        hashdown0 = Hashdown(url, allow_unknown_files=True)
        for file in [
            "examples/toydata.5chrom.bed",
            "examples/toydata.bim",
            "examples/toydata.fam",
            "util/util.py",
        ]:
            hashdown0.file_exists(file)
        _ = hashdown0.file_to_hash
        hashdown = Hashdown(
            url, file_to_hash=hashdown0.file_to_hash, allow_unknown_files=False
        )

        # Clear the directory
        assert self._is_error(
            lambda: hashdown.rmtree()
        )  # raise error because this is read-only
        # Rule: After you clear a directory, nothing is in it
        assert 4 == self._len(
            hashdown.walk()
        )  # returns the files in the file_to_hash_file
        assert not hashdown.file_exists("test.txt")
        assert not hashdown.file_exists("main.txt/test.txt")
        assert not hashdown.file_exists(r"main.txt\test.txt")
        assert self._is_error(
            lambda: hashdown.file_exists("test.txt/")
        )  # Can't query something that can't be a file name
        assert self._is_error(
            lambda: hashdown.file_exists("../test.txt")
        )  # Can't leave the current directory
        if os.name == "nt":
            assert self._is_error(
                lambda: hashdown.file_exists(r"c:\test.txt")
            )  # Can't leave the current directory

        # Rule: '/' and '\' are both OK, but you can't use ':' or '..' to leave the current root.
        assert self._is_error(lambda: 0 == self._len(hashdown.walk("..")))
        assert 0 == self._len(hashdown.walk("..x"))
        assert 0 == self._len(
            hashdown.walk("test.txt")
        )  # This is ok, because test.txt doesn't exist and therefore isn't a file
        assert 0 == self._len(hashdown.walk("a/b"))
        assert 0 == self._len(hashdown.walk("a\\b"))  # Backslash or forward is fine
        assert self._is_error(lambda: len(hashdown.walk("/")))  # Can't start with '/'
        assert self._is_error(lambda: len(hashdown.walk(r"\\")))  # Can't start with '\'
        assert self._is_error(
            lambda: len(hashdown.walk(r"\\computer1\share\3"))
        )  # Can't start with UNC

        # It should be there and be a file
        assert hashdown.file_exists("examples/toydata.bim")
        file_list = list(hashdown.walk())
        assert len(file_list) == 4 and "examples/toydata.bim" in file_list
        file_list2 = list(hashdown.walk("examples"))
        assert len(file_list2) == 3 and "examples/toydata.bim" in file_list2
        assert self._is_error(
            lambda: hashdown.join("examples/toydata.bim")
        )  # Can't create a directory where a file exists
        assert self._is_error(
            lambda: list(hashdown.walk("examples/toydata.bim"))
        )  # Can't have directory where a file exists

        # Read it
        assert (
            hashdown.load("examples/toydata.bim").split("\n")[0]
            == "1\tnull_0\t0\t1\tL\tH"
        )
        assert hashdown.file_exists("examples/toydata.bim")
        assert self._is_error(
            lambda: hashdown.load("examples")
        )  # This is an error because examples is actually a directory and they can't be opened for reading

        # Can query modified time of file.
        assert self._is_error(
            lambda: hashdown.getmtime("a/b/c.txt")
        ), "Can't get mod time from file that doesn't exist"
        assert self._is_error(
            lambda: hashdown.getmtime("examples")
        ), "Can't get mod time from directory"
        assert (
            hashdown.getmtime("examples/toydata.bim") == 0.0
        ), "expect all mod times to be 0"

        # try to write
        assert self._is_error(
            lambda: hashdown.save("main.txt", "")
        ), "Can't write. It's read only"
        # try to remove
        assert self._is_error(
            lambda: hashdown.remove("examples/toydata.bim")
        ), "Can't remove. It's read only"
        assert self._is_error(lambda: hashdown.rmtree()), "Can't remove. It's read only"
        # what if files are not local?
        shutil.rmtree(hashdown.directory, ignore_errors=True)
        # It should be there and be a file
        assert hashdown.file_exists("examples/toydata.bim")
        file_list = list(hashdown.walk())
        assert len(file_list) == 4 and "examples/toydata.bim" in file_list
        file_list2 = list(hashdown.walk("examples"))
        assert len(file_list2) == 3 and "examples/toydata.bim" in file_list2
        assert (
            hashdown.load("examples/toydata.bim").split("\n")[0]
            == "1\tnull_0\t0\t1\tL\tH"
        )
        # What if put bad file locally?
        os.rename(
            hashdown.directory + "/examples/toydata.bim",
            hashdown.directory + "/examples/toydata.fam",
        )
        assert (
            hashdown.load("examples/toydata.fam").split("\n")[0] == "per0 per0 0 0 0 0"
        )
        # What if hash doesn't match web hash?
        hashdown.file_to_hash["examples/toydata.fam"] = "WRONG"
        assert self._is_error(
            lambda: hashdown.load("examples/toydata.fam")
        ), "unexpected hash"
        # writing out file_to_hash and read_file_to_hash

        hashdown.save_hashdown("ignore/toydata.hashdown.json")
        hashdown2 = Hashdown.load_hashdown("ignore/toydata.hashdown.json")
        # It should be there and be a file
        assert hashdown2.file_exists("examples/toydata.bim")
        file_list = list(hashdown2.walk())
        assert len(file_list) == 4 and "examples/toydata.bim" in file_list
        file_list2 = list(hashdown2.walk("examples"))
        assert len(file_list2) == 3 and "examples/toydata.bim" in file_list2
        assert (
            hashdown2.load("examples/toydata.bim").split("\n")[0]
            == "1\tnull_0\t0\t1\tL\tH"
        )

    def test_peer_to_peer(self):
        logging.info("test_peer_to_peer")

        temp_dir = self._temp_dir()

        def storage_closure():
            def id_and_path_function():
                self.count += 1
                return self.count, temp_dir + "/{0}".format(self.count)

            storage = PeerToPeer(
                common_directory=temp_dir + "/common",
                id_and_path_function=id_and_path_function,
            )
            return storage

        self._write_and_read(storage_closure())
        self._distribute(storage_closure)

    @staticmethod
    def file_name(self, testcase_name):
        temp_fn = os.path.join(self.tempout_dir, testcase_name + ".txt")
        if os.path.exists(temp_fn):
            os.remove(temp_fn)
        return temp_fn

    tempout_dir = "tempout/one_milp2p"

    @classmethod
    def setUpClass(self):
        self.temp_parent = os.path.join(tempfile.gettempdir(), format(hash(os.times())))
        self.count = 0

    def _temp_dir(self):
        self.count += 1
        return "{0}/{1}".format(self.temp_parent, self.count)

    @classmethod
    def tearDownClass(self):
        not os.path.exists(self.temp_parent) or shutil.rmtree(self.temp_parent)

    @staticmethod
    def _is_error(lambda0):
        try:
            lambda0()
        except Exception as e:
            logging.debug(str(e))
            return True
        return False

    @staticmethod
    def _len(sequence):
        len = 0
        for item in sequence:
            len += 1
        return len

    def _write_and_read(self, storage):
        test_storage = storage.join("test_snps")
        # !!!How to induce an error: create a 'test_snps' file at the top level then try to create
        # an empty directory with the same name

        # Clear the directory
        test_storage.rmtree()
        # Rule: After you clear a directory, nothing is in it
        assert 0 == self._len(test_storage.walk())
        assert not test_storage.file_exists("test.txt")
        assert not test_storage.file_exists("main.txt/test.txt")
        assert not test_storage.file_exists(r"main.txt\test.txt")
        assert self._is_error(
            lambda: test_storage.file_exists("test.txt/")
        )  # Can't query something that can't be a file name
        assert self._is_error(
            lambda: test_storage.file_exists("../test.txt")
        )  # Can't leave the current directory
        if os.name == "nt":
            assert self._is_error(
                lambda: test_storage.file_exists(r"c:\test.txt")
            )  # Can't leave the current directory

        # Rule: '/' and '\' are both OK, but you can't use ':' or '..' to leave the current root.
        assert 0 == self._len(test_storage.walk())
        assert self._is_error(lambda: 0 == self._len(test_storage.walk("..")))
        assert 0 == self._len(test_storage.walk("..x"))
        assert 0 == self._len(
            test_storage.walk("test.txt")
        )  # This is ok, because test.txt doesn't exist and therefore isn't a file
        assert 0 == self._len(test_storage.walk("a/b"))
        assert 0 == self._len(test_storage.walk("a\\b"))  # Backslash or forward is fine
        assert self._is_error(
            lambda: len(test_storage.walk("/"))
        )  # Can't start with '/'
        assert self._is_error(
            lambda: len(test_storage.walk(r"\\"))
        )  # Can't start with '\'
        assert self._is_error(
            lambda: len(test_storage.walk(r"\\computer1\share\3"))
        )  # Can't start with UNC

        # Clear the directory, again
        test_storage.rmtree()
        assert 0 == self._len(test_storage.walk())
        test_storage.rmtree("main.txt")
        assert 0 == self._len(test_storage.walk("main.txt"))
        assert 0 == self._len(test_storage.walk())

        # Write to it.
        assert self._is_error(lambda: test_storage.save("../test.txt"))
        test_storage.save("main.txt/test.txt", "test\n")
        # Rule: It's an error to write to a file or directory that already exists
        assert self._is_error(lambda: test_storage.save("main.txt"))
        assert self._is_error(lambda: test_storage.save("main.txt/test.txt"))
        assert self._is_error(
            lambda: list(test_storage.walk("main.txt/test.txt"))
        ), "Rule: It's an error to walk a file (but recall that it's OK to walk a folder that doesn't exist)"

        # It should be there and be a file
        assert test_storage.file_exists("main.txt/test.txt")
        file_list = list(test_storage.walk())
        assert len(file_list) == 1 and file_list[0] == "main.txt/test.txt"
        file_list2 = list(test_storage.walk("main.txt"))
        assert len(file_list2) == 1 and file_list2[0] == "main.txt/test.txt"
        assert self._is_error(
            lambda: test_storage.join("main.txt/test.txt")
        )  # Can't create a directory where a file exists
        assert self._is_error(
            lambda: list(test_storage.walk("main.txt/test.txt"))
        )  # Can't create a directory where a file exists
        assert self._is_error(
            lambda: test_storage.rmtree("main.txt/test.txt")
        )  # Can't delete a directory where a file exists

        # Read it
        assert test_storage.load("main.txt/test.txt") == "test\n"
        assert test_storage.file_exists("main.txt/test.txt")
        assert self._is_error(
            lambda: test_storage.load("main.txt")
        )  # This is an error because main.txt is actually a directory and they can't be opened for reading

        # Remove it
        test_storage.remove("main.txt/test.txt")
        assert self._is_error(
            lambda: test_storage.remove("main.txt/test.txt")
        )  # Can't delete a file that doesn't exist
        assert not test_storage.file_exists("main.txt/test.txt")
        assert 0 == self._len(test_storage.walk())
        assert 0 == self._len(test_storage.walk("main.txt"))
        assert 0 == self._len(test_storage.walk("main.txt/test.txt"))  # Now allowed.

        #  writing zero length files is OK
        #  File share has a special file called "main.txt". Can we mess things up by using 'main.txt' as a directory
        #  name, too. It's OK to write to a file in a directory that used to exist, but now has no files.
        test_storage.save("main.txt", "")
        assert test_storage.file_exists("main.txt")
        file_list = list(test_storage.walk())
        assert len(file_list) == 1 and file_list[0] == "main.txt"
        assert test_storage.load("main.txt") == ""
        assert test_storage.file_exists("main.txt")

        # Can query modified time of file. It will be later, later.
        assert self._is_error(
            lambda: test_storage.getmtime("a/b/c.txt")
        ), "Can't get mod time from file that doesn't exist"
        test_storage.save("a/b/c.txt", "")
        m1 = test_storage.getmtime("a/b/c.txt")
        assert self._is_error(
            lambda: test_storage.getmtime("a/b")
        ), "Can't get mod time from directory"
        assert (
            test_storage.getmtime("a/b/c.txt") == m1
        ), "expect mod time to stay the same"
        assert test_storage.load("a/b/c.txt") == ""
        assert (
            test_storage.getmtime("a/b/c.txt") == m1
        ), "reading a file doesn't change its mod time"
        test_storage.remove("a/b/c.txt")
        assert self._is_error(
            lambda: test_storage.getmtime("a/b/c.txt")
        ), "Can't get mod time from file that doesn't exist"
        time.sleep(1)  # Sleep one second
        test_storage.save("a/b/c.txt", "")
        test_storage.save("a/b/d.txt", "")
        assert (
            test_storage.getmtime("a/b/c.txt") > m1
        ), "A file created later (after a pause) will have a later mod time"
        assert (
            test_storage.getmtime("a/b/d.txt") > m1
        ), "A file created later (after a pause) will have a later mod time"
        assert test_storage.getmtime("a/b/d.txt") >= test_storage.getmtime(
            "a/b/c.txt"
        ), "A file created later (with no pause) will have a later or equal mod time"

    # Look for code that isn't covered and make test cases for it
    # Distributed testing
    # What happens if you don't close an open_write? Undefined

    def _distribute(self, storage_lambda):
        storage1 = storage_lambda()
        storage2 = storage_lambda()
        # clear everything on #1
        storage1.rmtree()
        # 1 and #2 agree nothing is there
        assert 0 == self._len(storage1.walk())
        assert 0 == self._len(storage2.walk())
        # write on #1
        storage1.save("a/b/c.txt", "Hello")
        # read on #2
        assert storage2.load("a/b/c.txt") == "Hello"
        # read on #2, again and manually see that it caches
        assert storage2.load("a/b/c.txt") == "Hello"
        # remove on #1
        storage1.remove("a/b/c.txt")
        # assert not exists on #2
        assert not storage2.file_exists("a/b/c.txt")
        # write on #1
        storage1.save("a/b/c.txt", "Hello")
        # read on #2
        assert storage2.load("a/b/c.txt") == "Hello"
        # remove on #1
        storage1.remove("a/b/c.txt")
        # write something different on #1
        storage1.save("a/b/c.txt", "There")
        # read on #2 and see that it is different.
        assert storage2.load("a/b/c.txt") == "There"

    def test_util_filecache_testmod(self):
        import doctest
        import pysnptools.util.filecache

        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        for mod in [
            pysnptools.util.filecache,
            pysnptools.util.filecache.filecache,
            pysnptools.util.filecache.localcache,
            pysnptools.util.filecache.peertopeer,
            pysnptools.util.filecache.hashdown,
        ]:
            buffer = StringIO()
            old_stdout = sys.stdout
            sys.stdout = buffer  # Redirect the stdout to the buffer

            try:
                result = doctest.testmod(
                    mod, optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE
                )
                if result.failed:
                    assert (
                        result.failed == 0
                    ), f"Doc test failed in {mod.__name__} with the following output:\n{buffer.getvalue()}"
            finally:
                sys.stdout = old_stdout
                buffer.close()
        os.chdir(old_dir)


def getTestSuite():
    """
    set up composite test suite
    """

    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestFileCache))
    return test_suite


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARN)
    suites = getTestSuite()

    r = unittest.TextTestRunner(failfast=True)  # !!! should be false
    ret = r.run(suites)
    assert ret.wasSuccessful()

    logging.info("done")
