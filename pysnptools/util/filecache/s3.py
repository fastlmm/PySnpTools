import os
import shutil

from contextlib import contextmanager
import pysnptools.util as pstutil
from pysnptools.util.filecache import FileCache
import tempfile
import json
import boto
import boto.s3.connection
from pysnptools.util import log_in_place

# Based on an idea see in Danilo Horta's Bgen-reader-py
class S3(FileCache):
    """#!!!cmkupdate
    A :class:`.FileCache` for downloading files with known MD5 hashes from a URL. Unloading is not supported.

    See :class:`.FileCache` for general examples of using FileCache.

    **Constructor:**
        :Parameters: * **url** (*string*) -- The URL from which to download any needed files.
                     * **file_to_hash** (*dictionary*) -- A dictionary from file names to MD5 hashes.
                     * **directory** (optional, *string*) -- Local location for files. If not given will be under the system temp directory
                            (typically controlled with the TEMP environment variable).
                     * **allow_unknown_files** (optional, *bool*) -- By default, all requested files must be in the dictionary. If True,
                            other files can be requested. If found under the URL, they will be downloaded and an entry will be added
                            to the dictionary.
                     * **trust_local_files** (optional, bool) -- By default, when **allow_unknown_files** is True, unknown files
                            will be download. If **trust_local_files** is also True, then any local files in **directory** will
                            be assumed to have the correct hash.
                     * **_relative_directory** (*string*) -- Used internally

        :Example:

        >>> from pysnptools.util.filecache import S3
        >>> file_to_hash= {'pysnptools/examples/toydata.5chrom.bed': '766f55aa716bc7bc97cad4de41a50ec3',
        ...                'pysnptools/examples/toydata.5chrom.bim': '6a07f96e521f9a86df7bfd7814eebcd6',
        ...                'pysnptools/examples/toydata.5chrom.fam': 'f4eb01f67e0738d4865fad2014af8537'}
        >>> S3 = S3('https://github.com/fastlmm/PySnpTools/raw/cf248cbf762516540470d693532590a77c76fba2',
        ...                      file_to_hash=file_to_hash)
        >>> S3.file_exists('pysnptools/examples/toydata.5chrom.bed')
        True
        >>> S3.load('pysnptools/examples/toydata.5chrom.fam').split('\\n')[0]
        'per0 per0 0 0 2 0.408848'

    """

    def __init__(
        self,
        s3_root,
        local_root=None
        credentials='~/.aws/credentials',
        subpath=".",#!!!cmk should hashdown rename to this?
    ):
        super(S3, self).__init__()
        self.s3_root = s3_root
        self.local_root = local_root
        self.credentials = credentials
        self.subpath = subpath

    def __repr__(self):
        return "{0}('{1}')".format(self.__class__.__name__, self.name)

    @property
    def name(self):
        """
        A path-like name for this `S3`.

        :rtype: string

        >>> from pysnptools.util.filecache import S3
        >>> file_to_hash= {'pysnptools/examples/toydata.5chrom.bed': '766f55aa716bc7bc97cad4de41a50ec3',
        ...                'pysnptools/examples/toydata.5chrom.bim': '6a07f96e521f9a86df7bfd7814eebcd6',
        ...                'pysnptools/examples/toydata.5chrom.fam': 'f4eb01f67e0738d4865fad2014af8537'}
        >>> S3 = S3('https://github.com/fastlmm/PySnpTools/raw/cf248cbf762516540470d693532590a77c76fba2',
        ...                      file_to_hash=file_to_hash)
        >>> S3.name
        'S3/9ac30da2bf589db947e91744dff0ec24'
        >>> S3.join('pysnptools').name
        'S3/9ac30da2bf589db947e91744dff0ec24/pysnptools'

        """
        return self.s3_root #!!!cmk should include subpath?

    def _simple_file_exists(self, simple_file_name):
        rel_part = (
            "" if self._relative_directory is None else self._relative_directory + "/"
        )
        full_file = self.directory + "/" + rel_part + simple_file_name
        relative_file = (
            simple_file_name
            if self._relative_directory is None
            else self._relative_directory + "/" + simple_file_name
        )
        full_url = self.url + "/" + simple_file_name
        hash = self.file_to_hash.get(relative_file)

        if hash is not None:
            return True

        if not self.allow_unknown_files:
            return False

        if self._get_large_file(full_url, full_file, self.trust_local_files):
            hash = self._get_hash(full_file)
            self.file_to_hash[relative_file] = hash
            return True
        else:
            return False

    def _example_file(self, pattern, endswith=None):
        return_file = None
        for filename in fnmatch.filter(self.file_to_hash, pattern):
            with self.open_read(filename) as local_file:
                if return_file is None and (
                    endswith is None or fnmatch.fnmatch(filename, endswith)
                ):
                    return_file = local_file
        assert return_file is not None, "Pattern not found '{0}'{1}".format(
            pattern, "" if (endswith is None) else "'{0}'".format(endswith)
        )
        return return_file

    @contextmanager
    def _simple_open_read(self, simple_file_name, updater=None):
        logging.info("open_read('{0}')".format(simple_file_name))

        relative_file = (
            simple_file_name
            if self._relative_directory is None
            else self._relative_directory + "/" + simple_file_name
        )
        full_file = self.directory + "/" + relative_file
        full_url = self.url + "/" + simple_file_name
        hash = self.file_to_hash.get(relative_file)
        assert self._simple_file_exists(
            simple_file_name
        ), "File doesn't exist ('{0}')".format(relative_file)
        if os.path.exists(full_file):
            local_hash = self._get_hash(full_file)
        else:
            local_hash = None

        if local_hash is None or local_hash != hash:
            assert self._get_large_file(
                full_url, full_file, trust_local_files=False
            ), "URL return 'no item' ('{0}')".format(full_url)
            local_hash = self._get_hash(full_file)
            if hash is None:
                assert self.allow_unknown_files, "real assert"
                self.file_to_hash[relative_file] = local_hash
            else:
                assert (
                    hash == local_hash
                ), 'URL file has unexpected hash ("{0}")'.format(full_url)

        yield full_file

        logging.info("close('{0}')".format(simple_file_name))

    @contextmanager
    def _simple_open_write(self, simple_file_name, size=0, updater=None):
        raise ValueError("S3 is read only. writing is not allowed.")

    def _simple_rmtree(self, updater=None):
        raise ValueError('S3 is read only. "rmtree" is not allowed.')

    def _simple_remove(self, simple_file_name, updater=None):
        raise ValueError('S3 is read only. "remove" is not allowed.')

    def _simple_getmtime(self, simple_file_name):
        assert self._simple_file_exists(
            simple_file_name
        ), "file doesn't exist ('{0}')".format(simple_file_name)
        return 0

    def _simple_join(self, path):
        directory = self.directory + "/" + path
        _relative_directory = (
            path
            if self._relative_directory is None
            else self._relative_directory + "/" + path
        )
        if not self.allow_unknown_files:
            assert not self.file_exists(
                _relative_directory
            ), "Can't treat an existing file as a directory"
        return S3(
            url=self.url + "/" + path,
            directory=directory,
            file_to_hash=self.file_to_hash,
            allow_unknown_files=self.allow_unknown_files,
            _relative_directory=_relative_directory,
        )

    def _simple_walk(self):
        for rel_file in self.file_to_hash:
            if self._relative_directory is None or rel_file.startswith(
                self._relative_directory + "/"
            ):
                file = (
                    rel_file
                    if self._relative_directory is None
                    else rel_file[len(self._relative_directory) + 1:]
                )
                if self.file_exists(file):
                    yield file

    def save_S3(self, filename):
        """
        Save a S3 object to a json file.

        :param filename: name of file to save to.
        :type path: string

        >>> from pysnptools.util.filecache import S3
        >>> file_to_hash= {'pysnptools/examples/toydata.5chrom.bed': '766f55aa716bc7bc97cad4de41a50ec3',
        ...                'pysnptools/examples/toydata.5chrom.bim': '6a07f96e521f9a86df7bfd7814eebcd6',
        ...                'pysnptools/examples/toydata.5chrom.fam': 'f4eb01f67e0738d4865fad2014af8537'}
        >>> S31 = S3('https://github.com/fastlmm/PySnpTools/raw/cf248cbf762516540470d693532590a77c76fba2',
        ...                      file_to_hash=file_to_hash)
        >>> S31.save_S3('tempdir/demo.S3.json')
        >>> S32 = S3.load_S3('tempdir/demo.S3.json')
        >>> S32.file_exists('pysnptools/examples/toydata.5chrom.bed')
        True
        >>> S32.load('pysnptools/examples/toydata.5chrom.fam').split('\\n')[0]
        'per0 per0 0 0 2 0.408848'

        """
        pstutil.create_directory_if_necessary(filename)
        dict0 = dict(self.__dict__)
        del dict0["directory"]
        del dict0["_relative_directory"]
        del dict0["allow_unknown_files"]
        del dict0["trust_local_files"]
        with open(filename, "w") as json_file:
            json.dump(dict0, json_file)

    @staticmethod
    def load_S3(
        filename, directory=None, allow_unknown_files=False, trust_local_files=False
    ):
        """
        Load a S3 object from a json file.

        :param filename: name of file to load from.
        :type path: string

        :param directory:  Local location for files. If not given will be under the system temp directory
                           (typically controlled with the TEMP environment variable).
        :type path: string

        :param allow_unknown_files: By default, all requested files must be in the dictionary. If True,
                                    other files can be requested. If found under the URL, they will be downloaded and an entry will be added
                                    to the dictionary.
        :type path: bool

        :param trust_local_files: By default, when **allow_unknown_files** is True, unknown files
                                  will be download. If **trust_local_files** is also True, then any local files in **directory** will
                                  be assumed to have the correct hash.
        :type path: bool

        :rtype: :class:`.S3`

        >>> from pysnptools.util.filecache import S3
        >>> file_to_hash= {'pysnptools/examples/toydata.5chrom.bed': '766f55aa716bc7bc97cad4de41a50ec3',
        ...                'pysnptools/examples/toydata.5chrom.bim': '6a07f96e521f9a86df7bfd7814eebcd6',
        ...                'pysnptools/examples/toydata.5chrom.fam': 'f4eb01f67e0738d4865fad2014af8537'}
        >>> S31 = S3('https://github.com/fastlmm/PySnpTools/raw/cf248cbf762516540470d693532590a77c76fba2',
        ...                      file_to_hash=file_to_hash)
        >>> S31.save_S3('tempdir/demo.S3.json')
        >>> S32 = S3.load_S3('tempdir/demo.S3.json')
        >>> S32.file_exists('pysnptools/examples/toydata.5chrom.bed')
        True
        >>> S32.load('pysnptools/examples/toydata.5chrom.fam').split('\\n')[0]
        'per0 per0 0 0 2 0.408848'

        """
        with open(filename) as json_file:
            dict0 = json.load(json_file)
        S3 = S3(
            url=dict0["url"],
            file_to_hash=dict0["file_to_hash"],
            directory=directory,
            allow_unknown_files=allow_unknown_files,
            trust_local_files=trust_local_files,
        )
        return S3

    @staticmethod
    def scan_local(local_directory, url=None, logging_level=logging.WARNING):
        '''
        Bootstrap a S3 by recursively walking a local directory and finding the local MD5 hashes.
        (A local hash might be wrong if the files are out of date or have OS-dependent line endings.)
        Typically, you'll then want to save the result to a JSON file and then edit that JSON file
        manually to remove uninteresting files.

        :param local_directory: Local directory to recursively walk
        :type path: string

        :param url:  URL to give to the S3. (It will not be checked.)
        :type path: string

        :param logging_level: Logging level for printing progress of the walk. Default
               is logging.WARNING)

        :rtype: :class:`.S3`

        '''
        from pysnptools.util.filecache import LocalCache

        file_to_hash = {}
        localcache = LocalCache(local_directory)
        with log_in_place("scanning", logging_level) as updater:
            for file in localcache.walk():
                updater(file)
                with localcache.open_read(file) as full_file:
                    hash = S3._get_hash(full_file)
                    file_to_hash[file] = hash
        return S3(url, file_to_hash=file_to_hash)


if __name__ == "__main__":
    import doctest
    logging.basicConfig(level=logging.INFO)

    if True:
        storage = S3(bucket='traviscibucket',credentials='~/.aws/credentials')
        test_storage = storage.join('test_snps') #!!!How to induce an error: create a 'test_snps' file at the top level then try to create an empty directory with the same name


        #Clear the directory
        test_storage.rmtree()
        #Rule: After you clear a directory, nothing is in it
        assert 0 == self._len(test_storage.walk())
        assert not test_storage.file_exists("test.txt")
        assert not test_storage.file_exists("main.txt/test.txt")
        assert not test_storage.file_exists(r"main.txt\test.txt")
        assert self._is_error(lambda : test_storage.file_exists("test.txt/")) #Can't query something that can't be a file name
        assert self._is_error(lambda : test_storage.file_exists("../test.txt")) #Can't leave the current directory
        if os.name == 'nt':
            assert self._is_error(lambda : test_storage.file_exists(r"c:\test.txt")) #Can't leave the current directory

        #Rule: '/' and '\' are both OK, but you can't use ':' or '..' to leave the current root.
        assert 0 == self._len(test_storage.walk())
        assert self._is_error(lambda : 0 == self._len(test_storage.walk("..")))
        assert 0 == self._len(test_storage.walk("..x"))
        assert 0 == self._len(test_storage.walk("test.txt")) #This is ok, because test.txt doesn't exist and therefore isn't a file
        assert 0 == self._len(test_storage.walk("a/b"))
        assert 0 == self._len(test_storage.walk("a\\b")) #Backslash or forward is fine
        assert self._is_error(lambda : len(test_storage.walk("/"))) #Can't start with '/'
        assert self._is_error(lambda : len(test_storage.walk(r"\\"))) #Can't start with '\'
        assert self._is_error(lambda : len(test_storage.walk(r"\\computer1\share\3"))) #Can't start with UNC

        #Clear the directory, again
        test_storage.rmtree()
        assert 0 == self._len(test_storage.walk())
        test_storage.rmtree("main.txt")
        assert 0 == self._len(test_storage.walk("main.txt"))
        assert 0 == self._len(test_storage.walk())


        #Write to it.
        assert self._is_error(lambda : test_storage.save("../test.txt"))
        test_storage.save("main.txt/test.txt","test\n")
        #Rule: It's an error to write to a file or directory that already exists
        assert self._is_error(lambda : test_storage.save("main.txt")) 
        assert self._is_error(lambda : test_storage.save("main.txt/test.txt")) 
        assert self._is_error(lambda : list(test_storage.walk("main.txt/test.txt"))), "Rule: It's an error to walk a file (but recall that it's OK to walk a folder that doesn't exist)"

        #It should be there and be a file
        assert test_storage.file_exists("main.txt/test.txt")
        file_list = list(test_storage.walk())
        assert len(file_list)==1 and file_list[0] == "main.txt/test.txt"
        file_list2 = list(test_storage.walk("main.txt"))
        assert len(file_list2)==1 and file_list2[0] == "main.txt/test.txt"
        assert self._is_error(lambda : test_storage.join("main.txt/test.txt")) #Can't create a directory where a file exists
        assert self._is_error(lambda : list(test_storage.walk("main.txt/test.txt"))) #Can't create a directory where a file exists
        assert self._is_error(lambda : test_storage.rmtree("main.txt/test.txt")) #Can't delete a directory where a file exists

        #Read it
        assert test_storage.load("main.txt/test.txt")=="test\n"
        assert test_storage.file_exists("main.txt/test.txt")
        assert self._is_error(lambda : test_storage.load("main.txt"))  #This is an error because main.txt is actually a directory and they can't be opened for reading

        #Remove it
        test_storage.remove("main.txt/test.txt")
        assert self._is_error(lambda : test_storage.remove("main.txt/test.txt")) #Can't delete a file that doesn't exist
        assert not test_storage.file_exists("main.txt/test.txt")
        assert 0 == self._len(test_storage.walk())
        assert 0 == self._len(test_storage.walk("main.txt"))
        assert 0 == self._len(test_storage.walk("main.txt/test.txt")) #Now allowed.

        #  writing zero length files is OK
        #  File share has a special file called "main.txt". Can we mess things up by using 'main.txt' as a directory name, too.
        #It's OK to write to a file in a directory that used to exist, but now has no files.
        test_storage.save("main.txt","")
        assert test_storage.file_exists("main.txt")
        file_list = list(test_storage.walk())
        assert len(file_list)==1 and file_list[0] == "main.txt"
        assert test_storage.load("main.txt") == ""
        assert test_storage.file_exists("main.txt")

        #Can query modified time of file. It will be later, later.
        assert self._is_error(lambda : test_storage.getmtime("a/b/c.txt")), "Can't get mod time from file that doesn't exist"
        test_storage.save("a/b/c.txt","")
        m1 = test_storage.getmtime("a/b/c.txt")
        assert self._is_error(lambda : test_storage.getmtime("a/b")), "Can't get mod time from directory"
        assert test_storage.getmtime("a/b/c.txt") == m1, "expect mod time to stay the same"
        assert test_storage.load("a/b/c.txt") == ""
        assert test_storage.getmtime("a/b/c.txt") == m1, "reading a file doesn't change its mod time"
        test_storage.remove("a/b/c.txt")
        assert self._is_error(lambda : test_storage.getmtime("a/b/c.txt")), "Can't get mod time from file that doesn't exist"
        time.sleep(1) #Sleep one second
        test_storage.save("a/b/c.txt","")
        test_storage.save("a/b/d.txt","")
        assert test_storage.getmtime("a/b/c.txt") > m1, "A file created later (after a pause) will have a later mod time"
        assert test_storage.getmtime("a/b/d.txt") > m1, "A file created later (after a pause) will have a later mod time"
        assert test_storage.getmtime("a/b/d.txt") >= test_storage.getmtime("a/b/c.txt"), "A file created later (with no pause) will have a later or equal mod time"

        logging.info("done")



    doctest.testmod(optionflags=doctest.ELLIPSIS)
