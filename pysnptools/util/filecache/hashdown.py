import os
import shutil
import hashlib
import fnmatch

from urllib.request import urlopen, HTTPError
import logging
from contextlib import contextmanager
from pysnptools.util.filecache import FileCache
import tempfile
import json
from pysnptools.util import log_in_place

# Based on an idea see in Danilo Horta's Bgen-reader-py
class Hashdown(FileCache):
    """
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

        >>> from pysnptools.util.filecache import Hashdown
        >>> file_to_hash= {'pysnptools/examples/toydata.5chrom.bed': '766f55aa716bc7bc97cad4de41a50ec3',
        ...                'pysnptools/examples/toydata.5chrom.bim': '6a07f96e521f9a86df7bfd7814eebcd6',
        ...                'pysnptools/examples/toydata.5chrom.fam': 'f4eb01f67e0738d4865fad2014af8537'}
        >>> hashdown = Hashdown('https://github.com/fastlmm/PySnpTools/raw/cf248cbf762516540470d693532590a77c76fba2',
        ...                      file_to_hash=file_to_hash)
        >>> hashdown.file_exists('pysnptools/examples/toydata.5chrom.bed')
        True
        >>> hashdown.load('pysnptools/examples/toydata.5chrom.fam').split('\\n')[0]
        'per0 per0 0 0 2 0.408848'

    """

    def __init__(
        self,
        url,
        file_to_hash={},
        directory=None,
        allow_unknown_files=False,
        trust_local_files=False,
        _relative_directory=None,
    ):
        super(Hashdown, self).__init__()
        self.url = url
        self.file_to_hash = file_to_hash
        self.allow_unknown_files = allow_unknown_files
        self.trust_local_files = trust_local_files
        base_url = (
            url if _relative_directory is None else url[: -len(_relative_directory) - 1]
        )
        url_hash = hashlib.md5(base_url.encode("utf-8")).hexdigest()
        self.directory = tempfile.gettempdir() + "/hashdown/{0}".format(url_hash)
        self._relative_directory = _relative_directory
        if os.path.exists(self.directory):
            assert not os.path.isfile(
                self.directory
            ), "A directory cannot exist where a file already exists."

    def __repr__(self):
        return "{0}('{1}')".format(self.__class__.__name__, self.url)

    @property
    def name(self):
        """
        A path-like name for this `Hashdown`.

        :rtype: string

        >>> from pysnptools.util.filecache import Hashdown
        >>> file_to_hash= {'pysnptools/examples/toydata.5chrom.bed': '766f55aa716bc7bc97cad4de41a50ec3',
        ...                'pysnptools/examples/toydata.5chrom.bim': '6a07f96e521f9a86df7bfd7814eebcd6',
        ...                'pysnptools/examples/toydata.5chrom.fam': 'f4eb01f67e0738d4865fad2014af8537'}
        >>> hashdown = Hashdown('https://github.com/fastlmm/PySnpTools/raw/cf248cbf762516540470d693532590a77c76fba2',
        ...                      file_to_hash=file_to_hash)
        >>> hashdown.name
        'hashdown/9ac30da2bf589db947e91744dff0ec24'
        >>> hashdown.join('pysnptools').name
        'hashdown/9ac30da2bf589db947e91744dff0ec24/pysnptools'

        """
        if self._relative_directory is None:
            url_hash = hashlib.md5(self.url.encode("utf-8")).hexdigest()
            result = "hashdown/{0}".format(url_hash)
        else:
            base_url = self.url[: -len(self._relative_directory) - 1]
            url_hash = hashlib.md5(base_url.encode("utf-8")).hexdigest()
            result = "hashdown/{0}/{1}".format(url_hash, self._relative_directory)

        return result

    @staticmethod
    def _get_large_file(url, file, trust_local_files, length=16 * 1024):
        """https://stackoverflow.com/questions/1517616/stream-large-binary-files-with-urllib2-to-file
        """
        import pysnptools.util as pstutil # put here to avoid recursive nesting

        logging.info("Downloading'{0}'".format(url))
        if trust_local_files and os.path.exists(file):
            logging.info("Trusting local file'{0}'".format(file))
            return True

        try:
            req = urlopen(url)
        except HTTPError as e:
            if e.code == 404:
                return False
            raise

        pstutil.create_directory_if_necessary(file)
        with open(file, "wb") as fp:
            shutil.copyfileobj(req, fp, length)
        return True

    @staticmethod
    def _get_hash(filename):
        """https://stackoverflow.com/questions/16874598/how-do-i-calculate-the-md5-checksum-of-a-file-in-python
        """
        logging.info("Find hash of '{0}'".format(filename))
        with open(filename, "rb") as f:
            file_hash = hashlib.md5()
            while True:
                chunk = f.read(8192)
                if chunk == b"":
                    break
                file_hash.update(chunk)
        return file_hash.hexdigest()

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
            hash = Hashdown._get_hash(full_file)
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
            local_hash = Hashdown._get_hash(full_file)
        else:
            local_hash = None

        if local_hash is None or local_hash != hash:
            assert self._get_large_file(
                full_url, full_file, trust_local_files=False
            ), "URL return 'no item' ('{0}')".format(full_url)
            local_hash = Hashdown._get_hash(full_file)
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
        raise ValueError("Hashdown is read only. writing is not allowed.")

    def _simple_rmtree(self, updater=None):
        raise ValueError('Hashdown is read only. "rmtree" is not allowed.')

    def _simple_remove(self, simple_file_name, updater=None):
        raise ValueError('Hashdown is read only. "remove" is not allowed.')

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
        return Hashdown(
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

    def save_hashdown(self, filename):
        """
        Save a Hashdown object to a json file.

        :param filename: name of file to save to.
        :type path: string

        >>> from pysnptools.util.filecache import Hashdown
        >>> file_to_hash= {'pysnptools/examples/toydata.5chrom.bed': '766f55aa716bc7bc97cad4de41a50ec3',
        ...                'pysnptools/examples/toydata.5chrom.bim': '6a07f96e521f9a86df7bfd7814eebcd6',
        ...                'pysnptools/examples/toydata.5chrom.fam': 'f4eb01f67e0738d4865fad2014af8537'}
        >>> hashdown1 = Hashdown('https://github.com/fastlmm/PySnpTools/raw/cf248cbf762516540470d693532590a77c76fba2',
        ...                      file_to_hash=file_to_hash)
        >>> hashdown1.save_hashdown('tempdir/demo.hashdown.json')
        >>> hashdown2 = Hashdown.load_hashdown('tempdir/demo.hashdown.json')
        >>> hashdown2.file_exists('pysnptools/examples/toydata.5chrom.bed')
        True
        >>> hashdown2.load('pysnptools/examples/toydata.5chrom.fam').split('\\n')[0]
        'per0 per0 0 0 2 0.408848'

        """
        import pysnptools.util as pstutil # put here to avoid recursive nesting

        pstutil.create_directory_if_necessary(filename)
        dict0 = dict(self.__dict__)
        del dict0["directory"]
        del dict0["_relative_directory"]
        del dict0["allow_unknown_files"]
        del dict0["trust_local_files"]
        with open(filename, "w") as json_file:
            json.dump(dict0, json_file)

    @staticmethod
    def load_hashdown(
        filename, directory=None, allow_unknown_files=False, trust_local_files=False
    ):
        """
        Load a Hashdown object from a json file.

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

        :rtype: :class:`.Hashdown`

        >>> from pysnptools.util.filecache import Hashdown
        >>> file_to_hash= {'pysnptools/examples/toydata.5chrom.bed': '766f55aa716bc7bc97cad4de41a50ec3',
        ...                'pysnptools/examples/toydata.5chrom.bim': '6a07f96e521f9a86df7bfd7814eebcd6',
        ...                'pysnptools/examples/toydata.5chrom.fam': 'f4eb01f67e0738d4865fad2014af8537'}
        >>> hashdown1 = Hashdown('https://github.com/fastlmm/PySnpTools/raw/cf248cbf762516540470d693532590a77c76fba2',
        ...                      file_to_hash=file_to_hash)
        >>> hashdown1.save_hashdown('tempdir/demo.hashdown.json')
        >>> hashdown2 = Hashdown.load_hashdown('tempdir/demo.hashdown.json')
        >>> hashdown2.file_exists('pysnptools/examples/toydata.5chrom.bed')
        True
        >>> hashdown2.load('pysnptools/examples/toydata.5chrom.fam').split('\\n')[0]
        'per0 per0 0 0 2 0.408848'

        """
        with open(filename) as json_file:
            dict0 = json.load(json_file)
        hashdown = Hashdown(
            url=dict0["url"],
            file_to_hash=dict0["file_to_hash"],
            directory=directory,
            allow_unknown_files=allow_unknown_files,
            trust_local_files=trust_local_files,
        )
        return hashdown

    @staticmethod
    def scan_local(local_directory, url=None, logging_level=logging.WARNING):
        '''
        Bootstrap a Hashdown by recursively walking a local directory and finding the local MD5 hashes.
        (A local hash might be wrong if the files are out of date or have OS-dependent line endings.)
        Typically, you'll then want to save the result to a JSON file and then edit that JSON file
        manually to remove uninteresting files.

        :param local_directory: Local directory to recursively walk
        :type path: string

        :param url:  URL to give to the Hashdown. (It will not be checked.)
        :type path: string

        :param logging_level: Logging level for printing progress of the walk. Default
               is logging.WARNING)

        :rtype: :class:`.Hashdown`

        '''
        from pysnptools.util.filecache import LocalCache

        file_to_hash = {}
        localcache = LocalCache(local_directory)
        with log_in_place("scanning", logging_level) as updater:
            for file in localcache.walk():
                updater(file)
                with localcache.open_read(file) as full_file:
                    hash = Hashdown._get_hash(full_file)
                    file_to_hash[file] = hash
        return Hashdown(url, file_to_hash=file_to_hash)


if __name__ == "__main__":
    import doctest
    if False:
        hashdown_local = Hashdown.scan_local(
            r"D:\OneDrive\programs\fastlmm",
            url="https://github.com/fastlmm/FaST-LMM/raw/ff183d3aa09c78cf5fdf1961e9241e8a9b9dd172",
        )
        hashdown_local.save_hashdown("deldir/fastlmm.hashdown.json")

    logging.basicConfig(level=logging.INFO)

    doctest.testmod(optionflags=doctest.ELLIPSIS)
