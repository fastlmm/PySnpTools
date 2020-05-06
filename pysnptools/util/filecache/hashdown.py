import os
import shutil
import hashlib
try:
    from urllib.request import urlopen, HTTPError # Python 3
except ImportError:
    from urllib2 import urlopen, HTTPError # Python 2
import logging
from contextlib import contextmanager
import pysnptools.util as pstutil
from pysnptools.util.filecache import FileCache
import tempfile

class Hashdown(FileCache):
    '''
    A :class:`.FileCache` for working with files locally.#!!!cmk update

    See :class:`.FileCache` for general examples of using FileCache.

    This is the simplest :class:`.FileCache` in that it stores everything on a local disk rather than storing things remotely.

    **Constructor:**
        :Parameters: * **directory** (*string*) -- The directory under which files should be written and read.

        :Example:

        >>> from pysnptools.util.filecache import Hashdown
        >>> file_cache = Hashdown('localcache1')
        >>> file_cache.rmtree()
        >>> file_cache.save('sub1/file1.txt','Hello')
        >>> file_cache.file_exists('sub1/file1.txt')
        True

    '''
    def __init__(self,url,file_to_hash={},directory=None,allow_unhashed_files=False,relative_dir=None):
        super(Hashdown, self).__init__()
        self.url = url
        self.file_to_hash = file_to_hash
        self.allow_unhashed_files = allow_unhashed_files
        self.directory =  tempfile.gettempdir() + '/hashdown/{0}'.format(hash(url))#!!!cmk0 use md5
        self.relative_dir = relative_dir
        if os.path.exists(self.directory): assert not os.path.isfile(self.directory), "A directory cannot exist where a file already exists."


    def __repr__(self): 
        return "{0}('{1}')".format(self.__class__.__name__,self.url)

    @property
    def name(self):
        '''
        A path-like name for this `LocalCache`.#!!!update

        :rtype: string
        '''
        return self.directory

    @staticmethod
    def _get_large_file(url, file, length=16*1024):
        '''https://stackoverflow.com/questions/1517616/stream-large-binary-files-with-urllib2-to-file
        '''
        logging.info("Downloading'{0}'".format(url))
        try:
            req = urlopen(url)
        except HTTPError as e:
            if e.code == 404:
                return False
            raise

        pstutil.create_directory_if_necessary(file)
        with open(file, 'wb') as fp:
            shutil.copyfileobj(req, fp, length)
        return True

    @staticmethod
    def _get_hash(filename):
        '''https://stackoverflow.com/questions/16874598/how-do-i-calculate-the-md5-checksum-of-a-file-in-python
        '''
        logging.info("Find hash of '{0}'".format(filename))
        with open(filename, "rb") as f:
            file_hash = hashlib.md5()
            while True:
                chunk = f.read(8192)
                if chunk is b'':
                    break
                file_hash.update(chunk)
        return file_hash.hexdigest()

    def _simple_file_exists(self,simple_file_name):
        full_file = self.directory + "/" + simple_file_name
        relative_file = simple_file_name if self.relative_dir is None else self.relative_dir + '/' + simple_file_name
        full_url = self.url + "/" + simple_file_name
        hash = self.file_to_hash.get(relative_file)

        if hash is not None:
            return True

        if not self.allow_unhashed_files:
            return False

        if self._get_large_file(full_url,full_file): #!!!cmk this should overwrite any local file. #!!!cmk when did directory get created?
            hash = self._get_hash(full_file)
            self.file_to_hash[relative_file]=hash
            return True
        else:
            return False

            #if os.path.exists(full_file):
            #    local_hash = self._get_hash(full_file)
            #    if local_hash == hash:
            #        return True
            #    else:
            #        logging.info('Local file has wrong hash. Will try to download ("{0}")'.format(relative_file))
            #        if not self._get_large_file(full_url,full_file): #!!!cmk this should overwrite any local file. #!!!cmk when did directory get created?
            #            return False
            #        local_hash = self._get_hash(full_file)
            #        if hash!=local_hash:
            #            logging.warn('URL file has unexpected hash ("{0}")'.format(full_url))
            #            return False
            #        return True
                        

    @contextmanager
    def _simple_open_read(self,simple_file_name,updater=None):
        logging.info("open_read('{0}')".format(simple_file_name))

        full_file = self.directory + "/" + simple_file_name
        relative_file = simple_file_name if self.relative_dir is None else self.relative_dir + '/' + simple_file_name
        full_url = self.url + "/" + simple_file_name
        hash = self.file_to_hash.get(relative_file)
        assert self._simple_file_exists(simple_file_name), "File doesn't exist ('{0}')".format(relative_file)
        if os.path.exists(full_file):
            local_hash = self._get_hash(full_file)
        else:
            local_hash = None

        if local_hash is None or local_hash != hash:
            assert self._get_large_file(full_url,full_file), "URL return 'no item' ('{0}')".format(full_url)
            local_hash = self._get_hash(full_file)
            assert hash==local_hash, 'URL file has unexpected hash ("{0}")'.format(full_url)

        yield full_file

        logging.info("close('{0}')".format(simple_file_name))

    @contextmanager
    def _simple_open_write(self,simple_file_name,size=0,updater=None):
        logging.info("open_write('{0}',size={1})".format(simple_file_name,size))

        #Register the file name in the directory
        file_name=self.directory + "/" + simple_file_name
        if os.path.exists(file_name): #This is only OK, if it is a directory containing no files (and thus doesn't exist)
            assert not os.path.isfile(file_name), "Can't open a file for write if it already exists."
            assert not self._at_least_one(self.walk(simple_file_name)), "Can't open a file for write if a directory with files already has the same name ({0},{1})".format(self,simple_file_name)
            shutil.rmtree(file_name)
        else:
            pstutil.create_directory_if_necessary(file_name,isfile=True)

        yield file_name

        logging.info("close('{0}')".format(simple_file_name))
        assert os.path.exists(file_name), "File doesn't exist in LocalCache. File is '{0}'".format(file_name)

    def _simple_rmtree(self,updater=None):
        raise ValueError('Hashdown is read only. "rmtree" is not allowed.')

    def _simple_remove(self,simple_file_name, updater=None):
        raise ValueError('Hashdown is read only. "remove" is not allowed.')

    def _simple_getmtime(self,simple_file_name):
        full_file = self.directory + "/" + simple_file_name
        assert os.path.exists(full_file) and os.path.isfile(full_file), "Expect file to exist (and be a file)"
        return os.path.getmtime(full_file)

    def _simple_join(self,path):
        directory = self.directory + "/" + path
        relative_dir = path if self.relative_dir is None else self.relative_dir + '/' + path
        if not self.allow_unhashed_files:
            assert not self.file_exists(relative_dir), "Can't treat an existing file as a directory"
        return Hashdown(url=self.url+'/'+path,directory=directory,file_to_hash=self.file_to_hash,allow_unhashed_files=self.allow_unhashed_files,
                        relative_dir=relative_dir)

            
            

    def _simple_walk(self):
        for rel_file in self.file_to_hash:
            if self.relative_dir is None or rel_file.startswith(self.relative_dir+"/"):
                file = rel_file if self.relative_dir is None else rel_file[len(self.relative_dir)+1:]
                if self.file_exists(file):
                    yield file



if __name__ == "__main__":
    if True:
        from pysnptools.util.filecache.test import TestFileCache as self
        logging.basicConfig(level=logging.INFO)
        logging.info("test_hashdown")

        url='https://github.com/fastlmm/PySnpTools/tree/9de8e93a91b330b064b482c918a38104904b45c0/pysnptools'
        hashdown_setup = Hashdown(url,
                            #file_to_hash_file ='' #!!!cmk a text file
                            #!!!cmk option to specify folder -- otherwise temp + hash of URL
                            #!!!cmk don't let it go above local folder
                            #!!!cmk do support folders
                            allow_unhashed_files=True #!!!cmk .create_file_to_hash_file
                            )
        for file in ['examples/toydata.bed','examples/toydata.bim','examples/toydata.fam','util/util.py']:
            hashdown_setup.file_exists(file)
        file_to_hash = hashdown_setup.file_to_hash

        hashdown_setup = Hashdown(url, file_to_hash = file_to_hash)

        hashdown = Hashdown(url, file_to_hash=hashdown_setup.file_to_hash, allow_unhashed_files=False)

        #Clear the directory
        assert self._is_error(lambda : hashdown.rmtree()) #!!!cmk raise error because this is read-only
        #Rule: After you clear a directory, nothing is in it
        assert 4 == self._len(hashdown.walk()) #returns the files in the file_to_hash_file
        assert not hashdown.file_exists("test.txt")
        assert not hashdown.file_exists("main.txt/test.txt")
        assert not hashdown.file_exists(r"main.txt\test.txt")
        assert self._is_error(lambda : hashdown.file_exists("test.txt/")) #Can't query something that can't be a file name
        assert self._is_error(lambda : hashdown.file_exists("../test.txt")) #Can't leave the current directory
        if os.name == 'nt':
            assert self._is_error(lambda : hashdown.file_exists(r"c:\test.txt")) #Can't leave the current directory

        #Rule: '/' and '\' are both OK, but you can't use ':' or '..' to leave the current root.
        assert self._is_error(lambda : 0 == self._len(hashdown.walk("..")))
        assert 0 == self._len(hashdown.walk("..x"))
        assert 0 == self._len(hashdown.walk("test.txt")) #This is ok, because test.txt doesn't exist and therefore isn't a file
        assert 0 == self._len(hashdown.walk("a/b"))
        assert 0 == self._len(hashdown.walk("a\\b")) #Backslash or forward is fine
        assert self._is_error(lambda : len(hashdown.walk("/"))) #Can't start with '/'
        assert self._is_error(lambda : len(hashdown.walk(r"\\"))) #Can't start with '\'
        assert self._is_error(lambda : len(hashdown.walk(r"\\computer1\share\3"))) #Can't start with UNC


        #It should be there and be a file
        assert hashdown.file_exists('examples/toydata.bim')
        file_list = list(hashdown.walk())
        assert len(file_list)==4 and 'examples/toydata.bim' in file_list
        file_list2 = list(hashdown.walk("examples"))
        assert len(file_list2)==3 and 'examples/toydata.bim' in file_list2
        assert self._is_error(lambda : hashdown.join('examples/toydata.bim')) #Can't create a directory where a file exists
        assert self._is_error(lambda : list(hashdown.walk('examples/toydata.bim'))) #Can't have directory where a file exists

        #Read it
        assert hashdown.load('examples/toydata.bim').split('\n')[0] =="1\tnull_0\t0\t1\tL\tH\n"
        assert hashdown.file_exists('examples/toydata.bim')
        assert self._is_error(lambda : hashdown.load("examples"))  #This is an error because examples is actually a directory and they can't be opened for reading


        #Can query modified time of file. It will be later, later.
        assert self._is_error(lambda : hashdown.getmtime("a/b/c.txt")), "Can't get mod time from file that doesn't exist"
        assert self._is_error(lambda : hashdown.getmtime("examples")), "Can't get mod time from directory"
        assert hashdown.getmtime('examples/toydata.bim') == hashdown.fixed_mtime, "expect all mod times to be the same"


    logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()