#cmk
#import os
#import logging
#from pysnptools.util.filecache import FileCache
#from contextlib import contextmanager
#from pysnptools.util.filecache import ip_address
#import smart_open
##cmk pip install smart_open[s3] 
##cmk import boto3

#def _path_join(*p): #!!!cmk understand this and see if needed
#    result = os.path.normpath("/".join(p)).replace('\\','/')
#    return result


#class S3so(FileCache):

#    '''
#    A class that subclasses :class:`FileCache` storage on Amazon S3 Storage.

#    **Constructor:**
#        :Parameters: * **folder** (*string*) -- The path on S3 storage under which data will be stored. The form of the path is:
#                           /BUCKET/morepath.
#                     * **id_and_path_function** (*a zero-augment lambda*) -- When called, tells were to store data locally. See cmk :func:`ip_address_local`.

#        cmk tell about creds

#    **All the methods of FileCache plus these:**
#    '''

#    def __init__(self, folder, id_and_path_function=lambda:(None,"."), process_count=None,upload_local=False):
#        super(S3so, self).__init__()
#        self.folder = folder
#        self.id_and_path_function = id_and_path_function
#        self.s3so_file_system = smart_open.S3FileSystem() #cmk do something with process_count
#        self._upload_local = upload_local #cmk?

#    def __repr__(self):
#        return "{0}('{1}')".format(self.__class__.__name__,self.name)

#    @property
#    def name(self):
#        return self.folder

#    def _simple_file_exists(self,simple_file_name,updater=None):
#        local = self._get_local(simple_file_name)
#        remote = self._get_remote(simple_file_name)
#        if self._upload_local and os.path.exists(local): #!!!cmk understand this code path
#            if not self.s3so_file_system.file_exists(remote):
#                self.s3so_file_system.upload(local,remote,do_sync_date=True,updater=updater)
#            return True
#        else:
#            return self.s3so_file_system.exists(remote) and self.s3so_file_system.isfile(remote)

#    def _get_remote(self,file_name):
#            return self.folder + "/" + file_name

#    def _get_local(self,file_name):
#        local = _path_join(self.id_and_path_function()[1], self.folder, file_name)
#        return local

#    @contextmanager
#    def _simple_open_read(self,simple_file_name,updater=None):
#        local = self._get_local(simple_file_name)
#        remote = self._get_remote(simple_file_name)
#        if self._upload_local and os.path.exists(local) and not self.s3so_file_system.exists(remote):
#            size = os.path.getsize(local)
#            with progress_reporter("Uploading local before read: "+os.path.basename(local), size, level=logging.INFO) as updater2:
#                self.s3so_file_system.upload(local,remote,do_sync_date=True,updater=updater2)
#        self.s3so_file_system.get(remote,local)#cmk,as_needed=True,updater=updater)
#        yield local

#    @contextmanager
#    def _simple_open_write(self,simple_file_name,size=0,updater=None):

#        assert not self._simple_file_exists(simple_file_name), "Can't open a file for write if it already exists. ({0},'{1}')".format(self,simple_file_name)
#        assert not self._at_least_one(self.walk(simple_file_name)), "Can't open a file for write if a directory with files already has the same name ({0},{1})".format(self,simple_file_name)

#        #later could remove other local files to make room for this one
#        local=self._get_local(simple_file_name)
#        #Anything in storage (file or directory) can be cleaned up.
#        self._create_directory(local)
#        remote = self._get_remote(simple_file_name)

#        yield local

#        self.s3so_file_system.put_file(local,remote)#!!!cmk,do_sync_date=True,updater=updater)

#    def _simple_remove(self,simple_file_name,updater=None):
#        assert self._simple_file_exists(simple_file_name), "Expect file to exist (and be a file) ({0},'{1}')".format(self,simple_file_name)
#        self.s3so_file_system.rm(self.folder + "/" + simple_file_name)#!!!cmk,updater=updater)
#        if self._upload_local:
#           local=self._get_local(simple_file_name)
#           if os.path.exists(local):
#                os.remove(local)

#    def _simple_rmtree(self,updater=None):
#        if self.s3so_file_system.exists(self.folder):
#            self.s3so_file_system.rm(self.folder,recursive=True)

#    def _simple_join(self,path):
#        assert not self.s3so_file_system.exists(self.folder + "/" + path) or self.s3so_file_system.isdir(self.folder + "/" + path), "Can't treat an existing file as a directory: '{0}'".format(self.folder + "/" + path)
#        return S3so(self.folder + "/" + path, id_and_path_function=self.id_and_path_function,upload_local=self._upload_local)

        
#    def _simple_walk(self):
#        for full in self.s3so_file_system.find(self.folder):
#            rel = os.path.relpath("/"+full,self.folder).replace('\\','/')
#            yield rel

#    def _simple_getmtime(self,simple_file_name):
#        assert self._simple_file_exists(simple_file_name),"Can only get mtime from files"
#        return self.s3so_file_system.modified(self.folder + "/" + simple_file_name).timestamp()



#if __name__ == '__main__':
#    logging.basicConfig(level=logging.INFO)

#    if True:
#        # Created S3so storage and got key and secret (could access with S3so 3rd party tool)
#        # Create ~/.aws/credentials as per https://stackoverflow.com/questions/33297172/boto3-error-botocore-exceptions-nocredentialserror-unable-to-locate-credential
#        import smart_open

#        with smart_open.open()

#        fs = s3fs.S3FileSystem(anon=False) 
#        #/traviscibucket/fastlmm/PySnpTools/master/136.bd53dc96f9213bba622f5cf21677c09b5c544037/win_py38/
#        fs.ls("/traviscibucket/fastlmm/PySnpTools/master/136.bd53dc96f9213bba622f5cf21677c09b5c544037/win_py38/pysnptools-0.4.21.tar.gz")

#    if False:
#        from pysnptools.util.filecache.xfsspec import S3so, ip_address
#        def id_and_path_function():
#            ip = ip_address()
#            return ip, 's3/{0}'.format(ip)
#        file_cache = S3so(folder='/traviscibucket/deldir',id_and_path_function=id_and_path_function)
#        print(file_cache)
#        #-> PeerToPeer('peertopeer1/common',id_and_path_function=...')

#        # Remove anything already in remote storage
#        file_cache.rmtree()
        
#        # Create a random SNP file and store it remotely.
#        from pysnptools.snpreader import SnpGen, Dense
#        snp_gen = SnpGen(seed=123,iid_count=100,sid_count=500)
#        with file_cache.open_write('r123.100x500.dense.txt') as local_filename:
#            dense1 = Dense.write(local_filename,snp_gen.read())
#        list(file_cache.walk())
#        #-> ['r123.100x500.dense.txt']

#        #Copy back from remote storage (if needed) and then read SNP data from local file.

#        with file_cache.open_read('r123.100x500.dense.txt') as local_filename:
#             dense2 = Dense(local_filename)
#             print(dense2[:3,:3].read().val) #Read 1st 3 individuals and SNPs
#        #[[ 0. nan nan]
#        # [ 0.  0. nan]
#        # [ 0.  0.  0.]]
