from __future__ import absolute_import
from pysnptools.util.mapreduce1.runner import Runner,_JustCheckExists, _run_one_task
import os
import logging
try:
    import dill as pickle
except:
    logging.warning("Can't import dill, so won't be able to clusterize lambda expressions. If you try, you'll get this error 'Can't pickle <type 'function'>: attribute lookup __builtin__.function failed'")
    import pickle
import subprocess, sys, os.path
import multiprocessing
from pysnptools.util import create_directory_if_necessary, _datestamp


class LocalMultiProc(Runner):
    '''
    A :class:`.Runner` that runs a :func:`.map_reduce` as multiple processes on a single machine.

    **Constructor:**
        :Parameters: * **taskcount** (*number*) -- The number of processes to run on.
                     * **mkl_num_threads** (*number*) -- (default None) Limit on the number threads used by the NumPy MKL library.
                     * **weights** (*array of integers*) -- (default None) If given, tells the relative amount of work assigned to
                            each task. The length of the array must be **taskcount**. If not given, all tasks are assigned the
                            same amount of work.
                     * **taskindex_to_environ** (*function from integers to dictionaries*) -- (default None). If given, this
                            should be function from taskindex to a dictionary. The dictionary tells how to temporarily set
                            environment variables while the task is running. The dictionary is a mapping of
                            variables and values (both strings).
                     * **just_one_process** (*bool*) -- (default False) Divide the work for multiple processes, but runs sequentially on one process. Can be useful for debugging.
                     * **logging_handler** (*stream*) --  (default stdout) Where to output logging messages.
        
        :Example:

        >>> from pysnptools.util.mapreduce1 import map_reduce
        >>> from pysnptools.util.mapreduce1.runner import LocalMultiProc
        >>> def holder1(n,runner):
        ...     def mapper1(x):
        ...         return x*x
        ...     def reducer1(sequence):
        ...        return sum(sequence)
        ...     return map_reduce(range(n),mapper=mapper1,reducer=reducer1,runner=runner)
        >>> holder1(100,LocalMultiProc(4))
        328350

    '''

    def __init__(self, taskcount, mkl_num_threads = None, 
                 weights = None,
                 taskindex_to_environ = None,
                 just_one_process = False, logging_handler=logging.StreamHandler(sys.stdout)):
            # !!! cmk make so not just mkl_num_threads
        self.just_one_process = just_one_process
        logger = logging.getLogger()
        if not logger.handlers:
            logger.setLevel(logging.INFO)
        for h in list(logger.handlers):
            logger.removeHandler(h)
        if logger.level == logging.NOTSET:
            logger.setLevel(logging.INFO)
        logger.addHandler(logging_handler)

        self.taskcount = taskcount
        self.mkl_num_threads = mkl_num_threads
        self.weights = weights
        self.taskindex_to_environ = taskindex_to_environ

    def run(self, distributable):
        _JustCheckExists().input(distributable)

        localpath = os.environ["PATH"]
        localwd = os.getcwd()

        import datetime
        now = datetime.datetime.now()
        run_dir_rel = os.path.join("runs",_datestamp(appendrandom=True))
        run_dir_abs = os.path.join(localwd,run_dir_rel)
        create_directory_if_necessary(run_dir_rel, isfile=False)

        distributablep_filename = os.path.join(run_dir_rel, "distributable.p")
        with open(distributablep_filename, mode='wb') as f:
            pickle.dump(distributable, f, pickle.HIGHEST_PROTOCOL)

        distributable_py_file = os.path.join(os.path.dirname(__file__),"..","distributable.py")
        if not os.path.exists(distributable_py_file): raise Exception("Expect file at " + distributable_py_file + ", but it doesn't exist.")

        if not self.just_one_process:
            proc_list = []
            for taskindex in range(self.taskcount):
                environ = self.taskindex_to_environ(taskindex) if self.taskindex_to_environ is not None else None
                command_string_list = [sys.executable,
                                       distributable_py_file,
                                       distributablep_filename, 
                                       f"LocalInParts({taskindex},{self.taskcount},mkl_num_threads={self.mkl_num_threads},weights={self.weights},environ={environ})".replace(" ","")
                                       ]

                #logging.info(command_string_list)
                proc = subprocess.Popen(command_string_list, cwd=os.getcwd())
                proc_list.append(proc)

            for taskindex, proc in enumerate(proc_list):            
                rc = proc.wait()
                if not 0 == rc : raise Exception("Running python in python results in non-zero return code in task#{0}".format(taskindex))
        else:
            from pysnptools.util.mapreduce1.runner import LocalInParts
            for taskindex in range(self.taskcount):
                environ = self.taskindex_to_environ(taskindex) if self.taskindex_to_environ is not None else None
                LocalInParts(taskindex,self.taskcount, mkl_num_threads=self.mkl_num_threads, weights=self.weights, environ=environ).run(distributable)

        environ = self.taskindex_to_environ(self.taskcount) if self.taskindex_to_environ is not None else None
        result = _run_one_task(distributable, self.taskcount, self.taskcount, distributable.tempdirectory, weights=self.weights, environ=environ)

        _JustCheckExists().output(distributable)
        return result

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    from pysnptools.util.mapreduce1 import map_reduce #Needed to work around thread local variable issue

    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
