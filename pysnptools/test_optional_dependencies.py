"""Exercise optional formats without hiding broken backend imports."""

import builtins
from importlib.util import find_spec
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from pysnptools.pstreader import PstData, PstHdf5


_has_hdf5 = find_spec("h5py") is not None
_requires_hdf5 = unittest.skipUnless(_has_hdf5, "requires the pysnptools[hdf5] extra")


class TestOptionalHdf5(unittest.TestCase):
    def test_core_imports_without_h5py(self):
        subprocess.run(
            [
                sys.executable,
                "-c",
                "import sys; sys.modules['h5py'] = None; sys.modules['bgen_reader'] = None; "
                "from pysnptools.snpreader import Bed, SnpData, SnpHdf5; "
                "from pysnptools.kernelreader import KernelHdf5; "
                "from pysnptools.distreader import DistHdf5; "
                "data = SnpData(iid=[['a','a'],['b','b']], sid=['s'], val=[[0.0],[2.0]]); "
                "assert data.standardize().val.shape == (2, 1)",
            ],
            check=True,
        )

    def test_missing_h5py_read(self):
        from pysnptools.snpreader import SnpHdf5
        from pysnptools.kernelreader import KernelHdf5
        from pysnptools.distreader import DistHdf5

        for reader in [PstHdf5, SnpHdf5, KernelHdf5, DistHdf5]:
            with self.subTest(reader=reader), patch.dict(sys.modules, {"h5py": None}):
                with self.assertRaisesRegex(ModuleNotFoundError, r"pysnptools\[hdf5\]") as caught:
                    reader("unused.hdf5").read()
                self.assertEqual(caught.exception.name, "h5py")
                self.assertIsInstance(caught.exception.__cause__, ModuleNotFoundError)

    def test_missing_h5py_write(self):
        data = PstData(row=["a"], col=["b"], val=[[1.0]])
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "unused.hdf5"
            with patch.dict(sys.modules, {"h5py": None}):
                with self.assertRaisesRegex(ModuleNotFoundError, r"pysnptools\[hdf5\]"):
                    PstHdf5.write(filename, data)
            self.assertFalse(filename.exists())

    def test_broken_backend_preserves_exception(self):
        original_import = builtins.__import__
        errors = [
            ModuleNotFoundError("missing backend dependency", name="backend_dependency"),
            ImportError("incompatible h5py binary"),
            OSError("missing HDF5 shared library"),
        ]
        for error in errors:
            def broken_import(name, *args, **kwargs):
                if name == "h5py":
                    raise error
                return original_import(name, *args, **kwargs)

            with self.subTest(error=error), patch("builtins.__import__", broken_import):
                with self.assertRaises(type(error)) as caught:
                    PstHdf5("unused.hdf5").read()
                self.assertIs(caught.exception, error)


class TestOptionalBgen(unittest.TestCase):
    def test_missing_backend_read(self):
        from pysnptools.distreader import Bgen

        with patch.dict(sys.modules, {"bgen_reader": None}):
            with self.assertRaisesRegex(ModuleNotFoundError, r"pysnptools\[bgen\]") as caught:
                Bgen("unused.bgen").read()
            self.assertEqual(caught.exception.name, "bgen_reader")
            self.assertIsInstance(caught.exception.__cause__, ModuleNotFoundError)

    def test_missing_backend_write(self):
        from pysnptools.distreader import Bgen

        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "unused.bgen"
            with patch.dict(sys.modules, {"bgen_reader": None}):
                with self.assertRaisesRegex(ModuleNotFoundError, r"pysnptools\[bgen\]"):
                    Bgen.write(filename, None)
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_broken_backend_preserves_exception(self):
        from pysnptools.distreader import Bgen

        original_import = builtins.__import__
        error = ModuleNotFoundError("missing native BGEN dependency", name="cbgen")

        def broken_import(name, *args, **kwargs):
            if name == "bgen_reader":
                raise error
            return original_import(name, *args, **kwargs)

        with patch("builtins.__import__", broken_import):
            with self.assertRaises(ModuleNotFoundError) as caught:
                Bgen("unused.bgen").read()
            self.assertIs(caught.exception, error)
