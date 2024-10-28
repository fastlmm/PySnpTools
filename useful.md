# Useful

```bash
export QCTOOLPATH="/mnt/c/Users/carlk/OneDrive/Projects/Science/virginia/qctool/qctool"
set QCTOOLPATH=ubuntu run /mnt/c/Users/carlk/OneDrive/Projects/Science/virginia/qctool/qctool

uv sync --all-extras
pytest
```

```bash
set PROMP
set pythonpath=%cd%
cd tests
python tests.py
cd ..
pytest cmk/cmk/test_bgen2.py
```

Create a local distribution"

```bash
pip install twine

python setup.py sdist
twine upload dist/pysnptools-0.5.12b3.tar.gz

```
