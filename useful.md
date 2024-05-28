# Useful

```bash
set PROMPT=$P$G
set pythonpath=%cd%
cd tests
python tests.py
cd ..
pytest cmk/cmk/test_bgen2.py
```

Create a local distribution"

```bash
python setup.py sdist
pip install twine
twine upload dist/pysnptools-0.5.12b1.tar.gz

```