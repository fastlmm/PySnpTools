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
pip install twine

python setup.py sdist
twine upload dist/pysnptools-0.5.12b3.tar.gz

```