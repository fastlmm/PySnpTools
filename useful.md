# Useful

```bash
export QCTOOLPATH="/mnt/c/Users/carlk/OneDrive/Projects/Science/virginia/qctool/qctool"
set QCTOOLPATH=ubuntu run /mnt/c/Users/carlk/OneDrive/Projects/Science/virginia/qctool/qctool

uv sync --all-extras
python tests/test.py
pytest
```

```bash
cd doc
make.bat html
build\html\index.html
xcopy /c /e /s /h build\html ..\docs
cd ..
```

* Download and extract wheel artifacts from GitHub.

```bash
cd /d "C:\Users\carlk\Downloads\wheels (43)"
twine upload pysnptools*
```

Create a local distribution"

```bash
uv build
twine upload dist/pysnptools-0.????.tar.gz

```
