from setuptools import setup

# Version number
version = "0.5.13"


def readme():
    with open("README.md") as f:
        return f.read()


install_requires = [
    "scipy>=1.1.0",
    "numpy>=1.17.0",
    "pandas>=1.1.1",
    "psutil>=5.6.3",
    "h5py>=2.10.0",
    "cloudpickle>=2.2.0",
    "more-itertools>=8.7.0",
    "bed-reader[samples]>=1.0.5",
]

# Optional dependencies
extras_require = {
    "bgen": [
        "cbgen>=1.0.6",
        "bgen-reader>=4.0.8",
    ],
}

setup(
    name="pysnptools",
    version=version,
    description="PySnpTools",
    long_description=readme(),
    long_description_content_type="text/markdown",
    keywords="gwas bioinformatics sets intervals ranges regions plink genomics "
    + "file-format reader genotype bed-format writer python snps",
    url="https://fastlmm.github.io/",
    author="FaST-LMM Team",
    author_email="fastlmm-dev@python.org",
    project_urls={
        "Bug Tracker": "https://github.com/fastlmm/PySnpTools/issues",
        "Documentation": "http://fastlmm.github.io/PySnpTools",
        "Source Code": "https://github.com/fastlmm/PySnpTools",
    },
    license="Apache 2.0",
    classifiers=[
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python",
    ],
    packages=[  # basically everything with a __init__.py
        "pysnptools",
        "pysnptools/kernelreader",
        "pysnptools/kernelstandardizer",
        "pysnptools/pstreader",
        "pysnptools/snpreader",
        "pysnptools/distreader",
        "pysnptools/standardizer",
        "pysnptools/util",
        "pysnptools/util/filecache",
        "pysnptools/util/mapreduce1",
        "pysnptools/util/mapreduce1/runner",
    ],
    package_data={
        "pysnptools": [
            "util/pysnptools.hashdown.json",
            "util/bgen.hashdown.json",
            "tests/mintest.py",
        ]
    },
    install_requires=install_requires,
    extras_require=extras_require,
)
