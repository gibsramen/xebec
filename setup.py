# Much of this page is taken from the gemelli setup file
import ast
import re
from setuptools import find_packages, setup

# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r"__version__\s+=\s+(.*)")

with open("xebec/__init__.py", "rb") as f:
    hit = _version_re.search(f.read().decode("utf-8")).group(1)
    version = str(ast.literal_eval(hit))

with open("README.md") as f:
    long_description = f.read()

classes = """
    Development Status :: 3 - Alpha
    License :: OSI Approved :: BSD License
    Topic :: Software Development :: Libraries
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3 :: Only
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split("\n") if s]

description = (
    "Snakemake pipeline for microbiome diversity effect size "
    "benchmarking."
)

standalone = ["xebec=xebec.cli.cli:xebec"]

setup(
    name="xebec",
    author="Gibraan Rahman",
    author_email="grahman@eng.ucsd.edu",
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gibsramen/xebec",
    version=version,
    license="BSD-3-Clause",
    packages=["xebec/src", "xebec/cli", "xebec"],
    install_requires=[
        "numpy",
        "h5py==3.1.0",
        "evident>=0.2.0",
        "gemelli>=0.0.8",
        "pandas>=1.0.0",
        "scikit-bio>=0.5.6",
        "unifrac",
        "snakemake",
        "cookiecutter",
        "seaborn",
        "bokeh",
        "click"
    ],
    include_package_data=True,
    package_data={"": ["xebec/js", "xebec/{{cookiecutter.project_name}}"]},
    entry_points={"console_scripts": standalone},
    classifiers=classifiers,
)
