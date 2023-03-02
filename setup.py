#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This file is used to create the package we'll publish to PyPI.

.. currentmodule:: setup.py
.. moduleauthor:: Benoit Valot <benoit.valot@univ-fcomte.fr>
"""

import importlib.util
import os
from pathlib import Path
from setuptools import setup, find_packages
from codecs import open  # Use a consistent encoding.
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Get the base version from the library.  (We'll find it in the `version.py`
# file in the src directory, but we'll bypass actually loading up the library.)
vspec = importlib.util.spec_from_file_location(
  "version",
  str(Path(__file__).resolve().parent /
      'pymlst'/"version.py")
)
vmod = importlib.util.module_from_spec(vspec)
vspec.loader.exec_module(vmod)
version = getattr(vmod, '__version__')

# If the environment has a build number set...
if os.getenv('buildnum') is not None:
    # ...append it to the version.
    version = "{version}.{buildnum}".format(
        version=version,
        buildnum=os.getenv('buildnum')
    )

setup(
    name='PyMLST',
    description="python Mlst Local Search Tool",
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(
        exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
    version=version,
    setup_requires=['wheel'],
    install_requires=[
        # Include dependencies here
        'biopython>=1.78',
        'click>=7.1',
        'pytest>=6.2',
        'sqlalchemy>=1.4,<2',
        'networkx>=2.5',
        'decorator>=4.4',
        'requests>=2.23',
        'pandas>=1.2',
        'numpy>=1.20.0',
        'beautifulsoup4>=4.9',
        'questionary>=1.9',
        'setuptools>=44.0',
        'alembic>=1.6'
    ],
    entry_points="""
    [console_scripts]
    pyMLST=pymlst.cmd:py
    wgMLST=pymlst.cmd:wg
    claMLST=pymlst.cmd:cla
    """,
    python_requires=">=3.7.0",
    license='GPLv3',  # noqa
    author='Benoit Valot',
    author_email='benoit.valot@univ-fcomte.fr',
    # Use the URL to the github repo.
    url='https://github.com/bvalot/pyMLST.git',
    download_url=(
        f'https://github.com/bvalot/pyMLST/archive/refs/tags/{version}.tar.gz'
    ),
    keywords=[
        'cgMLST', 'MLST', 'bacterial genome'
        # Add package keywords here.
    ],
    # See https://PyPI.python.org/PyPI?%3Aaction=list_classifiers
    classifiers=[
      # How mature is this project? Common values are
      #   3 - Alpha
      #   4 - Beta
      #   5 - Production/Stable
      'Development Status :: 5 - Production/Stable',

      # Indicate who your project is intended for.
      'Intended Audience :: Developers',
      'Topic :: Software Development :: Libraries',

      # Pick your license.  (It should match "license" above.)
        # noqa
      '''License :: OSI Approved :: GNU General Public License v3 (GPLv3)''',
        # noqa
      # Specify the Python versions you support here. In particular, ensure
      # that you indicate whether you support Python 2, Python 3 or both.
      'Programming Language :: Python :: 3.7',
    ],
    include_package_data=True
)
