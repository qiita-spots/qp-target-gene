#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2013, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from setuptools import setup
from glob import glob

__version__ = "0.1.0-dev"

classes = """
    Development Status :: 3 - Alpha
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering :: Bio-Informatics
    Topic :: Software Development :: Libraries :: Application Frameworks
    Topic :: Software Development :: Libraries :: Python Modules
    Programming Language :: Python
    Programming Language :: Python :: 2.7
    Programming Language :: Python :: Implementation :: CPython
    Operating System :: POSIX :: Linux
    Operating System :: MacOS :: MacOS X
"""

with open('README.rst') as f:
    long_description = f.read()

classifiers = [s.strip() for s in classes.split('\n') if s]

setup(name='qp-target-gene',
      version=__version__,
      long_description=long_description,
      license="BSD",
      description='Qiita Plugin: Target Gene',
      author="Qiita development team",
      author_email="qiita.help@gmail.com",
      url='https://github.com/qiita-spots/qp-target-gene',
      test_suite='nose.collector',
      packages=['qp_target_gene', 'qp_target_gene/split_libraries'],
      package_data={'qp_target_gene': ['support_files/config_file.cfg']},
      scripts=glob('scripts/*'),
      extras_require={'test': ["nose >= 0.10.1", "pep8"]},
      install_requires=['click >= 3.3', 'future', 'requests', 'pandas==0.20.3',
                        'biom-format <= 2.1.7',
                        'h5py >= 2.3.1', 'qiime >= 1.9.0, < 1.10.0',
                        'qiita-files @ https://github.com/'
                        'qiita-spots/qiita-files/archive/master.zip',
                        'qiita_client @ https://github.com/qiita-spots/'
                        'qiita_client/archive/master.zip'],
      classifiers=classifiers
      )
