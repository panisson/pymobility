#!/usr/bin/env python
# coding: utf-8
#from distutils.core import setup
from setuptools import setup

setup(name="pymobility",
      version="0.1.0",
      description="Mobility models for Python",
      license="GPL",
      author="André Panisson",
      author_email="panisson@gmail.com",
      url="http://github.com/panisson/pymobility",
      keywords= "mobility models",
      packages=["pymobility", "pymobility.models"],
      package_dir = {'': 'src'},
      include_package_data=True,
      install_requires=['numpy', 'matplotlib',],
      classifiers=[
                   'License :: OSI-Approved Open Source :: GNU General Public License version 2.0 (GPLv2)',
                   'Programming Language :: Python',
                   'Operating System :: OS Independent',
                   'Environment :: Console',
                   'User Interface :: Textual :: Command-line'
                   'Topic :: Scientific/Engineering :: Simulations',
                   'Topic :: Scientific/Engineering :: Information Analysis',
                   'Development Status :: 3 - Alpha',
                   'Intended Audience :: Science/Research',
                   ],
      )