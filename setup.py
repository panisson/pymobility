#!/usr/bin/env python
# coding: utf-8
from distutils.core import setup

setup(name="pymobility",
      version="0.1",
      description="Mobility models for python",
      license="GPL",
      author="André Panisson",
      author_email="panisson@gmail.com",
      url="http://github.com/panisson/pymobility",
      keywords= "mobility models",
      packages=["pymobility", "pymobility.models"],
      package_dir = {'': 'src'},
      include_package_data=True,
      install_requires=['numpy', 'matplotlib',],
      )