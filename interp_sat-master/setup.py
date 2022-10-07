#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Copyright © 2014 Martin Ueding <dev@martin-ueding.de>
# Licensed under The Lesser GNU Public License Version 2 (or later)

#from setuptools import setup
from numpy.distutils.core import setup

setup(
    author="Stacy Brodzik",
    description="Interpolation of Matlab ascii TRMM and GPM files",
    #license="LGPL2",
    name="interp_sat",
    packages=['interp_sat'],
    install_requires=[
        'numpy',
    ],
    #url="https://github.com/djpine/linfit",
    #download_url="https://github.com/djpine/linfit",
    #version="2014.09.03",
)
