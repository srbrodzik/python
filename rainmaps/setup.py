"""
UW Rain Maps
"""

import os
import sys
#from setuptools import setup
from numpy.distutils.core import setup

# Pull the header into a variable
doclines = __doc__.split("\n")

VERSION = '1.0'

# Set variables for setup
PACKAGES = ['uw_rainmaps']

# Run setup
setup(
      name='uw_rainmaps',
      version=VERSION,
      #url='http://www.atmos.washington.edu/MG/PDFs/JTECH16_Powell-etal_RainCat.pdf',
      author='Stacy Brodzik',
      author_email='brodzik@uw.edu',
      description=("Rain-rate map code in Cartesian coordinates for DYNAMO Legacy Products"),
      license='GNU',
      packages=PACKAGES,
      classifiers=["""
          Development Status :: V1.0,
          Programming Language :: Python",
          Topic :: Scientific/Engineering
          Topic :: Scientific/Engineering :: Atmospheric Science
          Operating System :: Unix
          Operating System :: POSIX :: Linux
          Operating System :: MacOS
          """],
      long_description="""
          Python tool for rainmap creation.
          To access, use the following in your analysis code:
          from uw_rainmaps import rainmaps
          """,
      )
