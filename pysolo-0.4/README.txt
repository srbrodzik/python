A Prototype version of PySolo radar display using PySide, which is a set of LGPL-licensed Python bindings 
for Qt.

PySolo has been tested with Python 2.7.3, matplotlib 1.2.0, PySide 1.1.2.
Some earlier versions may work, however there are known problems with using PyQt
under RHEL 6.3 (Python 2.6.6)

Instructions for running PySolo with the Anaconda Python environment:

Anaconda is a free Python distribution from Continuum Analytics  that includes numpy, 
matplotlib and pyside, among many other scientific computing packages.

Anaconda supports the following platforms:

Linux and Windows: 64-bit and 32-bit
Mac OS X: Intel 64-bit
Python versions 2.6, 2.7, 3.3

1) Download Anaconda from https://store.continuum.io/cshop/anaconda and
install according to the provided instructions.

2) Verify that Anaconda Python is first in your path (logout/login, or source ~/.bashrc)
$ which python
/usr/local/anaconda/bin/python

For Mac OS users, the "graphical" version of Anaconda must be first in your
path:
export PATH=/usr/local/anaconda/python.app/Contents/MacOS:/usr/local/anaconda/bin:$PATH

3) use the "conda" package manager to install the necessary packages

$ conda install netcdf4

4) install the PySolo package:
python setup.py install

5) For first time users, run "start_pysolo" to create a default version of ~/pysolo_config.xml.

Edit the 'data_dir' parameter to set the default directory for data files

6) run 'start_pysolo' to start the PySolo application


Note:
PySolo also works with Enthought Canopy Basic Edition or above, which can be
downloaded from https://www.enthought.com/store/.  Enthough Canopy is free for
students and staff of degree granting institutions.  The Enthought Canopy
Express edition will not work - it doesn't contain several packages needed by
PySolo.
