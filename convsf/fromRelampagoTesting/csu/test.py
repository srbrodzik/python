import numpy as np
import sys

missing = -99.

a = np.array([1.,2.,3.,4.])
print >>sys.stderr, 'a = ', a

a[a>3]=np.nan
print >>sys.stderr, 'a = ', a

a[a==np.nan]=missing
print >>sys.stderr, 'a = ', a

b = np.isnan(a)
print >>sys.stderr, 'b = ', b

a[np.isnan(a)]=missing
print >>sys.stderr, 'a = ', a



