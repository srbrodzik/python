import numpy as np
import numpy.ma as ma
from sharppy.sharptab import constants
from sharppy.sharptab.constants import MISSING
from sharppy.sharptab.profile import Profile
import numpy.testing as npt

sounding = """
 1000.00,    133.00,  -9999.00,  -9999.00,  -9999.00,  -9999.00
  976.00,    357.00,     22.20,     15.20,    160.00,     10.99
  971.00,    400.17,     21.20,     14.20,  -9999.00,  -9999.00
  946.89,    610.00,     18.79,     12.65,    165.00,     19.00
  943.00,    644.39,     18.40,     12.40,  -9999.00,  -9999.00
  925.00,    804.00,     16.80,     12.10,    165.00,     20.01
  913.06,    914.00,     15.75,     11.83,    165.00,     21.00
  880.76,   1219.00,     12.85,     11.07,    165.00,     22.01
  878.00,   1245.58,     12.60,     11.00,  -9999.00,  -9999.00
  850.00,   1517.00,     10.60,      9.50,    160.00,     21.00
  838.00,   1635.52,      9.80,      9.00,  -9999.00,  -9999.00
  818.70,   1829.00,      9.01,      7.88,    180.00,     14.01
  809.00,   1927.99,      8.60,      7.30,  -9999.00,  -9999.00
  789.06,   2134.00,      8.07,      4.56,    210.00,     19.00
  760.52,   2438.00,      7.29,      0.53,    240.00,     25.99
  750.00,   2553.00,      7.00,     -1.00,  -9999.00,  -9999.00
  732.74,   2743.00,      5.65,     -2.69,    255.00,     38.00
  700.00,   3116.00,      3.00,     -6.00,    260.00,     33.00
  671.00,   3456.95,      0.40,     -9.60,  -9999.00,  -9999.00
  653.95,   3658.00,     -1.42,    -10.46,    270.00,     25.99
  604.91,   4267.00,     -6.91,    -13.07,    290.00,     25.00
  581.75,   4572.00,     -9.67,    -14.38,    300.00,     25.00
  559.48,   4877.00,    -12.42,    -15.69,    300.00,     25.00
  551.00,   4996.30,    -13.50,    -16.20,  -9999.00,  -9999.00
  540.00,   5149.59,    -13.90,    -15.90,  -9999.00,  -9999.00
  537.69,   5182.00,    -14.07,    -16.43,    295.00,     24.01
  521.00,   5420.88,    -15.30,    -20.30,  -9999.00,  -9999.00
  500.00,   5730.00,    -17.90,    -22.90,    280.00,     29.99
  475.98,   6096.00,    -20.92,    -25.38,    280.00,     35.00
  473.00,   6142.70,    -21.30,    -25.70,  -9999.00,  -9999.00
  453.00,   6460.42,    -23.30,    -29.30,  -9999.00,  -9999.00
  400.00,   7360.00,    -29.70,    -35.70,    270.00,     42.00
  387.00,   7594.50,    -31.10,    -38.10,  -9999.00,  -9999.00
  385.60,   7620.00,    -31.31,    -38.43,    265.00,     46.00
  364.00,   8024.86,    -34.70,    -43.70,  -9999.00,  -9999.00
  337.75,   8534.00,    -39.19,    -48.19,    270.00,     56.99
  308.77,   9144.00,    -44.57,    -53.57,    280.00,     58.99
  300.00,   9340.00,    -46.30,    -55.30,    280.00,     54.99
  250.00,  10520.00,    -57.30,    -65.30,    285.00,     58.00
  249.00,  10545.31,    -57.50,    -65.50,  -9999.00,  -9999.00
  240.00,  10777.00,    -58.70,    -66.70,  -9999.00,  -9999.00
  214.00,  11487.91,    -63.90,    -70.90,  -9999.00,  -9999.00
  202.00,  11839.62,    -65.90,    -72.90,    280.00,     58.00
  200.43,  11887.00,    -65.74,    -72.74,    280.00,     58.00
  200.00,  11900.00,    -65.70,    -72.70,    275.00,     58.99
  187.00,  12308.99,    -65.90,    -72.90,  -9999.00,  -9999.00
  182.00,  12474.56,    -64.10,    -71.10,  -9999.00,  -9999.00
  179.00,  12577.22,    -61.30,    -69.30,  -9999.00,  -9999.00
  175.00,  12717.60,    -61.43,    -69.43,    280.00,     75.00
  172.64,  12802.00,    -61.51,    -69.51,    280.00,     75.00
  167.00,  13008.26,    -61.70,    -69.70,  -9999.00,  -9999.00
  164.40,  13106.00,    -60.46,    -68.46,    280.00,     56.99
  162.00,  13197.97,    -59.30,    -67.30,  -9999.00,  -9999.00
  158.00,  13354.96,    -59.10,    -67.10,  -9999.00,  -9999.00
  150.00,  13680.00,    -60.90,    -69.90,    280.00,     58.99
  145.00,  13890.13,    -61.70,    -70.70,  -9999.00,  -9999.00
  143.00,  13976.21,    -60.90,    -69.90,  -9999.00,  -9999.00
  131.00,  14521.26,    -60.30,    -69.30,  -9999.00,  -9999.00
  124.00,  14864.37,    -58.90,    -68.90,  -9999.00,  -9999.00
  118.00,  15174.39,    -60.10,    -70.10,  -9999.00,  -9999.00
  116.77,  15240.00,    -59.61,    -70.02,    275.00,     56.00
  115.00,  15335.36,    -58.90,    -69.90,  -9999.00,  -9999.00
  104.00,  15965.30,    -59.10,    -70.10,  -9999.00,  -9999.00
  100.00,  16210.00,    -60.70,    -71.70,    275.00,     49.01
   91.46,  16764.00,    -61.39,    -73.24,    285.00,     58.00
   90.10,  16856.89,    -61.50,    -73.50,  -9999.00,  -9999.00
   83.40,  17332.53,    -64.10,    -76.10,  -9999.00,  -9999.00
   74.50,  18026.50,    -61.90,    -74.90,  -9999.00,  -9999.00
   71.40,  18288.00,    -63.13,    -76.13,    285.00,     36.00
   70.00,  18410.00,    -63.70,    -76.70,    285.00,     36.00
   67.20,  18660.40,    -64.30,    -77.30,  -9999.00,  -9999.00
   61.55,  19202.00,    -62.05,    -75.92,    295.00,     39.01
   60.70,  19287.38,    -61.70,    -75.70,  -9999.00,  -9999.00
   58.57,  19507.00,    -62.35,    -76.35,    275.00,     38.00
   50.48,  20422.00,    -65.04,    -79.04,    315.00,     31.00
   50.30,  20443.48,    -65.10,    -79.10,  -9999.00,  -9999.00
   50.00,  20480.00,    -64.90,    -78.90,    305.00,     25.99
   48.20,  20703.50,    -64.10,    -78.10,  -9999.00,  -9999.00
   48.02,  20726.00,    -63.79,    -77.96,    295.00,     29.00
   46.20,  20964.56,    -60.50,    -76.50,  -9999.00,  -9999.00
   45.71,  21031.00,    -60.53,    -76.71,    305.00,     25.99
   43.52,  21336.00,    -60.69,    -77.66,    245.00,     25.99
   40.80,  21736.05,    -60.90,    -78.90,  -9999.00,  -9999.00
   38.60,  22083.01,    -56.90,    -78.90,  -9999.00,  -9999.00
   36.60,  22419.46,    -56.50,    -81.50,  -9999.00,  -9999.00
   35.82,  22555.00,    -57.40,    -83.31,    280.00,     19.00
   34.90,  22719.11,    -58.50,    -85.50,  -9999.00,  -9999.00
   34.13,  22860.00,    -56.96,    -87.45,    270.00,     27.00
   33.80,  22920.97,    -56.30,    -88.30,  -9999.00,  -9999.00
   33.10,  23053.37,    -56.90,    -88.90,  -9999.00,  -9999.00
   32.50,  23169.43,    -55.10,    -87.10,  -9999.00,  -9999.00
   31.00,  23470.00,    -54.83,    -86.83,    285.00,     41.01
   30.30,  23616.35,    -54.70,    -86.70,  -9999.00,  -9999.00
   30.00,  23680.00,    -53.70,    -86.70,    290.00,     31.99
   29.57,  23774.00,    -52.46,    -85.74,    295.00,     33.99
   28.50,  24013.41,    -49.30,    -83.30,  -9999.00,  -9999.00
   28.21,  24079.00,    -49.44,    -83.39,    300.00,     27.00
   26.93,  24384.00,    -50.09,    -83.81,    290.00,     22.01
   25.70,  24689.00,    -50.74,    -84.23,    255.00,     20.01
   24.52,  24994.00,    -51.39,    -84.64,    270.00,     15.00
   23.41,  25298.00,    -52.04,    -85.06,    245.00,     21.00
   23.30,  25327.59,    -52.10,    -85.10,  -9999.00,  -9999.00
   21.50,  25856.10,    -45.90,    -80.90,  -9999.00,  -9999.00
   21.33,  25908.00,    -45.66,    -80.81,    270.00,     31.99
   20.40,  26207.33,    -44.30,    -80.30,  -9999.00,  -9999.00
   20.38,  26213.00,    -44.33,    -80.29,    250.00,     27.00
   20.00,  26340.00,    -45.10,    -80.10,    250.00,     25.99
   18.59,  26822.00,    -46.96,    -81.50,    275.00,     31.00
   17.75,  27127.00,    -48.14,    -82.38,    250.00,     31.00
   17.10,  27375.31,    -49.10,    -83.10,  -9999.00,  -9999.00
   16.95,  27432.00,    -48.30,    -82.66,    255.00,     31.99
   16.70,  27531.12,    -46.90,    -81.90,  -9999.00,  -9999.00
   16.19,  27737.00,    -47.07,    -82.07,    260.00,     29.99
   16.10,  27773.04,    -47.10,    -82.10,  -9999.00,  -9999.00
   15.46,  28042.00,    -45.04,    -80.83,    220.00,     20.01
   15.30,  28111.81,    -44.50,    -80.50,  -9999.00,  -9999.00
   14.70,  28378.90,    -45.10,    -80.10,  -9999.00,  -9999.00
   14.12,  28651.00,    -41.29,    -77.95,    265.00,     24.01
   14.00,  28707.49,    -40.50,    -77.50,  -9999.00,  -9999.00
   13.50,  28955.39,    -39.50,    -76.50,  -9999.00,  -9999.00
   13.50,  28956.00,    -39.50,    -76.50,    275.00,     22.01
   13.00,  29212.33,    -41.10,    -78.10,  -9999.00,  -9999.00
   12.34,  29566.00,    -39.23,    -76.81,    260.00,     25.00
   11.90,  29816.29,    -37.90,    -75.90,  -9999.00,  -9999.00
   11.29,  30175.00,    -38.56,    -75.90,    285.00,     25.00
   11.00,  30356.07,    -38.90,    -75.90,  -9999.00,  -9999.00
   10.80,  30480.00,    -38.22,    -75.71,    265.00,     16.01
   10.60,  30610.53,    -37.50,    -75.50,  -9999.00,  -9999.00
   10.20,  30874.54,    -39.30,    -76.30,  -9999.00,  -9999.00
   10.00,  31010.00,    -39.10,    -76.10,    245.00,     20.01
    9.04,  31699.00,    -38.47,    -76.26,    255.00,     25.99
    8.80,  31887.35,    -38.30,    -76.30,  -9999.00,  -9999.00
"""
slines = sounding.split('\n')[1:-1]
pres = []
hght = []
tmpc = []
dwpc = []
wdir = []
wspd = []
for sline in slines:
    parts = sline.split(',')
    pres.append(float(parts[0]))
    hght.append(float(parts[1].strip()))
    tmpc.append(float(parts[2].strip()))
    dwpc.append(float(parts[3].strip()))
    wdir.append(float(parts[4].strip()))
    wspd.append(float(parts[5].strip()))
pres = ma.asarray(pres)
hght = ma.asarray(hght)
tmpc = ma.asarray(tmpc)
dwpc = ma.asarray(dwpc)
wdir = ma.asarray(wdir)
wspd = ma.asarray(wspd)


class TestProfile(object):
    def __init__(self):
        self.prof = Profile(pres=pres, hght=hght, tmpc=tmpc,
                            dwpc=dwpc, wdir=wdir, wspd=wspd)

    def test_prof_pres(self):
        pres[pres == MISSING] = ma.masked
        npt.assert_almost_equal(self.prof.pres, pres)

    def test_prof_hght(self):
        hght[hght == MISSING] = ma.masked
        npt.assert_almost_equal(self.prof.hght, hght)

    def test_prof_tmpc(self):
        tmpc[tmpc == MISSING] = ma.masked
        npt.assert_almost_equal(self.prof.tmpc, tmpc)

    def test_prof_dwpc(self):
        dwpc[dwpc == MISSING] = ma.masked
        npt.assert_almost_equal(self.prof.dwpc, dwpc)

    def test_prof_wdir(self):
        wdir[wdir == MISSING] = ma.masked
        npt.assert_almost_equal(self.prof.wdir, wdir)

    def test_prof_wspd(self):
        wspd[wspd == MISSING] = ma.masked
        npt.assert_almost_equal(self.prof.wspd, wspd)

    def test_find_sfc(self):
        correct_sfc_ind = 1
        returned_sfc_ind = self.prof.sfc
        npt.assert_equal(returned_sfc_ind, correct_sfc_ind)

    def test_find_sfc_missing_data(self):
        p = [1000, 925, 850, 700, 500, 400, 300, 250, 200]
        z = [MISSING, 769, 1487, 3099, 5750, 7400, 9430, 10660, 12090]
        t = [MISSING, 15.6, 13.6, 4.8, -14.7, -25.1, -38.1, -49.1, -59.1]
        d = [MISSING, 15.5, 11.5, 2.2, -17.8, -30.1, -46.1, -58.1, -67.1]
        wd = [MISSING, 155, 175, 225, 235, 240, 235, 240, 235]
        ws = [MISSING, 27, 24.01, 31.99, 44, 62.01, 85.01, 86, 94]
        prof = Profile(pres=p, hght=z, tmpc=t, dwpc=d, wdir=wd, wspd=ws)
        sfc_ind = 1
        npt.assert_almost_equal(prof.sfc, sfc_ind)


