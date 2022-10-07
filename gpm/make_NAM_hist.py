import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

%matplotlib

indir = '/home/disk/user_www/brodzik/public/nam_study/nam_data'

# convective length
dataPlains = [86,164,63,185,95,52,149,226,44,97,48,52,105,97,197,96,65,300,88,
              177,78,74,50,73,194,83,97,97,66,300,65,200,128,115,66,176,80,54,
              176,94,96,54,100,74,230,164,56,58,283,68,129,136,59,83,162,62,53,
              46,59,234,525,73,122,57,78,119,145,180,118,170,166,83,136,93,130,
              98,91,77,60,55,236,184,76,87,174,85,78,120,116,90,70,72,112,107,
              57,160,68,50,65]
dataGulf = [98,141,49,175,64,77,88,136,73,211,196,167,109,41,91,71,63,164,85,
            86,73,58,138,142,105,75,67,82,57,48,62,61,77,38,76,55,61,204,63]

plt.hist(dataPlains,alpha=0.3)
plt.hist(dataGulf,alpha=0.3)
plt.title('Convective Length')
plt.xlabel('Length (km)')
plt.ylabel('Counts')
plt.figure()


# stratiform major axis
dataPlains = [40,112,60,173,203,45,61,190,38,237,237,237,237,237,88,99,61,246,
              121,190,190,190,40,53,355,183,186,77,77,334,269,256,73,29,66,61,
              58,89,213,62,31,27,59,55,154,431,189,148,329,60,272,291,50,44,
              100,73,73,83,63,178,286,131,67,144,87,114,267,56,97,207,249,210,
              156,82,30,22,129,129,91,49,255,153,60,155,126,97,91,54,51,55,28,
              76,97,202,140,118,137,42,41]
dataGulf = [24,98,263,204,140,116,121,135,58,264,58,271,28,196,55,39,240,240,
            154,71,55,241,42,82,98,274,85,67,39,81,99,38,44,111,140,88,80,248,
            41]

plt.hist(dataPlains,alpha=0.3)
plt.hist(dataGulf,alpha=0.3)
plt.title('Stratiform Major Axis')
plt.xlabel('Length (km)')
plt.ylabel('Counts')
plt.figure()


# stratiform minor axis
dataPlains = []
dataGulf = [23,51,259,115,96,96,87,51,24,96,26,67,10,55,33,29,215,214,104,39,
            35,92,13,47,35,83,69,32,29,17,51,24,19,37,117,43,65,236,36]

plt.hist(dataPlains,alpha=0.3)
plt.hist(dataGulf,alpha=0.3)
plt.title('Stratiform Minor Axis')
plt.xlabel('Length (km)')
plt.ylabel('Counts')
plt.figure()
