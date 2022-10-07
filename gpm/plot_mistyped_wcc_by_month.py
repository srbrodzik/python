import numpy as np
import matplotlib.pyplot as plt
%matplotlib

my_xticks = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
x = np.array([0,1,2,3,4,5,6,7,8,9,10,11])
percent = np.array([58.55, 66.10, 57.66, 42.98, 38.40, 34.82, 29.10, 23.64, 30.57, 38.54, 40.28, 61.48])
title = "Percentage of Mis-typed WCC's by Month"

plt.xticks(x,my_xticks)
plt.plot(x,percent)
plt.title(title)
plt.ylabel('% Mistyped WCC')
#plt.show()
