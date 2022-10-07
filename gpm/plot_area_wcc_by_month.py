import numpy as np
import matplotlib.pyplot as plt
%matplotlib

my_xticks = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
x = np.array([0,1,2,3,4,5,6,7,8,9,10,11])
area_nasa = np.array([333001,287299,351934,523205,552701,599786,525567,542693,431230,424488,282238,258472])
area_uw2 = np.array([118793,78797,123853,294653,340941,387377,365618,407538,286517,252248,161272,87209])
title = "Total Area of WCC's by Month - NASA vs UW"

plt.xticks(x,my_xticks)
plt.plot(x,area_nasa)
plt.plot(x,area_uw2)
plt.legend(('NASA','UW'),loc='upper left')
plt.title(title)
plt.ylabel('km^2')
plt.tight_layout()
#plt.show()
