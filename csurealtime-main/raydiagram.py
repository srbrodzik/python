import numpy as np 
import matplotlib.pyplot as plt 
plt.switch_backend('agg')

colors = ['black', 'purple', 'blue', 'green', 'gold', 'darkorange', 'red', '0.5']
tile_colors = np.tile(colors, 10)


#all_angles = np.array([0.5, 0.8, 1.5, 2.0, 3.0, 3.75, 4.75, 6, 8, 10, 15, 20, 25, 30.0], np.float)

#  rain scan
#angles = np.array([0.5, 0.8, 1.5, 2.0, 3.0, 3.75, 4.75], np.float)
# full scan
#angles = np.array([0.5, 0.8, 1.5, 2.0, 3.0, 3.75, 4.75, 6.0, 7.0, 8, 9, 10, 11, 12, 13, 14, 15, 16.5, 18.0, 20, 22, 26, 29], np.float)
# far #scan
#angles = np.array([0.5, 0.8, 1.5, 2.0, 3.0, 3.75, 4.75, 6.0, 7.0, 8, 9, 10, 11, 12, 13], np.float)

# piston scan
#angles = np.array([0.8, 1.21, 1.62, 2.03, 2.44, 2.85, 3.25, 3.66, 4.07, 4.48, 4.93, 5.42, 5.96, 6.55,
				#7.21, 7.93, 8.72, 9.58, 10.54, 11.58, 12.73, 14.00, 15.38, 16.90])

# test one?
# angles = np.array([0.50, 0.91, 1.32, 1.73, 2.14, 2.55, 2.95, 3.36, 3.77, 4.18, 4.60, 5.06, 5.56, 6.12, 6.73, 7.40, 8.14, 
# 						8.95, 9.84, 10.82, 11.89, 13.07, 14.37, 15.79, 17.35])

# PPI far
# angles = np.array([0.80, 1.64, 2.48, 3.33, 4.17, 5.01, 5.85, 6.69, 7.52, 8.36, 9.19, 10.02, 10.85, 11.68, 12.50, 13.33, 
# 	14.15, 15.00, 15.89, 16.85, 17.85, 18.92, 20.04, 21.23])


# # PPI near
# angles = np.array([0.80, 1.99, 3.19, 4.38, 5.57, 6.76, 7.94, 9.12, 10.30, 11.48, 12.65, 13.81, 14.97, 
# 					16.12, 17.33, 18.62, 20.00, 21.49, 23.08, 24.78, 26.60, 28.54, 30.61, 32.82])

# RHIs??
angles = np.array([0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 18, 21, 24, 27, 30, 33, 36, 39])


radius = np.arange(0, 125)
lw = 0.8

scan_rate = 17.0

DPI = 150

#sin_angles = np.sin(np.radians(angles))
fill_flag = False



fig, ax = plt.subplots()

for a, sa in enumerate(angles):
	h_beams = radius*np.sin(np.radians(sa))+radius**2/(2*6370.0)
	ax.plot(radius, h_beams, label=angles[a], color=tile_colors[a])
	#ax.plot(radius, radius*np.sin(np.radians(sa)), linestyle='dashed', color=colors[a], linewidth=0.5, zorder=10-a)

	top = h_beams+(1.0)*np.pi*radius/(180.0*2)
	bottom = h_beams-(1.0)*np.pi*radius/(180.0*2)

	if fill_flag:
		ax.fill_between(radius, bottom, top, facecolor=tile_colors[a], alpha=0.3, zorder=10-a)


scan_time = 360.0/scan_rate + 3

total_scan_time = scan_time*len(angles)

#ax.axhline(y=1.0, color='black', linewidth=lw)
#ax.axhline(y=2.0, color='black', linewidth=lw)
#ax.axvline(x=10, color='black', linewidth=lw)
#ax.axvline(x=100, color='black', linewidth=lw)
#ax.axhline(y=4.5, color='red', linestyle='dashed', linewidth=lw)

xlim = ax.get_xlim()
#print xlim

lc = 'black'
lw = 2.5

# ax.axhline(y=0.5, xmin=(10-xlim[0])/(xlim[1]-xlim[0]), xmax=(40-xlim[0])/(xlim[1]-xlim[0]), color=lc, linewidth=lw)
# ax.axhline(y=1.0, xmin=(40-xlim[0])/(xlim[1]-xlim[0]), xmax=(60-xlim[0])/(xlim[1]-xlim[0]), color=lc, linewidth=lw)
# ax.axhline(y=1.5, xmin=(60-xlim[0])/(xlim[1]-xlim[0]), xmax=(80-xlim[0])/(xlim[1]-xlim[0]), color=lc, linewidth=lw)
# ax.axhline(y=2.0, xmin=(80-xlim[0])/(xlim[1]-xlim[0]), xmax=(120-xlim[0])/(xlim[1]-xlim[0]), color=lc, linewidth=lw)


ax.grid(True)
ax.set_xlabel('Horizontal range (km)')
ax.set_ylabel('Height (km)')

ax.legend(loc='best', prop={'size': 4.5})

ax.set_ylim(0, 22)

ax.set_title('Height vs slant range: %d elevations, %d seconds'%(len(angles), total_scan_time))

plt.tight_layout()


plt.savefig('ppi_scan_rays_rhi.png', dpi=DPI)

