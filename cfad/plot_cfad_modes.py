#matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
print("Dependencies loaded...")

indir_dow = '/home/disk/bob/olympex/zebra/moments/dow_lo_qc_to20km/rhi'
outdir = indir_dow
indir_npol = '/home/disk/bob/olympex/zebra/moments/npol_qc2_ctrOnDow_to20km/rhi'
#indir_gpm_npol = '/home/disk/bob/gpm/n1_oly_npol_ku/classify/ex_data'
indir_gpm_ocean = '/home/disk/bob/gpm/n1_oly_ocean_ku/classify/ex_data'

min_refl_for_plot = 10.0
refl_offset = -4.5

fig = plt.figure()
fig.set_size_inches(6,6)
ax = plt.axes()

xval_dow = np.loadtxt(indir_dow+'/dow_modes_'+str(int(min_refl_for_plot))+'db.txt')
xval_dow = xval_dow[0:13]
xval_npol_raw = np.loadtxt(indir_npol+'/npol_modes_dbz_tot_'+str(int(min_refl_for_plot))+'db_adj-0.0.txt')
xval_npol_raw = xval_npol_raw[0:13]
xval_npol = np.loadtxt(indir_npol+'/npol_modes_'+str(int(min_refl_for_plot))+'db.txt')
xval_npol = xval_npol[0:13]
xval_npol_adj = np.loadtxt(indir_npol+'/npol_modes_'+str(int(min_refl_for_plot))+'db_adj'+str(refl_offset)+'.txt')
xval_npol_adj = xval_npol_adj[0:13]
#xval_gpm_npol = np.loadtxt(indir_gpm_npol+'/gpm_modes_all.txt')
#xval_gpm_npol = xval_gpm_npol[0:13]
#xval_gpm_ocean = np.loadtxt(indir_gpm_ocean+'/gpm_modes_over_ocean_adj_all_levels.txt')
#xval_gpm_ocean = xval_gpm_ocean[0:49]
xval_gpm_ocean = np.loadtxt(indir_gpm_ocean+'/gpm_modes_over_ocean_adj_smoothed_'+str(int(min_refl_for_plot))+'db.txt')
xval_gpm_ocean = xval_gpm_ocean[0:13]

yval = np.linspace(2,8,13)
yval_gpm = np.linspace(2,8,49)

plt.plot( xval_dow, yval, color='blue', linestyle='solid', label='DOW' )
plt.plot( xval_npol_raw, yval, color='green', linestyle='solid', label='NPOL ZZ OverDow' )
plt.plot( xval_npol, yval, color='red', linestyle='solid', label='NPOL CZ OverDow' )
plt.plot( xval_npol_adj, yval, color='red', linestyle='dashed', label='NPOL CZ'+str(refl_offset)+'db OverDow' )
#plt.plot( xval_gpm_npol, yval, color='green', linestyle='dashed', label='GPM-OverNpol' )
#plt.plot( xval_gpm_ocean, yval_gpm, color='black', linestyle='solid', label='GPM-OverOcean' )
plt.plot( xval_gpm_ocean, yval, color='black', linestyle='solid', label='GPMKu+1.7db OverOcean' )

plt.xlim()
plt.ylim()

plt.title("Mode Values at each Height")
plt.xlabel("Reflectivity")
plt.ylabel("Height (km)")

plt.grid()

plt.legend()

#plt.show()
#plt.savefig(outdir+'/dow_npol_gpmAllLevels_modes_'+str(int(min_refl_for_plot))+'db.png', bbox_inches='tight')
#plt.savefig(outdir+'/dow_npol_gpmSmoothed_modes_'+str(int(min_refl_for_plot))+'db.png', bbox_inches='tight')
plt.savefig(outdir+'/dow_npolUncal_npolCal_gpmSmoothed_modes_'+str(int(min_refl_for_plot))+'db_adj'+str(refl_offset)+'.png', bbox_inches='tight')
#plt.savefig(outdir+'/dow_npol_modes_'+str(int(min_refl_for_plot))+'db_adj'+str(refl_offset)+'.png', bbox_inches='tight')
