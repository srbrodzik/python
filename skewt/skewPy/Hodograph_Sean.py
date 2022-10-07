#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 14:55:37 2019

@author: broneil
"""

    def make_skewt_axes(self, pmax=1100., pmin=50.):

        self.fig = figure(figsize = (9,8))
        self.fig.clf()

        rcParams.update({'font.size': 10,})

        P = linspace(pmax, pmin, 37)
        pres_levels = np.arange(1000, 0, -100)

        self.skewxaxis = self.fig.add_axes([.085,.1,.73,.8], projection='skewx')
        self.skewxaxis.set_yscale('log')

        xticklocs = arange(-80,45,10)
        T0 = xticklocs

        w = array([0.0001,0.0004,0.001, 0.002, 0.004, 0.007, 0.01, 0.016, 0.024, 0.032])
        self.skewxaxis.add_mixratio_isopleths(w,P[P>=550],color='purple',ls='--',alpha=1.,lw=0.5)
        self.skewxaxis.add_dry_adiabats(linspace(250,440,20)-273.15,P,color='red',ls='--',alpha=1.,lw=0.5)
        self.skewxaxis.add_moist_adiabats(linspace(8,32,7),P[P>=200],color='g',ls='--',alpha=1.,lw=0.5)
        self.skewxaxis.other_housekeeping()

        #print self.skewxaxis.get_xlabel()

        self.wbax = self.fig.add_axes([0.75,0.1,0.1,0.8],sharey=self.skewxaxis,frameon=False) # wind barb axis
        self.wbax.xaxis.set_ticks([],[])
        self.wbax.yaxis.grid(True,ls='-',color='y',lw=0.5)
        for tick in self.wbax.yaxis.get_major_ticks():
        # tick.label1On = False
            pass
        self.wbax.get_yaxis().set_tick_params(size=0,color='y')
        self.wbax.set_xlim(-1.5,1.5)
        self.wbax.get_yaxis().set_visible(False)

        # Set up standard atmosphere height scale on
        # LHS of plot. It's jus
        majorLocatorKM   = MultipleLocator(2)
        majorLocatorKFT  = MultipleLocator(5)
        minorLocator     = MultipleLocator(1)

        self.htax=self.fig.add_axes([0.1,0.1,1e-6,0.8], sharey=self.skewxaxis, frameon=False)
        self.htax.xaxis.set_ticks([],[])
        self.htax.spines['left'].set_color('k')
        self.htax.spines['right'].set_visible(False)
        self.htax.get_yaxis().set_tick_params(size=0,color='y')
        self.htax.get_yaxis().set_visible(False)
        pres_heights = np.array( [self.data['hght'][np.argmin(np.abs(_ - self.data['pres']))] for _ in pres_levels] )
        for ph in range(pres_levels.shape[0]):
            self.htax.text(0, pres_levels[ph], '%.0f m'%pres_heights[ph], fontsize = 8)
        #print pres_heights

        # This is a hack for the xlabel not working!
        self.fig.text(0.40, 0.05, 'Temperature (%s)'%(degC), fontsize=12)

        self.fig.text(0.83, 0.9, 'Sounding Params')
        #$\underline{Sounding Params}$')
        #h0 = interp(0, self.data['temp'], self.data['hght'])
        h0 = self.data['hght'][np.argmin(abs(self.data['temp']-0))]
        h10 = self.data['hght'][np.argmin(abs(self.data['temp']+10))]
        h20 = self.data['hght'][np.argmin(abs(self.data['temp']+20))]
        h30 = self.data['hght'][np.argmin(abs(self.data['temp']+30))]
        h40 = self.data['hght'][np.argmin(abs(self.data['temp']+40))]
        lcl_height = 216*(self.data['temp'][0]-self.data['dwpt'][0])
        wcd = h0 - lcl_height
        mr = MixRatio(SatVap(self.data['dwpt'][0]), self.data['pres']*100)
        pdiff = -1.0*np.diff(self.data['pres'])
        #print mr.shape, pdiff.shape
        # precipitable water
        pw = np.sum(np.array([mr[_]*100.0*pdiff[_]/9.8 for _ in range(pdiff.shape[0]) ]))
        # crude estimate of wet bulb temp
        ##tw = self.data['temp'][0]- (self.data['temp'][0]-self.data['dwpt'][0])/3.
        # now some shear calculations
        wspd6km = self.data['sknt'][np.argmin(abs(self.data['hght']-6000))]
        wdir6km = self.data['drct'][np.argmin(abs(self.data['hght']-6000))]

        udiff = wspd6km*np.cos(np.radians(270-wdir6km)) - self.data['sknt'][3]*np.cos(np.radians(270-self.data['drct'][3]))
        vdiff = wspd6km*np.sin(np.radians(270-wdir6km)) - self.data['sknt'][3]*np.sin(np.radians(270-self.data['drct'][3]))
        # print udiff, vdiff
        shear6km = np.sqrt(udiff**2 + vdiff**2)
        if (math.isnan(shear6km)):
            shear6km = 0

        # 850mb-200mb Shear
        wspd850mb = self.data['sknt'][np.argmin(abs(self.data['pres']-850))]
        wdir850mb = self.data['drct'][np.argmin(abs(self.data['pres']-850))]
        wspd200mb = self.data['sknt'][np.argmin(abs(self.data['pres']-200))]
        wdir200mb = self.data['drct'][np.argmin(abs(self.data['pres']-200))]

        udiff = wspd200mb*np.cos(np.radians(270-wdir200mb)) - wspd850mb*np.cos(np.radians(270-wdir850mb))
        vdiff = wspd200mb*np.sin(np.radians(270-wdir200mb)) - wspd850mb*np.sin(np.radians(270-wdir850mb))
        # print udiff, vdiff
        shear850200mb = np.sqrt(udiff**2 + vdiff**2)
        if (math.isnan(shear850200mb)):
            shear850200mb = 0

        # SFC-700mb Shear
        wspd700mb = self.data['sknt'][np.argmin(abs(self.data['pres']-700))]
        wdir700mb = self.data['drct'][np.argmin(abs(self.data['pres']-700))]

        udiff = wspd700mb*np.cos(np.radians(270-wdir700mb)) - self.data['sknt'][3]*np.cos(np.radians(270-self.data['drct'][3]))
        vdiff = wspd700mb*np.sin(np.radians(270-wdir700mb)) - self.data['sknt'][3]*np.sin(np.radians(270-self.data['drct'][3]))
        # print udiff, vdiff
        shear700mb = np.sqrt(udiff**2 + vdiff**2)
        if (math.isnan(shear700mb)):
            shear700mb = 0


        self.fig.text(0.84, 0.88, '0%s: %.0f m'%(degC, h0))
        self.fig.text(0.84, 0.86, '-10%s: %.0f m'%(degC, h10)) ### Move
        self.fig.text(0.84, 0.84, '-20%s: %.0f m'%(degC, h20))
        self.fig.text(0.84, 0.82, '-30%s: %.0f m'%(degC, h30)) ### Move
        self.fig.text(0.84, 0.80, '-40%s: %.0f m'%(degC, h40))
        self.fig.text(0.84, 0.78, 'sfc LCL: %.0f m'%lcl_height)
        self.fig.text(0.84, 0.76, 'WCD: %d m'%wcd)
        self.fig.text(0.84, 0.74, 'PW: %.1f mm'%pw)
        self.fig.text(0.84, 0.72, '6 km shear: %d kts'%shear6km)
        self.fig.text(0.84, 0.68, '850-200mb shear: \n       %d kts'%shear850200mb) ### Move
        self.fig.text(0.84, 0.64, 'SFC-700mb shear: \n       %d kts'%shear700mb) ### Move
        # self.fig.text(0.84, 0.62, 'Trop: %d hPa'%(tropopause_pressure))


        self.fig.text(0.84, 0.20, 'Surface')
        self.fig.text(0.85, 0.18, 'P: %.1f hPa'%self.data['pres'][0])
        self.fig.text(0.85, 0.16, 'H: %.0f m'%(self.data['hght'][0]))

        self.fig.text(0.85, 0.14, 'T: %.1f %s'%(self.data['temp'][0], degC))
        self.fig.text(0.85, 0.12, 'T$_{D}$: %.1f %s'%(self.data['dwpt'][0], degC))
        ##self.fig.text(0.85, 0.10, 'T$_W$: %.1f %s'%(tw, degC))


        # now do the hodograph?
        self.hodoax = self.fig.add_axes([0.565, 0.68, 0.21, 0.21], frameon=True, polar=True)
        self.hodoax.xaxis.set_ticks([], [])
        self.hodoax.yaxis.set_ticks([], [])
        speed_mask = np.ma.masked_where(self.data['sknt'] > 999, self.data['sknt'])
        angle = 270 - self.data['drct']
        rad_angle = np.radians(angle)
        pres_mask = np.bitwise_and(self.data['pres'] > 200., self.data['drct'] < 361.) # only look at valid winds below 200 mb

        # UWYO is far more sparse so must be interpolated to register connected lines on the scatter plotted hodograph
        if self.fmt == 'UWYO':
            #Making copies of the datasets where idx is the original indexes of the bounds of interpolated data
            sknt_idx = self.data['sknt'][pres_mask]
            #The final interpolated array with the original bounds and 8 interpolated data points between them
            sknt_interp = self.data['sknt'][pres_mask]
             
            rad_angle_idx = rad_angle[pres_mask]
            rad_angle_interp = rad_angle[pres_mask]
             
            height_idx = self.data['hght'][pres_mask]
            height_interp = self.data['hght'][pres_mask]
            
            j = 1 #Insert index of the first array of interpolated data
            for i in range(len(sknt_idx)-1):
                #Interpolates data using linspace, then inserts it into the dataset
                sknt_interp = np.insert(sknt_interp , j , np.linspace(sknt_idx[i],sknt_idx[i+1],10)[1:-1] )
                rad_angle_interp = np.insert(rad_angle_interp , j , np.linspace(rad_angle_idx[i],rad_angle_idx[i+1],10)[1:-1] )
                height_interp = np.insert(height_interp , j , np.around(np.linspace(height_idx[i],height_idx[i+1],10)[1:-1] ))
                j += 9 #8 values are added so the index will jump by 9 indices
                
            sknt = sknt_interp
            rad_ag = rad_angle_interp
            hght = height_interp
#            self.hodoax.plot(rad_angle, self.data['sknt'], c = 'red', linewidth = 3)

        else:
            sknt = self.data['sknt'][pres_mask]
            hght = self.data['hght'][pres_mask]
            rad_ag = rad_angle[pres_mask]
            
        round_val = 1000.0
        max_height = 12000
        
        rounded_height = round_val*np.round(hght/round_val)
        hodo_sc = self.hodoax.scatter(rad_ag, sknt, c=rounded_height, \
        edgecolors='none', s=3, cmap=plt.cm.hsv, vmin=0, vmax=max_height)
        cb_ax = self.fig.add_axes([0.585, 0.665, 0.17, 0.005])
        hodo_cb = plt.colorbar(hodo_sc, cax=cb_ax, orientation='horizontal', drawedges=False)
        #hodo_cb = plt.colorbar(hodo_sc, cax=cb_ax, orientation='horizontal', fraction=0.06, pad=0.01, drawedges=False)
        hodo_cb.outline.set_visible(False)
        hodo_cb.ax.set_xticklabels((np.arange(0, max_height+round_val, round_val*2)/1000.0).astype(int))
        hodo_cb.ax.tick_params(labelsize=6, axis='x', which='both',length=0)

            # for label in hodo_cb.ax.get_xticks():
            #     label.set_visible(False)



        # legend_heights = np.arange(0, 12000, round_val)

        # legend_radius = 120
        # legend_angle_range = [30, 80]

        # legend_angles = np.linspace(legend_angle_range[0], legend_angle_range[1], len(legend_heights))
        # print 'legend angles: {}'.format(legend angles)

        # self.hodoax.scatter(legend_angles, np.zeros(len(legend_angles), float) + legend_radius, c=legend_heights,
        #                                 edgecolors='none', s=10, cmap=plt.cm.jet_r, vmin=0, vmax=12000)



        self.hodoax.set_yticks(np.arange(0, 120, 25))
        self.hodoax.tick_params(labelsize=5)
        self.hodoax.set_xticks(np.arange(0, 2*np.pi, np.pi/2))
        self.hodoax.set_xticklabels([])
        self.hodoax.grid(True)
        try:
            self.hodoax.set_rlabel_position(180)
        except AttributeError:
            pass

    def calc_pw(self):
        mr = MixRatio(SatVap(self.data['dwpt']), self.data['pres']*100)
        pdiff = -1.0*np.diff(self.data['pres'])
        #print mr.shape, pdiff.shape
        # precipitable water
        pw = np.sum(np.array([mr[_]*100.0*pdiff[_]/9.8 for _ in range(pdiff.shape[0]) ]))
	return pw
