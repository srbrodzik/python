def raintype(self):
        print 'Setting up default C/S parameters. Can change in RadarData.py'
        minZdiff = 20 
        deepcoszero = 40
        shallowconvmin = 28
        #truncZconvthres = 43;
        #dBZformaxconvradius = 46;
        truncZconvthres = self.zconv
        dBZformaxconvradius = self.zconv+3;
        weakechothres = 7
        backgrndradius = 5       #(in km)
        maxConvRadius = 10       #(in km)
        minsize = 8              #(in km^2)
        startslope = 50          #(in km^2)
        maxsize = 2000           #(in km^2)

        if self.x_name == 'longitude':
            dx = self.dx*110.
        else:
            dx = self.dx
#        print dx
        zlev = self.get_ind(self.cs_z,self.data[self.z_name].data)
#        print zlev
        
        refl = np.squeeze(self.data[self.dz_name].sel(z=slice(zlev,zlev)).data)
        if len(np.shape(refl)) > 2:
            refl = np.squeeze(self.data[self.dz_name].sel(z=slice(zlev,zlev+1)).data)

        raintype,self.rtypes = rt.raintype(refl, refl_missing_val=self.data[self.dz_name].data.min(), 
                                   refl_dx=dx, minZdiff=minZdiff, deepcoszero=deepcoszero,
                                   shallowconvmin=shallowconvmin,truncZconvthres=truncZconvthres,
                                   dBZformaxconvradius=dBZformaxconvradius,
                                   weakechothres=weakechothres, backgrndradius=backgrndradius,
                                   maxConvRadius=maxConvRadius,minsize=minsize,
                                   startslope=startslope, maxsize=maxsize)
        nlevs = np.shape(self.data[self.z_name].data)[0]
#        print nlevs
        rpt = np.tile(raintype,(nlevs,1,1))
        self.raintype= rpt        
        self.def_convstrat()