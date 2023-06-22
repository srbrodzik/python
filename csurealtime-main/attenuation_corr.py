# Brody Fuchs, Oct 2018
# brfuchs@atmos.colostate.edu

# Rewriting the Carey et al 2000 attenuation correction method
# in python, because I can't find anything else that's very good
# at doing Zh and ZDR attenuation correction

import numpy as np 
import os
import glob
from copy import deepcopy
import scipy



def carey_atten_corr(dz, dp, dr=None, kd=None, rh=None, height=None, a_default=-0.08, deva=0.03, 
				bad=-999.9, force_default=False, phase_diff_calc=True,
				start_ind=None, end_ind=None, hfrz=10.0, min_pts=200):

	# going to set this up so it should be able to handle pyart data field objects
	# set the defaults up here


	nrays = dz.shape[0]
	#print 'nrays: {}'.format(nrays)

	zhcorr = deepcopy(dz)
	ah = np.zeros_like(dz) + bad
	a_array = np.zeros((nrays), float) # this will return the a values used for each ray


	for nr in range(nrays):

		# want this to be checked and calculated each time
		a_coeff = -1.0*a_default

		# this assumes the data is laid out in azimuth/elevation, gates
		dz_loop = dz[nr]
		dp_loop = dp[nr]
		kd_loop = kd[nr]

		if phase_diff_calc:
			# need to find first non masked value?
			phase_comp = dp_loop.compressed()
			phase_comp_valid = phase_comp[ (phase_comp > -360.0) & (phase_comp < 361) ]
			first_phase = np.average(phase_comp_valid[0:4])
			#print first_phase
			dp_loop = dp_loop - first_phase
		#	print dp_loop.compressed()

		#print 'dp loop: {}'.format(dp_loop)

		#print 'dz loop shape: {}'.format(dz_loop.shape)

		if not force_default:

#			print 'height: {}'.format(height[nr])
#			print 'dz: {}'.format(dz_loop)
#			print 'dp: {}'.format(dp_loop)
#			print 'kd: {}'.format(kd_loop)

#			print 'going to calculate a and b for each ray'
			# first find where the data is good/bad
			wh_gd = np.where((dz_loop != bad) & (dp_loop != bad) & (kd_loop >= 0.0) & (kd_loop <= 6.0) & 
							(height[nr] <= hfrz) & (height[nr] >= 0.5) & (dz_loop <= 60.0))

			# defining the corrected Z and the specific attenuation

			if len(wh_gd[0]) > 0:

				# here we're grabbing the good dz values and differential phase values
				zhdat = dz_loop[wh_gd]
				dpdat = dp_loop[wh_gd]
				# I'm guessing this is a flag to keep going with respect to something
				continue_run = True
				# I guess there needs to be a decent number of good 
				if len(wh_gd[0]) > min_pts:
					# ******* need to figure out what is returned from the "regress" function 
					# fit is the slope, everything else is returned from the keywords in the function call
					#fit = regress(dpdat,zhdat,chisq=chisq,const=const,correlation=corr,sigma=sigma,yfit=yfit)

					fit = scipy.stats.linregress(dpdat, zhdat)
					#print 'fit: {}'.format(fit)
					zhslope = fit[0]
					intercept = fit[1]
					corr = fit[2]
					sigma = fit[4]
					yfit = zhslope*dpdat + intercept

		#			print,zhslope
					numpts = len(dpdat)

					corr_check = np.abs(corr) < 0.90
					numpts_check = numpts > min_pts
	
	

					#print 'ray: {}, corr check: {}, num pts check: {}'.format(nr, corr_check, numpts_check)

					# corr gets returned from the regress function, but IDL's syntax is different and weird
					while corr_check and numpts_check and continue_run:
						#print 'checking again'
						# this is the difference between the actual values and the fit
						fit_diff = np.abs(zhdat-yfit)
						# Need to figure out what exactly the "moment" IDL function does and returns
						# I think sdev gets returned by the function

						# diffmom is just a dummy returned variable, what you really want is the standard deviation
						# that's easy to replace in python
						#diffmom = moment(diff, sdev=sdev)
						sdev = np.std(fit_diff)

						# checking where the difference is less than 2 standard deviations
						# probably to just filter out the outliers
						wh_good = np.where(np.abs(fit_diff) < 2*sdev)    

						# If there are some values, then go and do stuff
						# ******** I think this will basically run continuosly (and filter out outliers) until the fit is good enough
						if len(wh_good[0]):

						# select the Z and phase that satisfy the 2 sigma criteria
							dp_new = dpdat[wh_good]
							zh_new = zhdat[wh_good]
						# go and do another regression
						#regress2 = regress(dp_new, zh_new,chisq=chisq,const=const,correlation=corr,sigma=sigma,yfit=yfit)
							fit2 = scipy.stats.linregress(dp_new, zh_new)
							#print 'fit2: {}'.format(fit2)

						# overwrite the dpdat (phase data) and the zhdat (dbz data)
							dpdat = dp_new
							zhdat = zh_new
							numpts = len(dpdat)
							zhslope = fit2[0]
							intercept = fit2[1]
							corr = fit2[2]
							yfit = zhslope*dpdat + intercept

							corr_check = np.abs(corr) < 0.90
							numpts_check = numpts > min_pts

							#print 'ray: {}, corr2 check: {}, num pts2 check: {}'.format(nr, corr_check, numpts_check)

						else:
						    continue_run = False

					# Here we're done with the a calculation

					unc = sigma/(np.sqrt(numpts))
					#print 'ray: {}, unc: {}, max dp: {}, slope: {}'.format(nr, unc, np.max(dpdat), zhslope)
					if ( (corr**2 >= 0.25) and (numpts >= min_pts) and (unc <= 5.5) and (np.max(dpdat) >= 15) and 
									(zhslope >= (a_default-deva)) and (zhslope <= (a_default+deva)) ):
					# if we do get a good fit, then set the a coefficient as the slope
						a_coeff = -1*zhslope
					else:
						# if not, then use a default value
						#print 'Slope is not representative of the overall slope! ', zhslope, 'Using default a value ', a_default
						a_coeff = -1*a_default	


		else:
			#print 'just doing the default a and b'
			pass



		a_array[nr] = a_coeff

		# regardless of how the a_coeff was attained, do the filtering and actual attenuation
		# correction down here

		# wh_gdh looks to be where there are some valid points
		# dz and dp are defined in the function call, so no need to do anything additional

		wh_gdh = np.where( (dz_loop != bad) & (dp_loop > 0.0))[0]

		#print dp_loop[wh_gdh]
		#print a_coeff*(dp_loop[wh_gdh])

		if len(wh_gdh): # add the original dbz to the a coeff times the differential phase
			zhcorr[nr, wh_gdh] += a_coeff*(dp_loop[wh_gdh])
		else:
			#print 'no good dbz values: {}'.format(nr)
			pass

		#print zhcorr[nr, wh_gdh], dz_loop

		# only do this part if kdp is defined
		if kd is not None:

			# this is looking for Kdp here, wherever it's not bad
			wh_gkdh = np.where( (kd_loop != bad) & (kd_loop > 0.0))[0]

			#print 'doing specific attenuation'
			# specific attenuation just the Kdp multiplied by the a coefficient
			#print wh_gkdh
			if len(wh_gkdh):
				ah[nr, wh_gkdh] = kd_loop[wh_gkdh]*a_coeff



	return zhcorr, ah, a_array




def big_drop_correction(dz, dr, kd, dp, rh, height, a_star, b_star, a_coeff, b_coeff, bad):



    pass


# FUNCTION big_drop,dz,dr,kd,dp,rh,height,a_star,b_star,a_coeff,b_coeff,bad
# ############################################################################
# #Apply a big-drop correction as described by Carey et al. 2000.
# # Here, a* and b* are derived from scattering simulations (as given in Carey et al. 2000)
# #We also have to limit the data where delta is large and rhohv is small but there is some significant amount of Kdp.
# #This technique requires identifiying the "big drop core" by looking for points wthat meet certain requirements, and then
# #applying a correct to all range gates down range of the core.
# ############################################################################

# print,'Doing Big drop correction!'
# bigdrop_flag=dz*0
# zhcorr=dz
# drcorr=dr

# #print,'Max of dr: ', max(drcorr)

# nsweeps=n_elements(dz[0,0,*])
# nrays=n_elements(dz[0,*,0])
# ngates=n_elements(dz[*,0,0])

# tnumpts=0
# 	FOR pp=0,nsweeps-1 DO BEGIN
# 		ray_correction=0
# 		FOR qq=0,nrays-1 DO BEGIN
# #			FOR kk=0,ngates-1 DO BEGIN
# 			zhdat=dz[*,qq,pp]
# 			rhdat=rh[*,qq,pp]
# 			kddat=kd[*,qq,pp]
# 			drdat=dr[*,qq,pp]
# 			dpnewdat=dp[*,qq,pp]
# 			hgt=height(*,qq,pp)
# 			wh_bigdrops=where(rhdat < 0.80 and rhdat != bad and kddat > 1.0 and dpnewdat != bad and hgt < 5.0 and zhdat != bad,complement=wh_notbigdrops)
	
# 			IF wh_bigdrops[0] != -1 THEN BEGIN
# #				print,'Found big drops!'
# #				print,radar.volume[0].sweep[pp].ray[qq].h.azimuth
# 				nbd=1	
# 				bd_gates=fltarr(1)
# 			FOR chk=1,n_elements(wh_bigdrops)-1 DO BEGIN
# 				IF wh_bigdrops(chk) == wh_bigdrops(chk-1)+1 THEN BEGIN
# 					nbd=nbd+1
# 					bd_gates=[bd_gates,wh_bigdrops(chk)]
# 				ENDIF ELSE BEGIN
# 					IF nbd GE 4 THEN BEGIN
# #						print,'big drop region found with ', nbd, ' Correcting'
# 						tnumpts=tnumpts+nbd
# 					#Now use the a_star and b_star corrections from Carey et al. 2000
# 					#Zhcorr(r)=Zh(r)+aphidp(r)+(a*-a)[phi(r2)-phi(r1)]
# 					#Need to figure out delta phi between the start and end of the big drop column
# 						delphi=dpnewdat[bd_gates[nbd-1]]-dpnewdat[bd_gates[1]]
# #						print,'delphi ', delphi
# 					#Now add the extra correction of (a*-a)[phi(r2)-phi(r1)] since we already applied the first correction.
# 						zhdatcorr=zhdat(bd_gates[1]:ngates-1)
# 						wh_good=where(zhdatcorr != bad)
# #						IF wh_good[0] != -1 THEN BEGIN
# 						zdat_dum=zhdatcorr(wh_good)
# 						delphi_dum=fltarr(n_elements(wh_good))+delphi
# 						zhdatcorr(wh_good)=zdat_dum[*]+(a_star-a_coeff)*delphi_dum[*]
# 						zhdat(bd_gates[1]:ngates-1)=zhdatcorr
# 						zhcorr[*,qq,pp]=zhdat

# 						drdatcorr=drdat(bd_gates[1]:ngates-1)
# 						#STOP
# 						wh_gooddr=where(drdatcorr != bad)
# 						drdat_dum=drdatcorr(wh_gooddr)
# 						delphi_dum=fltarr(n_elements(wh_gooddr))+delphi
# 						drdatcorr(wh_gooddr)=drdat_dum[*]+(b_star-b_coeff)*delphi_dum[*]
# 						drdat(bd_gates[1]:ngates-1)=drdatcorr
# 						drcorr[*,qq,pp]=drdat
# 						#STOP
# 						bigdrop_flag[bd_gates[1]:ngates-1,qq,pp]=1
# 						bigdrop_flag[bd_gates[1]:bd_gates[nbd-1],qq,pp]=2
# 						nbd=0
# 						bd_gates=fltarr(1)
# 					ENDIF ELSE BEGIN
# #						print,'not enough points for big drop core'
# 						nbd=0
# 						bd_gates=fltarr(1)					
# 					ENDELSE # END IF / ELSE OVER number of points
# 				ENDELSE # END IF / ELSE over if /where there are bigdrops
# 			ENDFOR # END loop over finding locaitons of bigdrop cores
# 		ENDIF # End loop over if wh_bigdrops was successful in finding potential points
# 	ENDFOR # END loop over numrays		
# ENDFOR	#ENd loop over nsweeps

# 	print,'Max dr after BG correction: ', max(drcorr)

# 		dat={dz:zhcorr,$
# 			dr:drcorr,$
# 			flag:bigdrop_flag,$
# 			pts:tnumpts}

# RETURN,dat
# END






# # Below is the original IDL code, just going to keep there

# # function declaration here


# FUNCTION lc_atten_corr,dz,dr,dp, kd,rh,height, a_default,deva,b_default,devb, bad

# 	# first find where the data is good/bad
# 	wh_gd=where(dz != bad and dp != bad and kd ge 0.0 and kd le 2.0 and height le 2.5 and height ge 0.5 and dz lt 60.0)

# #First do DZ
# 	# defining the corrected Z and the specific attenuation
# 	zhcorr=dz*0.0+bad
# 	ah=dz*0.0+bad
	
# 	#STOP
# 	IF wh_gd[0] ne -1 THEN BEGIN
# 		zhdat=dz[wh_gd]
# 		dpdat=dp[wh_gd]
# 		continue_run=1
# 		IF n_elements(wh_gd) gt 200 THEN BEGIN
# 			fit=regress(dpdat,zhdat,chisq=chisq,const=const,correlation=corr,sigma=sigma,yfit=yfit)
# 			zhslope=fit
# #			print,zhslope
# 			numpts=n_elements(dpdat)
# 			WHILE abs(corr) lt 0.90 and numpts gt 500 and continue_run == 1 DO BEGIN
# 				diff=abs(zhdat-yfit)
# 				diffmom=moment(diff,sdev=sdev)
# 				wh_good=where(abs(diff) lt 2*sdev)    
# 				IF wh_good[0] NE -1 THEN BEGIN
# 				dp_new=dpdat(wh_good)
# 				zh_new=zhdat(wh_good)
# 				regress2=regress(dp_new,zh_new,chisq=chisq,const=const,correlation=corr,sigma=sigma,yfit=yfit)
# 				dpdat=dp_new
# 				zhdat=zh_new
# 				numpts=n_elements(dpdat)
# 				zhslope=regress2

# 	#			plot, dp_new,zh_new,psym=2

# 				ENDIF ELSE BEGIN
# 				continue_run=0
# 				ENDELSE
# 			ENDWHILE

# 			IF (corr^2 ge 0.25 AND numpts ge 200 AND sigma/(SQRT(numpts)) le 5.5 AND max(dpdat) ge 15 and zhslope[0] ge (a_default-deva) and zhslope[0] le (a_default+deva)) THEN BEGIN
# 				a_coeff=(-1)*(zhslope[0])
# 			ENDIF ELSE BEGIN
# 				print,'Slope is not representative of the overall slope! ', zhslope[0], 'Using default a value ', a_default
# 				a_coeff=(-1)*(a_default)				
# 			ENDELSE

# 		ENDIF ELSE BEGIN
# #			print,'Using default value due to not enough points!'
# 			a_coeff=(-1)*(a_default)
# 		ENDELSE
		
# 			wh_gdh=where(dz ne bad and dp gt 0.0)
# #			stop
# 			#
# 			IF wh_gdh[0] ne -1 THEN BEGIN
# 				zhcorr[wh_gdh]=dz[wh_gdh]+a_coeff*(dp[wh_gdh])
# 			ENDIF
			
# 			wh_gkdh=where(kd ne bad and kd gt 0.0)
# 			IF wh_gkdh[0] ne -1 THEN BEGIN
# 				ah[wh_gkdh]=kd[wh_gkdh]*a_coeff
# 			ENDIF 
# 	ENDIF ELSE BEGIN
# 		zhcorr=dz*0.0+bad
# 		ah=dz*0.0+bad
# 		a_coeff=!values.f_nan
# 	ENDELSE
# #stop
# ###Now for ZDR

# 	wh_gd=where(dz ne bad and dp ne bad and kd ge 1.0 and kd le 2.0 and height le 2.5 and height ge 0.5 and dz lt 60.0 and dr ne bad)

# 	zrcorr=dr*0.0+bad
# 	adp=dr*0.0+bad
# 	IF wh_gd[0] ne -1 THEN BEGIN
# 		zhdat=dr[wh_gd]
# 		dpdat=dp[wh_gd]
# 		continue_run=1
# 		IF n_elements(wh_gd) gt 200 THEN BEGIN
# 			fit=regress(dpdat,zhdat,chisq=chisq,const=const,correlation=corr,sigma=sigma,yfit=yfit)
# 			zhslope=fit
# #			print,zhslope
# 			numpts=n_elements(dpdat)
# 			WHILE abs(corr) lt 0.90 and numpts gt 500 and continue_run eq 1 DO BEGIN
# 				diff=abs(zhdat-yfit)
# 				diffmom=moment(diff,sdev=sdev)
# 				wh_good=where(abs(diff) lt 2*sdev)    
# 				IF wh_good[0] NE -1 THEN BEGIN
# 				dp_new=dpdat(wh_good)
# 				zh_new=zhdat(wh_good)
# 				regress2=regress(dp_new,zh_new,chisq=chisq,const=const,correlation=corr,sigma=sigma,yfit=yfit)
# 				dpdat=dp_new
# 				zhdat=zh_new
# 				numpts=n_elements(dpdat)
# 				zhslope=regress2

# 	#			plot, dp_new,zh_new,psym=2

# 				ENDIF ELSE BEGIN
# 				continue_run=0
# 				ENDELSE
# 			ENDWHILE

# 			IF (corr^2 ge 0.25 AND numpts ge 200 AND sigma/(SQRT(numpts)) le 5.5 AND max(dpdat) ge 15 and zhslope[0] ge (b_default-devb) and zhslope[0] le (b_default+devb)) THEN BEGIN
# 				b_coeff=(-1)*(zhslope[0])
# 			ENDIF ELSE BEGIN
# #				print,'Slope is not representative of the overall slope! ', zhslope[0], 'Using default a value ', b_default
# 				b_coeff=(-1)*(b_default)				
# 			ENDELSE

# 		ENDIF ELSE BEGIN
# #			print,'Using default value due to not enough points!'
# 			b_coeff=(-1)*(b_default)
# 		ENDELSE
		
# 			wh_gdh=where(dr ne bad and dp gt 0.0)
# 			IF wh_gdh[0] ne -1 THEN BEGIN
# 				zrcorr[wh_gdh]=dr[wh_gdh]+b_coeff*(dp[wh_gdh])
# 			ENDIF
			
# 			wh_gkdh=where(kd ne bad and kd gt 0.0)
# 			IF wh_gkdh[0] ne -1 THEN BEGIN
# 				adp[wh_gkdh]=kd[wh_gkdh]*b_coeff
# 			ENDIF 
# 	ENDIF ELSE BEGIN
# 		zrcorr=dr*0.0+bad
# 		adp=dr*0.0+bad
# 		b_coeff=!values.f_nan
# 	ENDELSE

# 	dat={dz:zhcorr,dr:zrcorr,ah:ah,adp:adp,acoeff:a_coeff,bcoeff:b_coeff}
# return,dat
# END


# FUNCTION big_drop,dz,dr,kd,dp,rh,height,a_star,b_star,a_coeff,b_coeff,bad
# ############################################################################
# #Apply a big-drop correction as described by Carey et al. 2000.
# # Here, a* and b* are derived from scattering simulations (as given in Carey et al. 2000)
# #We also have to limit the data where delta is large and rhohv is small but there is some significant amount of Kdp.
# #This technique requires identifiying the "big drop core" by looking for points wthat meet certain requirements, and then
# #applying a correct to all range gates down range of the core.
# ############################################################################

# print,'Doing Big drop correction!'
# bigdrop_flag=dz*0
# zhcorr=dz
# drcorr=dr

# #print,'Max of dr: ', max(drcorr)

# nsweeps=n_elements(dz[0,0,*])
# nrays=n_elements(dz[0,*,0])
# ngates=n_elements(dz[*,0,0])

# tnumpts=0
# 	FOR pp=0,nsweeps-1 DO BEGIN
# 		ray_correction=0
# 		FOR qq=0,nrays-1 DO BEGIN
# #			FOR kk=0,ngates-1 DO BEGIN
# 			zhdat=dz[*,qq,pp]
# 			rhdat=rh[*,qq,pp]
# 			kddat=kd[*,qq,pp]
# 			drdat=dr[*,qq,pp]
# 			dpnewdat=dp[*,qq,pp]
# 			hgt=height(*,qq,pp)
# 			wh_bigdrops=where(rhdat lt 0.80 and rhdat ne bad and kddat gt 1.0 and dpnewdat ne bad and hgt lt 5.0 and zhdat ne bad,complement=wh_notbigdrops)
	
# 			IF wh_bigdrops[0] ne -1 THEN BEGIN
# #				print,'Found big drops!'
# #				print,radar.volume[0].sweep[pp].ray[qq].h.azimuth
# 				nbd=1	
# 				bd_gates=fltarr(1)
# 			FOR chk=1,n_elements(wh_bigdrops)-1 DO BEGIN
# 				IF wh_bigdrops(chk) eq wh_bigdrops(chk-1)+1 THEN BEGIN
# 					nbd=nbd+1
# 					bd_gates=[bd_gates,wh_bigdrops(chk)]
# 				ENDIF ELSE BEGIN
# 					IF nbd GE 4 THEN BEGIN
# #						print,'big drop region found with ', nbd, ' Correcting'
# 						tnumpts=tnumpts+nbd
# 					#Now use the a_star and b_star corrections from Carey et al. 2000
# 					#Zhcorr(r)=Zh(r)+aphidp(r)+(a*-a)[phi(r2)-phi(r1)]
# 					#Need to figure out delta phi between the start and end of the big drop column
# 						delphi=dpnewdat[bd_gates[nbd-1]]-dpnewdat[bd_gates[1]]
# #						print,'delphi ', delphi
# 					#Now add the extra correction of (a*-a)[phi(r2)-phi(r1)] since we already applied the first correction.
# 						zhdatcorr=zhdat(bd_gates[1]:ngates-1)
# 						wh_good=where(zhdatcorr ne bad)
# #						IF wh_good[0] ne -1 THEN BEGIN
# 						zdat_dum=zhdatcorr(wh_good)
# 						delphi_dum=fltarr(n_elements(wh_good))+delphi
# 						zhdatcorr(wh_good)=zdat_dum[*]+(a_star-a_coeff)*delphi_dum[*]
# 						zhdat(bd_gates[1]:ngates-1)=zhdatcorr
# 						zhcorr[*,qq,pp]=zhdat

# 						drdatcorr=drdat(bd_gates[1]:ngates-1)
# 						#STOP
# 						wh_gooddr=where(drdatcorr ne bad)
# 						drdat_dum=drdatcorr(wh_gooddr)
# 						delphi_dum=fltarr(n_elements(wh_gooddr))+delphi
# 						drdatcorr(wh_gooddr)=drdat_dum[*]+(b_star-b_coeff)*delphi_dum[*]
# 						drdat(bd_gates[1]:ngates-1)=drdatcorr
# 						drcorr[*,qq,pp]=drdat
# 						#STOP
# 						bigdrop_flag[bd_gates[1]:ngates-1,qq,pp]=1
# 						bigdrop_flag[bd_gates[1]:bd_gates[nbd-1],qq,pp]=2
# 						nbd=0
# 						bd_gates=fltarr(1)
# 					ENDIF ELSE BEGIN
# #						print,'not enough points for big drop core'
# 						nbd=0
# 						bd_gates=fltarr(1)					
# 					ENDELSE # END IF / ELSE OVER number of points
# 				ENDELSE # END IF / ELSE over if /where there are bigdrops
# 			ENDFOR # END loop over finding locaitons of bigdrop cores
# 		ENDIF # End loop over if wh_bigdrops was successful in finding potential points
# 	ENDFOR # END loop over numrays		
# ENDFOR	#ENd loop over nsweeps

# 	print,'Max dr after BG correction: ', max(drcorr)

# 		dat={dz:zhcorr,$
# 			dr:drcorr,$
# 			flag:bigdrop_flag,$
# 			pts:tnumpts}

# RETURN,dat
# END



