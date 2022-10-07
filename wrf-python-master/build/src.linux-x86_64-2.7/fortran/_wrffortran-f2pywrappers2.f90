!     -*- f90 -*-
!     This file is autogenerated with f2py (version:2)
!     It contains Fortran 90 wrappers to fortran functions.

      
      subroutine f2pyinitwrf_constants(f2pysetupfunc)
      use wrf_constants, only : wrf_earth_radius
      use wrf_constants, only : rhowat
      use wrf_constants, only : t_base
      use wrf_constants, only : rad_per_deg
      use wrf_constants, only : thtecon1
      use wrf_constants, only : thtecon3
      use wrf_constants, only : thtecon2
      use wrf_constants, only : p1000mb
      use wrf_constants, only : rv
      use wrf_constants, only : cp
      use wrf_constants, only : rd
      use wrf_constants, only : abscoef
      use wrf_constants, only : celkel
      use wrf_constants, only : abscoefi
      use wrf_constants, only : eslcon2
      use wrf_constants, only : eslcon1
      use wrf_constants, only : pi
      use wrf_constants, only : tlclc2
      use wrf_constants, only : tlclc3
      use wrf_constants, only : rho_g
      use wrf_constants, only : tlclc1
      use wrf_constants, only : tlclc4
      use wrf_constants, only : deg_per_rad
      use wrf_constants, only : cpmd
      use wrf_constants, only : sclht
      use wrf_constants, only : ussalr
      use wrf_constants, only : default_fill
      use wrf_constants, only : rho_s
      use wrf_constants, only : rho_r
      use wrf_constants, only : alpha
      use wrf_constants, only : celkel_triple
      use wrf_constants, only : ezero
      use wrf_constants, only : expon
      use wrf_constants, only : rgasmd
      use wrf_constants, only : g
      use wrf_constants, only : errlen
      use wrf_constants, only : eps
      use wrf_constants, only : gamma_seven
      use wrf_constants, only : gammamd
      use wrf_constants, only : algerr
      use wrf_constants, only : gamma
      use wrf_constants, only : exponi
      external f2pysetupfunc
      call f2pysetupfunc(wrf_earth_radius,rhowat,t_base,rad_per_deg,thte&
     &con1,thtecon3,thtecon2,p1000mb,rv,cp,rd,abscoef,celkel,abscoefi,es&
     &lcon2,eslcon1,pi,tlclc2,tlclc3,rho_g,tlclc1,tlclc4,deg_per_rad,cpm&
     &d,sclht,ussalr,default_fill,rho_s,rho_r,alpha,celkel_triple,ezero,&
     &expon,rgasmd,g,errlen,eps,gamma_seven,gammamd,algerr,gamma,exponi)
      end subroutine f2pyinitwrf_constants

