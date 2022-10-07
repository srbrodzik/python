c  HP 9000 wind synthesis program

C     Written by mostly by Dave Jorgensen by Feb 1994
c        Brad Smull wrote the "old" vertical integration and O'Brien
c        vertical divergence error adjustment subroutines
C        Tom Matejka wrote the subroutines that calculate the
c        overdetermined u,v wind calculations and the 
c        vertical w calculation (variational O'Brien routines)
C     NOAA/ERL
C     National Severe Storms Laboratory
C     Mesoscale Research Division
C     Boulder, Colorado  80303
      
c   Written originally for an HP 1000/A900 in Sept 1984 by DJ
c     Various modules added to enhance vertical velocity calc by BS
c     April, 1990, changes from Idres to Windsyn made by RH
C     June, 1993 modified for an HP 9000 by DJ (Windsyn)
c     Added Tom Matejka's routines Jan 1994 by DJ
C
C     Windsyn can work with single-Doppler (one input file) or
C     dual-Doppler synthesis of winds (two input files). Divergence
C     fields and resulting vertical velocities may be calculated using
C     upward or downward vertical integration of an elastic continuity.
C     If desired, boundary conditions may be enforced using an O'Brien-
C     type adjustment which distributes error uniformly with height.
C     The program can handle any combination of P-3 legs and ground
c     -based
C     observations.
C     
C       The input files to this program are those files written by programs
C     fast_Interp and Grnd_Interp.  The output files (with .HDR and .DPW
C     extensions) can be read and plotted with programs Pls_Retr and
C     Pls_Retr_Xsecn (among others).
C
c       Listings are written to a .log file

C * * * * * * * * * The contents of the .prm file: * * * * * * * * * * *

C Line 1: Istyle, Nrdrs, ThresH, ThresV     (*)
C   Istyle is an integer which determines which kind of dual-doppler
C       or single-doppler analysis will be done.

c Istyle  Radar1  Radar2 #  u,v calcul.        w calculation         O'Brien?
c ------  ------ ------- - -------------      ---------------       ----------
c   1      P-3     P-3   2 "old" L-tracks    "old" div integration    0=n 1=y
c   2      P-3     P-3   2 "old" FAST        "old" div integration    0=n 1=y
c   3      P-3     P-3  >2 Matejka triple        directly               N/A
c   4      P-3     P-3  >1 Matejka dual      "new" div integration    0=n 1=y
c   5      P-3     P-3  >2 Hybrid (2 or 3)   "new" div integration    0=n 1=y
c   6      P-3     P-3  >2 Hybrid (3 else 2) "new" div integration    0=n 1=y
c   7      P-3     P-3  >2 Full Variational     variational           0=n 1=y
c   8      P-3     P-3  >2 Hybrid (2 or 3)   "old" div integration    0=n 1=y
c   9      P-3     P-3  >1 Matejka dual      "old" div integration    0=n 1=y
c  10    Ground    P-3   2 Conven. dual      "old" div integration    0=n 1=y
c  11    P-3     Ground >2 Matejka dual      "old" div integration    0=n 1=y
c -11     P-3    P-3    >2 like 11 but only uses ground-radar reflectivity
c  12     P-3    P-3    >2 like 2 but uses ground-radar reflectivity
c  20    Ground  Ground  2 Conven. dual      "old" div integration    0=n 1=y
c  30    Single Grnd Rdr 1    N/A                  N/A                  N/A
c  40    Single P-3 Leg  1    N/A                  N/A               N/A

c   Nrdrs is the number of radars used in the analysis
c    each radar needs a corresponding file

c   ThresH is the error variance threshold used to get rid
c   of bad horizontal wind data
c   ThresV is the error variance threshold used to get rid
c   of bad vertical wind data

c Line 2: Iw_at_top, Iw_0

c   Iw_at_Top is a flag that is used in the variational adjustment routine
c      to specify what is done when the routine fails to produce a 
c      column of vertical velocities.  If Iw_at_top is 1 then a second
c      attempt is made to integrate the failed column using the hybrid approach. 
c      If Iw_at_top is 0 then no hybrid approach is tried if a column fails
c      the adjustment procedure.

c   Iw_0 is a flag that is used in both the variational and hybrid routines
c      to specify what is done with those routines fail to produce a column of
c      vertical velocities.  If Iw_0 = 1, then a second attempt is made to 
c      integrate the failed column using w=0 as a top boundry condition.  If both
c      Iw_at_top and Iw_0 are 1 then the hybrid approach is tried first and if
c      that attempt also failed then the w=0 approach is tried.

c   Top_Hgt is the lowest acceptable level for cloud top height used in the hybrid
c      scheme (km)

C Line 3: Idbz_Parm, VT_Snow, VT_Rain, IHole_Fill, Vert Interp

C   Idbz_parm is a flag which determines how the dBZ will be calculated:
C         -1: Maximum of radars 1 to Nrdrs
C          n: Use radar n
C         -2: Average of radars 1 to Nrdrs
C         -3: Difference radar 1 - radar 2
C       If only a single radar analysis is being done, then
C       Idbz_parm will forced to be 1 (only radar available)

C   VT_Snow, VT_Rain: optional floating point values for the heights (in km)
C       of the transition zones where snow starts to melt and where the
C       snow is completely melted to rain.  These values are important in
C       calculating terminal fall velocities.  If VT_Snow and VT_Rain are
C       not included, then default values will be used.

c   Ihole_Fill:  Horizontal hole filling flag.
c                = 1 means fill holes, =0 means no hole filling.  Horizontal 
c                Hole filling parameters are controlled in line 4.

c   Vert Interp: Number of vertical levels to be interpolated if holes exist
c                in the vertical column.

C Line 4: Nsmth, IsmTyp, Zsfc, Int_Dir, Kup, Iobr, Iters

c   Nsmth = number of smooth iterations =0 no smoothing

c   IsmTyp = Type of smoothing desired:
c         0 = No smoothing done
c         1 = Leise
c         2 = Binomial

C   Zsfc = Height (in km) of the surface.  If over the sea Zsfc should be 0.
c          This is the height at which the lower boundary of w=0 is applied
c          if an O'Brien divergence error correction is applied.

c   Int_Dir = Direction of integration if no variational or O'Brien correction
c             is desired.  =1 means upwards, -1 means downwards.

c   Kup = highest height level that will be used to extrapolate divergence
c           downward if divergence is missing on level #1.

c   Iobr = O'Brien divergence correction =0 no; =1 yes

c   Iters = number of iterations on the wind calculation

c Line 5: Z1, Tv1, Z2, Tv2, Z3, Rho3
c         Z1 is an arbitrary height (m) where TV1 (virtual temperature
c            in K) is defined
c         TV1 = virtual temperature (K) at Z1
c         Z2  = height in meters
c         TV2 = virtual temperature (K) at Z2
c         Z3  = height in meters
c         Rho3= density (g/kg) at Z3

c Line 6: Max_Cons_Unocc_Octants, Max_Search_Radius,
c         Min_Values_Per_Octant, Max_Values_Per_Octant

c Line 7: Output file name (a)

C Lines 8 to 8+Nrdrs: file                              (a50)
C   file is the input file name (full path name) of the ith radar

c Example .prm file:
c  7 4 4.00 3.00      ! solution type, # of radars, hor thres, ver thres
c  0 0 1.5            ! Iw_at_top, Iw_0, Top_Hgt
c  -1 4.75 4.25 1 2   ! dbz, snow hgt, rain hgt, fill holes?, Vert Interp
c  2 2 0.0 -1 2 1 1   ! # smooths, sm type, Zsfc, Idir, Kup, O'Brien?, iters
c  1500.0 296.8 5550.0 269.3 3000.0 0.909  ! Z1, Tv1, Z2, Tv2, Z3, Rho3  
c  3 5 1 2 ! max_cons_unoc_oct, max_srch_rad, min_val_per_oct, max_val_per_oct
c  test                                    ! Output file name
c  /disk2/doppler/toga/feb22/2057h.for     ! 1st input file
c  /disk2/doppler/toga/feb22/2057h.aft     ! 2nd input file
c  /disk2/doppler/toga/feb22/2052i.for     ! 3rd input file
c  /disk2/doppler/toga/feb22/2052i.aft     ! 4th input file

C  The contents of Data are as follows:
C          Data(Imax,Jmax,Kmax,1) = East-West Velocity [m/s] (u)
C          Data(Imax,Jmax,Kmax,2) = North-South Velocity [m/s] (v)
C          Data(Imax,Jmax,Kmax,3) = Max Reflectivity from any radar [dBZ]
C          Data(Imax,Jmax,Kmax,4) = Vertical Velocity [m/s]
C          Data(Imax,Jmax,Kmax,5) = Horizontal Divergence [per s]
C          Data(Imax,Jmax,Kmax,6) = u standard error [m/s] 
C          Data(Imax,Jmax,Kmax,7) = v standard error [m/s]
C          Data(Imax,Jmax,Kmax,8) = w standard error [m/s]
C          Data(Imax,Jmax,Kmax,9) = terminal fall velocity [m/s]
C          Data(Imax,Jmax,Kmax,10)= Div standard error [s-1]
C          Data(Imax,Jmax,Kmax,11)= Vt standard error [m/s]
C          Data(Imax,Jmax,Kmax,12)= Max time difference between any radars
C          Data(Imax,Jmax,Kmax,13)= Average time of all radars

c  Note that Divergence and terminal fall velocity are adjusted values
c  if the variational adjustment technique is used

c  The contents of Radar are as follows:
c          Radar(1,Imax,Jmax,Kmax,Nrdrs) = Azimuth from North [deg]
c          Radar(2,Imax,Jmax,Kmax,Nrdrs) = Radial Velocity [m/s]
c          Radar(3,Imax,Jmax,Kmax,Nrdrs) = Elevation angle [deg]
c          Radar(4,Imax,Jmax,Kmax,Nrdrs) = Reflectivity [dBZ]
c          Radar(5,Imax,Jmax,Kmax,Nrdrs) = Time difference from nominal
c          Radar(6,Imax,Jmax,Kmax,Nrdrs) = Range to target (m)

c     Dynamic memory for big arrays

C      Subroutine windsyn(Istyle, Nrdrs, ThresH, ThresV, 
C     #          Iw_at_top, Iw_0, 
C     #          Idbz_Parm, VT_Snow, 
C     #          VT_Rain, IHole_Fill, Vert Interp, 
C     #          Nsmth, IsmTyp, Zsfc, Int_Dir, Kup, Iobr, Iters, 
C     #          Z1, Tv1, Z2, Tv2, Z3, Rho3, 
C     #          Max_Cons_Unocc_Octants, Max_Search_Radius, 
C     #          Min_Values_Per_Octant, Max_Values_Per_Octant, 
C     #          file, Radar, SynthDat)
c      Subroutine windsyn(Imax, Jmax, Kmax, Nrdrs)
      Subroutine (Radar, SynthDat)
      
      Real,Dimension(:,:,:,:,:)::Radar
      Real,intent(inout),Dimension(:,:,:,:)::SynthDat
c      Real,Dimension(7,Imax,Jmax,Kmax,Nrdrs)::Radar
c      Real,intent(inout),Dimension(Imax,Jmax,Kmax,13)::SynthDat
Cf2py intent(out) SynthDat
Cf2py intent(hide) Radar
      Parameter (LuPar = 90)    ! lu to use for parameter file
      Parameter (MxR = 8)       ! maximum number of grnd radars + flight legs

      Dimension Elapsed(2)      ! for elapsed cpu time

      Dimension Imax1(MxR), Jmax1(MxR), Kmax1(MxR), Sx1(MxR), Sy1(MxR),
     #          Sz1(MxR), Olat1(MxR), Olon1(MxR), Z01(MxR), Nmosm1(MxR),
     #          Su1(MxR), Sv1(MxR), Xniq1(MxR), Rlat1(MxR), Rlon1(MxR),
     #          Rotcr(MxR), Eloff(MxR), Cant(MxR), RangeMax(MxR),
     #          Alt_Flag(MxR)
      Dimension Itime_Limits1(12,MxR), Init_Time1(6,MxR)

      Character Time_Current1(MxR)*24, Flid1(MxR)*8, Project1(MxR)*16,
     #          Namdf1(MxR)*100, Type1(MxR)*3, Nameif(MxR)*100

      Character File*100, Namof*100, Line*100
      Character*3 Typex1, Typex2, Typex3
      Character Real_Time*24, Smoothing(0:2)*10

      Integer*4 FirstChar
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, 
     #     Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /Sounding/ Z1, Tv1, Z2, Tv2, Z3, Rho3
      Common /Limits/ Zsfc, Iters, Kup, Klim, Top_Hgt
      Common /Print/ LuMes
      Common /HoleFillParms/ Max_Cons_Unoccupied_Octants,
     #                       Max_Search_Radius,
     #                       Min_Values_Per_Octant, 
     #                       Max_Values_Per_Octant
      Common /Debug/ Idebug, Ix1, Ix2, Jy1, Jy2
      Common /w_Techniques/ Iw_at_top, I_w_0


c  Defaults for some variables:

      Data VT_SnowDef /4.75/
      Data VT_RainDef /4.25/
      Data Smoothing /'None','Leise', 'Binomial'/

c  Number of fields in the array "SynthDat"

      Nflds = 13

C    Open the data file containing the parameter information.

      n = IargC()

      If (n .lt. 1) Then
         Write (6,'("Enter prm file name:",$)')
         Read (5,'(a)') File
      Else
         Call GetArg(1,File)
      End If

C  Open the parameter file
  
      Open (LuPar,Err=1020,File=File,Iostat=Ierr,Status='Old')

      Write (6,'(/" P-3 Doppler Wind Synthsis program - .prm File:",a)') File
 
c  Create and open the .log file with the same name as the .prm file

      nc = LastChar(File,'.')
      File = File(1:nc) // 'log'

      If (LuMes .ne. 6) Then
         Open (LuMes,File=File,Iostat=Ierr,Err=1020)
      End If

c  Create and open the .err file with the same name as the .prm file

      File = File(1:nc) // 'err'
      Open (7,File=File,Iostat=Ierr,Err=1020)

      Call Fdate(Real_Time)

      Write (LuMes,'(/20x," Program Windsyn"
     > /"Batch Doppler Radar Wind Synthesis started at ",a)') Real_Time
      Write (7,'(/20x," Program Windsyn"
     > /"Batch Doppler Radar Wind Synthesis started at ",a)') Real_Time

      Read (LuPar,'(a)') Line
      Read (Line,*,End=999,Err=888,Iostat=Ierr)
     #     Istyle, Nrdrs, ThresH, ThresV

      Read (LuPar,'(a)') Line
      Read (Line,*,End=999,Err=888,Iostat=Ierr)
     #     Iw_at_top, I_w_0, Top_Hgt

      If (Top_Hgt .lt. 0.0) Top_Hgt = 0.0
      If (Iw_at_top .lt. 0) Iw_at_top = 0
      If (Iw_at_top .gt. 1) Iw_at_top = 1
      If (I_w_0 .lt. 0) I_w_0 = 0
      If (I_w_0 .gt. 1) I_w_0 = 1

      Idebug = 0

      Write (LuMes,'(/"Style=",i3," Number of radars=",i2,
     #     " Threshold (horiz)=",f4.1," Threshold (vert)=",f4.1)')
     #      Istyle, Nrdrs, ThresH, ThresV

      If (Istyle .eq. 1) Then   ! two perpendicular flight legs
         Typex1 = 'LG1'
         Typex2 = 'LG2'

         If (Nrdrs .ne. 2) Then
            Write (LuMes,100) Istyle, Nrdrs
 100        Format ('# of radars needs to be 2 for Style=',i3,
     #              ' specified:',i4)
            Stop
         End If

      Else If (Istyle .eq. 2) Then     ! single FAST leg
         Typex1 = 'for'
         Typex2 = 'aft'

         If (Nrdrs .ne. 2) Then
            Write (LuMes,100) Istyle, Nrdrs
            Stop
         End If

      Else If (Istyle .eq. 3) Then     ! QUAD (two AC)
         Typex1 = 'for'
         Typex2 = 'aft'

         If (Nrdrs .lt. 3) Then
            Write (LuMes,101) Istyle, Nrdrs
 101        Format ('# of radars needs to be >2 for Style=',i3,
     #              ' you specified:',i2)
            Stop
         End If

      Else If (Istyle .eq. 4) Then     ! Matejka 2-component solution
         Typex1 = 'for'
         Typex2 = 'aft'

         If (Nrdrs .lt. 2) Then
            Write (LuMes,102) Istyle, Nrdrs
 102        Format ('# of radars needs to be at least 2 for Style=',i3,
     #              ' you specified:',i4)
            Stop
         End If

      Else If (Istyle .eq. 5) Then     ! Matejka 2 or 3 "quad" solution
         Typex1 = 'for'
         Typex2 = 'aft'

         If (Nrdrs .lt. 3) Then
            Write (LuMes,101) Istyle, Nrdrs
            Stop
         End If

      Else If (Istyle .eq. 6) Then     ! Matejka 2 else 3 "quad" solution
         Typex1 = 'for'
         Typex2 = 'aft'

         If (Nrdrs .lt. 3) Then
            Write (LuMes,101) Istyle, Nrdrs
            Stop
         End If

      Else If (Istyle .eq. 7) Then     ! Matejka variational adjustment
         Typex1 = 'for'
         Typex2 = 'aft'

         If (Nrdrs .lt. 2) Then
            Write (LuMes,101) Istyle, Nrdrs
            Stop
         End If

      Else If (Istyle .eq. 8) Then     ! Matejka 2 or 3 "quad" solution
         Typex1 = 'for'                ! but with "old" vertical velocity
         Typex2 = 'aft'

         If (Nrdrs .lt. 3) Then
            Write (LuMes,101) Istyle, Nrdrs
            Stop
         End If

      Else If (Istyle .eq. 9) Then     ! Matejka 2-component solution
         Typex1 = 'for'                ! but with "old" vertical velcity
         Typex2 = 'aft'

         If (Nrdrs .lt. 2) Then
            Write (LuMes,102) Istyle, Nrdrs
            Stop
         End If

      Else If (Istyle .eq. 10) Then     ! P-3 + ground radar volume
         Typex1 = 'p3'
         Typex2 = 'gnd'

         If (Nrdrs .ne. 2) Then
            Write (LuMes,100) Istyle, Nrdrs
            Stop
         End If

      Else If (Abs(Istyle) .eq. 11) Then     ! 2 P-3 beams + ground radar volume
         Typex1 = 'for'
         Typex2 = 'aft'
         Typex3 = 'gnd'

         If (Nrdrs .lt. 3) Then
            Write (LuMes,100) Istyle, Nrdrs
            Stop
         End If

      Else If (Istyle .eq. 12) Then     ! 2 P-3 beams + ground radar volume
         Typex1 = 'for'
         Typex2 = 'aft'
         Typex3 = 'gnd'

         If (Nrdrs .lt. 2) Then
            Write (LuMes,100) Istyle, Nrdrs
            Stop
         End If

      Else If (Istyle .eq. 20) Then     ! two ground radar volumes
         If (Nrdrs .ne. 2) Then
            Write (LuMes,100) Istyle, Nrdrs
            Stop
         End If

         Typex1 = 'gnd'
         Typex2 = 'gnd'
         Typex3 = 'gnd'
      Else If (Istyle .gt. 30) Then
         If (Nrdrs .gt. 1) Then
            Write (LuMes,103) Istyle, Nrdrs
 103        Format ('# of radars should be 1 for style:',i4,
     #              ' you specified:',i4)
            Stop
         End If
      Else
         Write (LuMes,'("Istyle out of range!!!",i5)') Istyle
         Stop
      End If

c  Input the control parameters:

      Read (LuPar,'(a)') Line
      Read (Line,*,End=999,Err=888,Iostat=Ierr) 
     #     IdBZ_Parm, Vt_Snow, Vt_Rain, Ihole_Fill, Klim

      Read (LuPar,'(a)') Line
      Read (Line,*,End=999,Err=888,Iostat=Ierr)
     #     Nsmth, IsmTyp, Zsfc, Idir_Int, Kup, Iobr, Iters

c     At least 1 iteration?

      If (Iters .lt. 1) Then
         Write (LuMes,'(/" Must have at least 1 iteration")')
         Iters = 1
      End If

c  Suiable Vt relationships?

      If (VT_Snow .lt. VT_Rain) Then
          Write (LuMes,'(/" Your snow and rain heights will be switched"
     >   ,"!!" /)')
          Dummy = VT_Snow
          VT_Snow = VT_Rain
          VT_Rain = Dummy
      End If

c  Use default relationships?

      If (Vt_Snow .le. 0.0) VT_Snow = VT_SnowDef
      If (Vt_Rain .le. 0.0) VT_Rain = VT_RainDef
      If (IsmTyp .lt. 0 .or. IsmTyp .gt. 2) IsmTyp = 0

c  direction of integration?

      If (Idir_Int .le. 0) Idir_Int = -1
      If (Idir_Int .gt. 0) Idir_Int = 1

c  Force dBZ parameter for single radar

      If (Istyle .ge. 30) Then
        IdBZ_Parm = 1
      End If

c  Horizontal hole filling?

      If (Ihole_Fill .lt. 0 .or. Ihole_Fill .gt. 1) Ihole_Fill = 0

c  Vertical hole filling?

      If (Klim .lt. 0) Klim = 0

c  Get the sounding information for the variatinal approach

      Read (LuPar,'(a)') Line
      Read (Line,*,End=999,Err=888,Iostat=Ierr)
     #     Z1, Tv1, Z2, Tv2, Z3, Rho3

c  Use the defaults (GATE mean sounding?)

      If (Z1 .le. 0.0) Then
         Z1 = 1500.0            ! 850 mb
         Tv1 = 23.7 + 273.16
         Z2 = 5550.0            ! 500 mb
         Tv2 = -3.9 + 273.16
         Z3 = 3000.0            ! 700 mb
         Rho3 = 0.909           ! kg/m3
      End If

c  Input the hole filling parameters

      Read (LuPar,'(a)') Line
      Read (Line,*,End=999,Err=888,Iostat=Ierr)
     #     Max_Cons_Unoccupied_Octants,
     #     Max_Search_Radius,
     #     Min_Values_Per_Octant,
     #     Max_Values_Per_Octant

      If (Max_Cons_Unoccupied_Octants .lt. 0) Then
         Max_Cons_Unoccupied_Octants = 3
         Max_Search_Radius = 5
         Min_Values_Per_Octant = 1
         Max_Values_Per_Octant = 2
      End If

      Write (LuMes,'(//"For the hole filling routine:"/
     #     "Max_Cons_Unoccupied_Octants=",i3/
     #     "Max_Search_Radius=",i3/
     #     "Min_Values_Per_Octant=",i3/
     #     "Max_Values_Per_Octant=",i3)') 
     #     Max_Cons_Unoccupied_Octants,
     #     Max_Search_Radius, Min_Values_Per_Octant, 
     #     Max_Values_Per_Octant

c  Read the name of the output wind and header files

      Read (LuPar,'(a)') Line
      Read (Line,'(a)',End=999,Err=888,Iostat=Ierr) Namof

      nch = FirstChar(Namof,' ') - 1
      File = Namof(1:nch)
      Namof = File

      Write (LuMes,'(/" IdbZ_Parm:",i2/,
     #                " Vt_Snow:",f7.2/,
     #                " Vt_Rain:",f7.2/,
     #                " Horizontal hole filling:",i3/,
     #                " Vertical hole filling (max holes):",i3/,
     #                " Modes for alternative w calc:"/,8x,
     #                    "  Use hybrid (for Variational):",i2/,8x,
     #                    "  Use w=0 at top (for Var and Hybrid):",i2/,
     #                " Lowest height for cloud top (Hybrid):",F6.1/,
     #                " # of smooths:",i3/,
     #                " Smoothing Type:",i3,1x,a/,
     #                " Zsfc, Idir_Int, Kup:",F7.2,1x,2i2/,
     #                " OBrien correction:",i3/,
     #                " Iterations on wind calculation:",i3/,
     #                " Output file name: ",a/
     #                " Z1: ",F6.1/
     #                " Tv1: ",F6.1/
     #                " Z2 :",F6.1/
     #                " Tv2: ",F6.1/
     #                " Z3: ",F6.1/
     #                " Rho3: ",F6.1)')
     #                        IdBZ_Parm, Vt_Snow, Vt_Rain, Ihole_Fill,
     #                        Klim, Iw_at_top, I_w_0, Top_Hgt,
     #                        Nsmth, IsmTyp, Smoothing(IsmTyp),
     #                        Zsfc, Idir_Int, Kup,
     #                        Iobr, Iters, Namof,
     #                        Z1, Tv1, Z2, Tv2, Z3, Rho3

c  Open and read the headers for the files  

      Do i = 1, Nrdrs
         Lu = 40 + i
         Read (LuPar,'(a)') Line
         Read (Line,'(a)',End=999,Err=888,Iostat=Ierr) Nameif(i)

         nch = FirstChar(Nameif(i),' ') - 1
         File = Nameif(i)(1:nch)
         Nameif(i) = File(1:nch)

         Open (Lu,Err=999,File=File,Iostat=Ierr,
     >         Form='UNFORMATTED',Status='Old')

         Write (LuMes,'(/" Input File Name:",a)') Nameif(i)

         Read (Lu) Time_Current1(i), Imax1(i), Jmax1(i), Kmax1(i),
     #       Sx1(i), Sy1(i), Sz1(i), Flid1(i), Olat1(i), Olon1(i),
     #       Z01(i), (Itime_Limits1(l,i),l=1,12), Nmosm1(i),
     #       (Init_Time1(l,i),l=1,6), Su1(i), Sv1(i), Project1(i),
     #       Rotcr(i), Eloff(i), Cant(i), Namdf1(i), Xniq1(i),
     #       Rlat1(i), Rlon1(i), Type1(i), RangeMax(i), Alt_Flag(i)

         Write (LuMes,12) Time_Current1(i), Imax1(i), Jmax1(i),
     #     Kmax1(i),
     #     Sx1(i), Sy1(i), Sz1(i), Flid1(i), Olat1(i), Olon1(i),
     #     Z01(i),
     #     (Itime_Limits1(l,i),l=1,12), Nmosm1(i),
     #     (Init_Time1(l,i),l=1,6),
     #     Su1(i), Sv1(i), Project1(i), Namdf1(i),
     #     Xniq1(i), Rlat1(i), Rlon1(i), Type1(i), RangeMax(i),
     #     Alt_Flag(i)
    
  12  Format (/' File created on: ',a/
     #         ' Imax Jmax Kmax: ',3i4/
     1         ' Sx Sy Sz [km]: ',3F7.3/
     #         ' Flight Id: ',a/
     #         ' Olat Olon Z0 :',2F8.3,F5.1/
     #         ' Start, End Times: ',i4.4,"/",2i2.2,1x,3i2.2,1x,i4.4,"/",2i2.2,1x,3i2.2/
     #         ' Nmosm:',I1,1x,'Init Time: ',i4.4,"/",2i2.2,1x,3i2.2/
     #         ' Storm u,v components: ', 2f7.2/
     #         ' Project: ',a/
     #         ' Name of .cor file: ',a/
     #         ' Nyquist: ',f7.2/
     #         ' Radar Position: ',2f8.2,' Type: ',a/
     #         ' Maximum Radar Range:', F8.2/
     #         ' Altitude Flag:', F5.1)
      End Do

c  Any debug columns to print out?

      If (Idebug .eq. 1) Then
         Read (LuPar,'(a)') Line
         Read (Line,*,End=999,Err=888,Iostat=Ierr)
     #     Ix1, Ix2, Jy1, Jy2
      End If

      Close (LuPar)

C  Since by definition grids must match, adopt values from R.1 file
  
      Imax = Imax1(1)
      Jmax = Jmax1(1)
      Kmax = Kmax1(1)
      Olat = Olat1(1)
      Olon = Olon1(1)
      Z0 = Z01(1)
      Sx = Sx1(1)
      Sy = Sy1(1)
      Sz = Sz1(1)
      Nmosm = Nmosm1(1)

      If (Nmosm .eq. 1) Then
         Su = Su1(1)
         Sv = Sv1(1)
      Else
         Su = 0.
         Sv = 0.
      End If

      If (Type1(1) .ne. Typex1 .and. Type1(1) .ne. Typex2 .and. Type1(1)
     $     .ne. Typex3) Then
         Write (LuMes,'("Input file is:",a,", not ",a," or ",
     #        a," or ",a)') Type1(1), Typex1, Typex2, Typex3
         Stop 777
      End If

      If (Type1(2) .ne. Typex1 .and. Type1(2) .ne. Typex2 .and. Type1(2)
     $     .ne. Typex3) Then
         Write (LuMes,'("Input file is:",a,", not ",a," or ",
     #        a," or ",a)') Type1(2), Typex1, Typex2, Typex3
         Stop 777
      End If

c  Compare the grids  

      If (Istyle .lt. 30) Then
         Do i = 1, Nrdrs-1
            If (Imax1(i) .ne. Imax1(i+1) .or.
     #          Jmax1(i) .ne. Jmax1(i+1) .or.
     #          Kmax1(i) .ne. Kmax1(i+1) .or.
     #          z01(i)   .ne. z01(i+1)   .or.
     1          Olat1(i) .ne. Olat1(i+1) .or. 
     #          Olon1(i) .ne. Olon1(i+1) .or.
     #          Alt_Flag(i) .ne. Alt_Flag(i+1)) Then
               Write (LuMes,'(/"Grids are not compatable !")')
               Stop
            End If
            
            Do k = 1, 6
               If (Init_Time1(k,i) .ne. Init_Time1(k,i+1)) Then
                  Write (LuMes,'(/"Storm motions not compatable !")')
                  Stop
               End If
            End Do
         End Do
      End If

c  Vertical velocity integration parameters within range?

      If (Zsfc .gt. Z0) Then
         Write (7,'(/"Zsfc must be less than Z0!!",
     #        " Zsfc=",F7.2," Z0=",F7.2"  Abort!")') Zsfc, Z0
         Stop
      End If

c  Create the header (.hdr) file

      nc = FirstChar(Namof,' ') - 1
      Namof = Namof(1:nc) // '.hdr'  

      Write (LuMes,'(/" Header (.hdr) File Name:",a)') Namof(1:LenTrim(Namof))

      File = Namof
      Open (3,Err=999,File=File,Iostat=Ierr,Form='UNFORMATTED', Convert='Big_Endian')
  
c  Write out the header info for the wind (.hdr) file
  
      Call Fdate(Real_Time)
  
c  Build the name of the .dpw file

      Idum = -999
      nc = LenTrim(Namof)
      Namof = Namof(1:nc-3) // 'dpw'

      Write (LuMes,'(" Output (.dpw) File Name:",a)') Namof(1:LenTrim(Namof))

      nc = LastChar(Namof,'/')
      Line = Namof(nc+1:Lentrim(Namof))

c  Write out the header information
  
      Write (3) Real_Time, Imax, Jmax, Kmax, Sx, Sy, Sz,
     #          Olat, Olon, Z0, Nrdrs, Nmosm, Su, Sv, IdBZ_Parm,
     #          Idir_Int, Istyle, Nsmth, Iobr, Line(1:50),
     #          Baddata, Vt_Snow, Vt_Rain, ThresH, ThresV, IsmTyp,
     #          Ihole_Fill, Klim, Iw_at_top, I_w_0,
     #          Top_Hgt, IDum, IDum, IDum, IDum, IDum

      Do i = 1, Nrdrs

c Straighten out the year

         If (Itime_Limits1(1,i) .gt. 1999) Then
            Itime_Limits1(1,i) = Itime_Limits1(1,i) - 2000
            Itime_Limits1(7,i) = Itime_Limits1(7,i) - 2000
         End If

         nc = LastChar(Nameif(i),'/')
         Nameif(i) = Nameif(i)(nc+1:Lentrim(Nameif(i)))
         nc = LastChar(Namdf1(i),'/')
         Namdf1(i) = Namdf1(i)(nc+1:Lentrim(Namdf1(i)))
         
         Write (3) Time_Current1(i), Nameif(i)(1:50),
     #        Imax1(i), Jmax1(i), Kmax1(i),
     #        Sx1(i), Sy1(i), Sz1(i), Flid1(i), Olat1(i), Olon1(i),
     #        Z01(i), (Itime_Limits1(l,i),l=1,12), Nmosm1(i),
     #        (Init_Time1(l,i),l=1,6), Su1(i), Sv1(i), Project1(i),
     #        Rotcr(i), Eloff(i), Cant(i), Namdf1(i)(1:50), Xniq1(i),
     #        Rlat1(i), Rlon1(i), Type1(i), RangeMax(i), Alt_flag(i),
     #        IDum, IDum 
      End Do

      Close (3)
  
c  This file will be direct access to facilitate reading and plotting
c  the SynthDat. Six fields are specified:
c     1 - U(Imax,Jmax,Kmax)   East - West relative velocity [m/s]
c     2 - V(Imax,Jmax,Kmax)   North - South relative velocity [m/s]
c     3 - W(Imax,Jmax,Kmax)   Vertical velocity [m/s]
c     4 - D(Imax,Jmax,Kmax)   Divergence x 1000 [s-1]
c     5 - Z(Imax,Jmax,Kmax)   Reflectivity [dBz]
c     6 - B(Imax,Jmax,Kmax)   Terminal fall velocity [m/s]
  
C  Open the output file

      Lrec = Imax * 4 
      File=Namof

      Open (3, Err=999, File=File, Iostat=Ierr, Access='Direct',
     #      Recl=Lrec, Convert='Big_Endian')

c  Load up the original radar parameters

      Do n = 1, Nrdrs
         Lu = 40 + n
         Call Ldpln (Lu, Radar, n)
         Close (Lu)
      End Do

c  Initialize the data arrays to "Baddata"

      Call Init_Field (SynthDat)

c  Perform the wind synthesis

      If (Istyle .le. 2) Then                          ! "old" dual approach
         Call Synthesis (SynthDat, Radar, Type1)

      Else If (Istyle .eq. 3) Then                     ! Matejka triple
         Call Triple (SynthDat, Radar, Type1)

      Else If (Istyle .ge. 4 .and. Istyle .le. 7) Then ! Matejka approach
         Call Synthesis_New (SynthDat, Radar, Type1)

      Else
         Call Synthesis (SynthDat, Radar, Type1)                  ! "old" approach
      End If

c  Create the time field

      Call Delta_Time (SynthDat, Radar)
  
c  Write the data to the files
  
      Call Wrpln (SynthDat)

 1000 Call Fdate(Real_Time)
      Total = Etime(Elapsed)
      Write (LuMes,'(/"Normal exit",1x,
     #     "of Program WINDSYN at ",a,1x," Elapsed CPU Time(s):",
     #     1x, 3F7.2)') 
     #     Real_Time, Total, Elapsed
      Write (7,'(/"Normal exit"
     > ," of Program WINDSYN at ",a)') Real_Time
      Stop
  
c  Error in reading .prm file so exit.
  
 1020 Write (6,'(/"Error in opening .par, .log, or .err files")')
      Stop

c  An error of some kind reading the parameters from the .prm file

 888  Write (7,'(//"ABORT!! Error:",i5," reading the .prm file"/
     #     "Last Line Read:"/a)') Ierr, Line
      Stop

c  Disc file error exit
  
 999  Write (7,'(/"Abort! Error",i7," with file ",a)')
     #   Ierr, File(1:LenTrim(File))
      Stop

      End Subroutine windsyn
      
c     ***********************************************

      Subroutine Wrpln (SynthDat)  

c     Routine that writes out the data to a disc file
c     Format is:
c     U(i,j,k), V(i,j,k), W(i,j,k), Div(i,j,k), dBZ(i,j,k), Vt(i,j,k),
c     Tdif(i,j,k), Tave(i,j,k)  

      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /Print/ LuMes
  
      Dimension SynthDat (Imax,Jmax,Kmax,Nflds)

      Dimension Out(Imax)

c  Write out the complete U array
  
      Write (LuMes,'(/"Writing U field")')

      Numrec = 0
  
      Do k = 1, Kmax
         Do j = 1, Jmax
            Do i = 1, Imax
               Out(i) = SynthDat(i,j,k,1)
            End Do

            NumRec = NumRec + 1
            Write (3,Rec=NumRec,Err=1001,Iostat=Ierr) (Out(i),i=1,Imax)
         End Do
      End Do
  
c  Write out the complete V array
  
      Write (LuMes,'("Writing V field")')
  
      Do k = 1, Kmax
         Do j = 1, Jmax
            Do i = 1, Imax
               Out(i) = SynthDat(i,j,k,2)
            End Do

            NumRec = NumRec + 1
            Write (3,Rec=NumRec,Err=1001,Iostat=Ierr) (Out(i),i=1,Imax)
         End Do
      End Do
  
c  Write out the complete W array
  
      Write (LuMes,'("Writing W field")')
  
      Do k = 1, Kmax
         Do j = 1, Jmax
            Do i = 1, Imax
               Out(i) = SynthDat(i,j,k,4)
            End Do

            NumRec = NumRec + 1
            Write (3,Rec=NumRec,Err=1001,Iostat=Ierr) (Out(i),i=1,Imax)
         End Do
      End Do
  
c  Write out the divergence array
  
      Write (LuMes,'("Writing divergence field")')
  
      Do k = 1, Kmax
         Do j = 1, Jmax
            Do i = 1, Imax
               Out(i) = SynthDat(i,j,k,5)
            End Do

            NumRec = NumRec + 1
            Write (3,Rec=NumRec,Err=1001,Iostat=Ierr) (Out(i),i=1,Imax)
         End Do
      End Do
  
c  Write out the complete dBZ array
  
      Write (LuMes,'("Writing Z field")')
  
      Do k = 1, Kmax
         Do j = 1, Jmax
            Do i = 1, Imax
               Out(i) = SynthDat(i,j,k,3)
            End Do

            NumRec = NumRec + 1
            Write (3,Rec=NumRec,Err=1001,Iostat=Ierr) (Out(i),i=1,Imax)
         End Do
      End Do
  
c  Write out the terminal velocity field
 
      Write (LuMes,'("Writing terminal velocity field")')
  
      Do k = 1, Kmax
         Do j = 1, Jmax
            Do i = 1, Imax
               Out(i) = SynthDat(i,j,k,9)
            End Do

            NumRec = NumRec + 1
            Write (3,Rec=NumRec,Err=1001,Iostat=Ierr) (Out(i),i=1,Imax)
         End Do
      End Do
  
c  Write out the time difference array
 
      Write (LuMes,'("Writing time difference field")')
  
      Do k = 1, Kmax
         Do j = 1, Jmax
            Do i = 1, Imax
               Out(i) = SynthDat(i,j,k,12)
            End Do

            NumRec = NumRec + 1
            Write (3,Rec=NumRec,Err=1001,Iostat=Ierr) (Out(i),i=1,Imax)
         End Do
      End Do

c  Write out the mean time field
 
      Write (LuMes,'("Writing mean time field")')
  
      Do k = 1, Kmax
         Do j = 1, Jmax
            Do i = 1, Imax
               Out(i) = SynthDat(i,j,k,13)
            End Do

            NumRec = NumRec + 1
            Write (3,Rec=NumRec,Err=1001,Iostat=Ierr) (Out(i),i=1,Imax)
         End Do
      End Do
  
      Write (LuMes,'(/i7," records writing to output file")') NumRec
      Close (3)
  
      Return

 1001 Write (7,'(/"Error: ",I4," when writing output file")') Ierr
      Return

      End Subroutine Wrpln
      
c     ***********************************************

      Subroutine Ldpln (Lu, Radar, n)
 
C  This Subroutine Reads the disc files containing the radar data

c  The contents of the radar data array are as follows:
C          Radar(1,Imax,Jmax,Kmax,Nrdrs) = Azimuth from North [deg]
C          Radar(2,Imax,Jmax,Kmax,Nrdrs) = Radial Velocity [m/s]
C          Radar(3,Imax,Jmax,Kmax,Nrdrs) = Elevation angle [deg]
C          Radar(4,Imax,Jmax,Kmax,Nrdrs) = Reflectivity [dBZ]
c          Radar(5,Imax,Jmax,Kmax,Nrdrs) = Time difference from nominal
c          Radar(6,Imax,Jmax,Kmax,Nrdrs) = Range to bin (m)

      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
  
      Common /Print/ LuMes

      Dimension Radar(6,Imax,Jmax,Kmax,Nrdrs)

      Integer*2 Iang(Imax), Iele(Imax), Ivel(Imax), Iref(Imax),
     #          Irng(Imax), Itme(Imax)

      Do k = 1, Kmax
         Do j = 1, Jmax
            Read (Lu) (Iang(i),Iele(i),Ivel(i),Iref(i),
     #                 Irng(i),Itme(i),i=1,Imax)
  
            Do i = 1, Imax
               Radar(1,i,j,k,n) = Float(Iang(i))/10.0
               Radar(2,i,j,k,n) = Float(Ivel(i))/10.0
               Radar(3,i,j,k,n) = Float(Iele(i))/10.0
               Radar(4,i,j,k,n) = Float(Iref(i))/10.0
               Radar(5,i,j,k,n) = Float(Itme(i))
               Radar(6,i,j,k,n) = Float(Irng(i))/10.0
            End Do
         End Do
      End Do

      Return
      End Subroutine Ldpln
      
c     ***********************************************

      Subroutine ChooZ (SynthDat, Radar)

c  Routine to create the final reflectivity field from the input radars
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /Print/ LuMes

      Dimension SynthDat (Imax,Jmax,Kmax,Nflds)
      Dimension Radar(6,Imax,Jmax,Kmax,Nrdrs)

C   Idbz_parm is a flag which determines how the dBZ will be calculated:
C         -1: Maximum of radars 1 to Nrdrs
C          n: Use radar n
C         -2: Average of radars 1 to Nrdrs
C         -3: Difference radar 1 - radar 2
C       If only a single radar analysis is being done, then
C       Idbz_parm will forced to be 1 (the only radar available)

      Bdbz = 0.0

      Do k=1, Kmax
         Do j=1, Jmax
            Do i=1, Imax
               If (IdBZ_Parm .eq. -1) Then   ! Maximum of all radars
                  Bdbz = Radar(4,i,j,k,1)

                  Do n = 2, Nrdrs
                     If (Radar(4,i,j,k,n) .gt. Bdbz)
     #                    Bdbz = Radar(4,i,j,k,n)
                  End Do

               Else If (IdBZ_Parm .gt. 0) Then ! Use the nth radar
                  Bdbz = Radar(4,i,j,k,IdBZ_Parm)

               Else If (IdBZ_Parm .eq. -2) Then    ! Average of all radars
                  Tot = 0.0
                  Ntot = 0

                  Do n = 1, Nrdrs
                     dbz = Radar(4,i,j,k,n)

                     If (dbz .ne. Baddata) Then
                        Tot = Tot + 10.0**(dbz/10.0)
                        ntot = ntot + 1
                     End If
                  End Do

                  If (Ntot .gt. 0) Then
                     Bdbz = 10.0*Alog10(Tot/Float(Ntot))
                  Else
                     Bdbz = Baddata
                  End If
               Else If (IdBZ_Parm .eq. -3) Then    ! Difference between 1&2
                  dbz1 = Radar(4,i,j,k,1)
                  dbz2 = Radar(4,i,j,k,2)

                  If (dbz1 .ne. Baddata .and. dbz2 .ne. Baddata) Then
                     Bdbz = dbz1 - dbz2
                  Else
                     Bdbz = Baddata
                  End If
               End If

               SynthDat(i,j,k,3)  = Bdbz
            End Do
         End Do
      End Do
  
      Return
      End Subroutine ChooZ
      
c     ***********************************************

      Subroutine Delta_Time (SynthDat, Radar)

c  Routine to create the delta time and mean time fields
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata

      Dimension SynthDat (Imax,Jmax,Kmax,Nflds)
      Dimension Radar(6,Imax,Jmax,Kmax,Nrdrs)

      Do k = 1, Kmax
         Do j = 1, Jmax
            Do i = 1, Imax
               Time_Big = -99999.0

c     Find the largest time

               Do n = 1, Nrdrs
                  If ( Radar(5,i,j,k,n) .ne. Baddata .and.
     #                 Radar(5,i,j,k,n) .gt. Time_Big)
     #                 Time_Big = Radar(5,i,j,k,n)
               End Do

c     Find the smallest time

               Time_Small = 99999.0

               Do n = 1, Nrdrs
                  If ( Radar(5,i,j,k,n) .ne. Baddata .and.
     #                 Radar(5,i,j,k,n) .lt. Time_Small)
     #                 Time_Small = Radar(5,i,j,k,n)
               End Do

c     Calculate the largest time difference

               If (Time_Big .ne. -99999.0 .and.
     #              Time_Small .ne. 99999.0) Then
                  Time_Dif = Abs(Time_Big - Time_Small)
                  If (Nrdrs .eq. 1) Time_Dif = Radar(5,i,j,k,1)
               Else
                  Time_Dif = Baddata
               End If

c     Calculate the mean time from all of the radar beams

               Time_Total = 0.0
               Ntot = 0

               Do n = 1, Nrdrs
                  If (Radar(5,i,j,k,n) .ne. Baddata) Then
                     Time_Total = Time_Total + Radar(5,i,j,k,n)
                     Ntot = Ntot + 1
                  End If
               End Do

               If (Ntot .ne. 0) Then
                  Time_Ave = Time_Total/Ntot
               Else
                  Time_Ave = Baddata
               End If

               SynthDat(i,j,k,12) = Time_Dif
               SynthDat(i,j,k,13) = Time_Ave
            End Do
         End Do
      End Do

      Return
      End Subroutine Delta_Time
      
c     ***********************************************

      Subroutine Comp_Wind (SynthDat, Radar, Type)
 
c  Routine to decide what kind of wind computation is needed
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata

      Common /Print/ LuMes

      Dimension SynthDat (Imax,Jmax,Kmax,Nflds)
      Dimension Radar(6,Imax,Jmax,Kmax,Nrdrs)

      Character*3 Type(Nrdrs)

      Do k = 1, Kmax
         Do j = 1, Jmax
           Do i = 1, Imax
              w = SynthDat(i,j,k,4)
              If (w .eq. Baddata) w = 0.0

c  Two P-3 flight legs (quasi-perpendicular)

              If (Istyle .le. 2 .or. Istyle .eq. 40
     #             .or. Istyle .eq. 12)
     #             Call P3_P3(SynthDat,Radar,i,j,k,w)

c  Two ground-based radar volume scans

              If (Istyle .eq. 20 .or. Istyle .eq. 30)
     #               Call Grnd_Grnd(SynthDat,Radar,i,j,k,w)

c  Ground radar + a P-3 leg

              If (Istyle .eq. 10) Call Grnd_P3(SynthDat,Radar,i,j,k,w)

c  Matejka overdetermined solutions

              If (Istyle .ge. 3 .and. Istyle .le. 9)
     #               Call Over_Det_Sol(SynthDat,Radar,Type,i,j,k,w)

c Two P-3 legs plus a ground based radar volume scan

              If (Abs(Istyle) .eq. 11)
     #             Call Over_Det_Sol(SynthDat,Radar,Type,i,j,k,w)
           End Do
         End Do
      End Do
 
      Return
      End Subroutine Comp_Wind
      
c     ***********************************************

      Subroutine Over_Det_Sol (SynthDat, Radar, Type, i, j, k, w)

c  Routine to compute winds using Matejka's overdetermined technique
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /Debug/ Idebug, Ix1, Ix2, Jy1, Jy2
      Common /Print/ LuMes

      Dimension SynthDat (Imax,Jmax,Kmax,Nflds)
      Dimension Radar(6,Imax,Jmax,Kmax,Nrdrs)

      Character*3 Type(Nrdrs)

      Logical Sub_Vp3, Cmp_Vp3, Use(Nrdrs)

      Dimension Dir_X1(Nrdrs), 
     #          Dir_X2(Nrdrs), 
     #          Dir_X3(Nrdrs), 
     #          Vr(Nrdrs), 
     #          SE_Vr(Nrdrs)

      Call Calc_Vt(i,j,k,SynthDat)
      Term_Vel = SynthDat(i,j,k,9)
      Ref =      SynthDat(i,j,k,3)

c  Get the info from each radar

      Do n = 1, Nrdrs
         Ang = Radar(1,i,j,k,n)
         Vel = Radar(2,i,j,k,n)
         Ele = Radar(3,i,j,k,n)
         Rng = Radar(6,i,j,k,n) * 1000.0
         Se_Vr(n) = 1.0

c  If any of the components are missing return bad data flag
c  The airborne P-3 radial velocities are defined opposite to that of the
c  ground based radars (i.e., negative velocities are receeding from the
c  radar, so need to multiply them by -1 before calling the Matejka
c  routines.
c  Also, since the radial velocity convention is reversed for this routine
c  the terminal fallspeeds must be SUBTRACTED from the cap w values rather
c  than added like in the P-3 routines.

c     Note that the elevation angles are corrected for the effect of the
c     4/3 earth curvature by subtracting the term 1.17E-07 * Rsin(Ele)

         If (Ang .eq. Baddata .or. Ref .eq. Baddata) Then
            Use(n) = .false.
            Vr(n) = Baddata
            Dir_X1(n) = Baddata
            Dir_X2(n) = Baddata
            Dir_X3(n) = Baddata
         Else
            Use(n) = .true.
            If (Istyle .eq. -11) Use(n) = .false.
            If (Type(n) .ne. 'gnd') Vr(n) = -Vel
            r = Rng * Sin(Ele*Ctr)
            Ele_Cor = Ele - (1.17E-07 * r)/Ctr
            Dir_X1(n) = Sin(Ang*Ctr) * Sin(Ele_Cor*Ctr)
            Dir_X2(n) = Cos(Ang*Ctr) * Sin(Ele_Cor*Ctr)
            Dir_X3(n) = Cos(Ele_Cor*Ctr)
         End If
      End Do

c  Print out some debug information?

      If ( i .ge. Ix1 .and. i .le. Ix2 .and.
     #     j .ge. Jy1 .and. j .le. Jy2) Then
         Write (LuMes,'("Grid Point:",3i3," Reflect:",E11.4,
     #        " Term Velocity:",E11.4)') i,j,k, Ref, Term_Vel
         Write (LuMes,'("Rdr  Use    Dir_X1      Dir_X2      ",
     #        "Dir_X3       Vr         SD_Vr")')

         Do n = 1, Nrdrs
            Write (LuMes,'(i2,2x,L1,2x,5E12.5)') n, Use(n),
     #           Dir_X1(n),
     #           Dir_X2(n),
     #           Dir_X3(n),
     #           Vr(n),
     #           SE_Vr(n) 
         End Do
      End If 

c  Perform triple-Doppler solution?

      If (Istyle .eq. 3) Then
         Call Dop_3_Comp_Soln(Nrdrs, Use, 
     #        Dir_X1, Dir_X2, Dir_X3,
     #        Vr, Se_Vr, Baddata, 
     #        Vp1, Se_Vp1, Vp2, Se_Vp2, Vp3, Se_Vp3)

c  Bad horizontal wind data?

         Varen = Sqrt(Se_Vp1*Se_Vp1 + Se_Vp2*Se_Vp2)
         SynthDat(i,j,k,6) = Se_Vp1
         SynthDat(i,j,k,7) = Se_Vp2
         SynthDat(i,j,k,8) = Se_Vp3

         If  (Vp1 .eq. Baddata .or. 
     #        Vp2. eq. Baddata .or. 
     #        Varen .gt. ThresH) Then
            SynthDat(i,j,k,1) = Baddata
            SynthDat(i,j,k,2) = Baddata
            SynthDat(i,j,k,4) = Baddata
            SynthDat(i,j,k,6) = Baddata
            SynthDat(i,j,k,7) = Baddata
         Else
            SynthDat(i,j,k,1) = Vp1 - Su
            SynthDat(i,j,k,2) = Vp2 - Sv
         End If

c  Bad vertical wind data?

         If (Vp3 .ne. Baddata .and. Term_Vel .ne. Baddata) Then
            SynthDat(i,j,k,4) = Vp3 - Term_Vel

            If (Se_Vp3 .gt. ThresV) Then
               SynthDat(i,j,k,4) = Baddata
               SynthDat(i,j,k,8) = Baddata
            End If
         Else
            SynthDat(i,j,k,4) = Baddata
            SynthDat(i,j,k,8) = Baddata
         End If
      End If

c  Overdetermined dual-Doppler solution?

      If (Istyle .eq. 4 .or. Istyle .eq. 9 .or. Abs(Istyle) .eq. 11)
     $     Then 
         Sub_Vp3 = .false.

         If (Term_Vel .ne. Baddata) Then
            Assumed_Vp3 = w + Term_Vel
            Assumed_Se_Vp3 = 1.0
            Call Dop_Soln_Driver_2 (Sub_Vp3, Nrdrs, Use,
     #           Dir_X1, Dir_X2, Dir_X3, Vr, Se_Vr, 
     #           Assumed_Vp3, Assumed_Se_Vp3, Baddata,
     #           Vp1, Se_Vp1, Vp2, Se_Vp2, Vp3, Se_Vp3)

c  Bad data?

            Varen = Sqrt(Se_Vp1*Se_Vp1 + Se_Vp2*Se_Vp2)
            SynthDat(i,j,k,6) = Se_Vp1
            SynthDat(i,j,k,7) = Se_Vp2
            SynthDat(i,j,k,8) = Se_Vp3
            
            If  (Vp1 .eq. Baddata .or. 
     #           Vp2. eq. Baddata .or. 
     #           Varen .gt. ThresH) Then
               SynthDat(i,j,k,1) = Baddata
               SynthDat(i,j,k,2) = Baddata
               SynthDat(i,j,k,6) = Baddata
               SynthDat(i,j,k,7) = Baddata
            Else
               SynthDat(i,j,k,1) = Vp1 - Su
               SynthDat(i,j,k,2) = Vp2 - Sv
            End If

c  No vertical velocity is calculated, i.e., only 2 component solution

            SynthDat(i,j,k,4) = Baddata
            SynthDat(i,j,k,8) = Baddata
         Else
            SynthDat(i,j,k,1) = Baddata
            SynthDat(i,j,k,2) = Baddata
            SynthDat(i,j,k,4) = Baddata
            SynthDat(i,j,k,6) = Baddata
            SynthDat(i,j,k,7) = Baddata
            SynthDat(i,j,k,8) = Baddata
         End If
      End If

c  2 or 3 component solution depending on what is more accurate

      If (Istyle .eq. 5 .or. Istyle .eq. 8) Then
         Sub_Vp3 = .false.

         If (Term_Vel .ne. Baddata) Then
            Assumed_Vp3 = w + Term_Vel
            Assumed_Se_Vp3 = 1.0
            Cmp_Vp3 = .false.

            Call Dop_Soln_Driver_3or2 (
     #           Sub_Vp3, Cmp_Vp3, Nrdrs, Use,
     #           Dir_X1, Dir_X2, Dir_X3, Vr, Se_Vr, 
     #           Assumed_Vp3, Assumed_Se_Vp3, Baddata,
     #           Vp1, Se_Vp1, Vp2, Se_Vp2, Vp3, Se_Vp3)
            
c  Bad horizontal data?

            Varen = Sqrt(Se_Vp1*Se_Vp1 + Se_Vp2*Se_Vp2)
            SynthDat(i,j,k,6) = Se_Vp1
            SynthDat(i,j,k,7) = Se_Vp2
            SynthDat(i,j,k,8) = Se_Vp3

            If ( Vp1 .eq. Baddata .or. 
     #           Vp2. eq. Baddata .or. 
     #           Varen .gt. ThresH) Then
               SynthDat(i,j,k,1) = Baddata
               SynthDat(i,j,k,2) = Baddata
               SynthDat(i,j,k,6) = Baddata
               SynthDat(i,j,k,7) = Baddata
            Else
               SynthDat(i,j,k,1) = Vp1 - Su
               SynthDat(i,j,k,2) = Vp2 - Sv
            End If

c  Bad cap W data?

            If (Vp3 .ne. Baddata) Then
               SynthDat(i,j,k,4) = Vp3 - Term_Vel

               If (Se_Vp3 .gt. ThresV) Then
                  SynthDat(i,j,k,4) = Baddata
                  SynthDat(i,j,k,8) = Baddata
               End If
            Else
               SynthDat(i,j,k,4) = Baddata
               SynthDat(i,j,k,8) = Baddata
            End If
         Else
            SynthDat(i,j,k,1) = Baddata
            SynthDat(i,j,k,2) = Baddata
            SynthDat(i,j,k,4) = Baddata
            SynthDat(i,j,k,6) = Baddata
            SynthDat(i,j,k,7) = Baddata
            SynthDat(i,j,k,8) = Baddata
         End If
      End If

c  3 component solution where possible otherwise 2 comp sol?

      If (Istyle .eq. 6) Then
         Sub_Vp3 = .false.

         If (Term_Vel .ne. Baddata) Then
            Assumed_Vp3 = w + Term_Vel
            Assumed_Se_Vp3 = 1.0

            Call Dop_Soln_Driver_3Else2(Sub_Vp3, Nrdrs, Use,
     #           Dir_X1, Dir_X2, Dir_X3, Vr, Se_Vr, 
     #           Assumed_Vp3, Assumed_Se_Vp3, Baddata,
     #           Vp1, Se_Vp1, Vp2, Se_Vp2, Vp3, Se_Vp3)

c  Bad horizontal data?

            Varen = Sqrt(Se_Vp1*Se_Vp1 + Se_Vp2*Se_Vp2)
            SynthDat(i,j,k,6) = Se_Vp1
            SynthDat(i,j,k,7) = Se_Vp2
            SynthDat(i,j,k,8) = Se_Vp3

            If ( Vp1 .eq. Baddata .or. 
     #           Vp2. eq. Baddata .or. 
     #           Varen .gt. ThresH) Then
               SynthDat(i,j,k,1) = Baddata
               SynthDat(i,j,k,2) = Baddata
               SynthDat(i,j,k,6) = Baddata
               SynthDat(i,j,k,7) = Baddata
            Else
               SynthDat(i,j,k,1) = Vp1 - Su
               SynthDat(i,j,k,2) = Vp2 - Sv
            End If

c  Bad cap W data?

            If (Vp3 .ne. Baddata) Then
               SynthDat(i,j,k,4) = Vp3 - Term_Vel

               If (Se_Vp3 .gt. ThresV) Then
                  SynthDat(i,j,k,4) = Baddata
                  SynthDat(i,j,k,8) = Baddata
               End If
            Else
               SynthDat(i,j,k,4) = Baddata
               SynthDat(i,j,k,8) = Baddata
            End If
         Else
            SynthDat(i,j,k,1) = Baddata
            SynthDat(i,j,k,2) = Baddata
            SynthDat(i,j,k,4) = Baddata
            SynthDat(i,j,k,6) = Baddata
            SynthDat(i,j,k,7) = Baddata
            SynthDat(i,j,k,8) = Baddata
         End If
      End If

c  Variational adjustment technique uses the 2 or 3 but no thresholding

      If (Istyle .eq. 7) Then
         Sub_Vp3 = .False.
         Cmp_Vp3 = .False.

         If (Term_Vel .ne. Baddata) Then
            Assumed_Vp3 = w + Term_Vel
            SE_Wp = SynthDat(i,j,k,8)

            If (SE_Wp .eq. Baddata) Then
               Assumed_Se_Vp3 = 1.0
            Else
               SE_F = SynthDat(i,j,k,11)
               Assumed_Se_Vp3 = Sqrt(SE_Wp*SE_Wp + SE_F*SE_F)
            End If
            
            Call Dop_Soln_Driver_3or2 (
     #           Sub_Vp3, Cmp_Vp3, Nrdrs, Use,
     #           Dir_X1, Dir_X2, Dir_X3, Vr, Se_Vr, 
     #           Assumed_Vp3, Assumed_Se_Vp3, Baddata,
     #           Vp1, Se_Vp1, Vp2, Se_Vp2, Vp3, Se_Vp3)
            
            Varen = Sqrt(Se_Vp1*Se_Vp1 + Se_Vp2*Se_Vp2)

            SynthDat(i,j,k,6) = Se_Vp1         ! SE of u
            SynthDat(i,j,k,7) = Se_Vp2         ! SE of v
            SynthDat(i,j,k,8) = Se_Vp3         ! SE of w

c  Threshold the (u,v) winds if they are bad, otherwise store relative winds

            If ( Vp1 .eq. Baddata .or. 
     #           Vp2. eq. Baddata .or. 
     #           Varen .gt. ThresH) Then
               SynthDat(i,j,k,1) = Baddata
               SynthDat(i,j,k,2) = Baddata
               SynthDat(i,j,k,6) = Baddata
               SynthDat(i,j,k,7) = Baddata
            Else
               SynthDat(i,j,k,1) = Vp1 - Su
               SynthDat(i,j,k,2) = Vp2 - Sv
            End If

c  Store w

            If (Vp3 .ne. Baddata) Then
               SynthDat(i,j,k,4) = Vp3 - Term_Vel ! W
            Else
               SynthDat(i,j,k,4) = Baddata
               SynthDat(i,j,k,8) = Baddata
            End If
         Else
            SynthDat(i,j,k,1) = Baddata
            SynthDat(i,j,k,2) = Baddata
            SynthDat(i,j,k,4) = Baddata
            SynthDat(i,j,k,6) = Baddata
            SynthDat(i,j,k,7) = Baddata
            SynthDat(i,j,k,8) = Baddata
         End If
      End If

c  Print out some debug information?

      If ( i .ge. Ix1 .and. i .le. Ix2 .and.
     #     j .ge. Jy1 .and. j .le. Jy2) Then
         Write (LuMes,'(/
     #        "      U          V          Z    ",
     #        "      W        Div        SD U ")')
         Write (LuMes,'(6E11.4)') (SynthDat(i,j,k,n),n=1,6)
         Write (LuMes,'(
     #        "    SD V       SD W        Vt    ",
     #        "   Div SD      Vt SD        Wp")')
         Write (LuMes,'(6E11.4)') (SynthDat(i,j,k,n),n=7,11), Vp3
         Write (LuMes,'(//)')
      End If

      Return
      End Subroutine Over_Det_Sol
      
c     ***********************************************

      Subroutine P3_P3 (SynthDat, Radar, i, j, k, w)

c  Routine to compute winds from two P-3 Doppler radar estimates
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /Print/ LuMes

      Dimension SynthDat (Imax,Jmax,Kmax,Nflds)
      Dimension Radar(6,Imax,Jmax,Kmax,Nrdrs)

      Ang1   = Radar(1,i,j,k,1)
      Ang2   = Radar(1,i,j,k,2)
      Velt1  = Radar(2,i,j,k,1)
      Velt2  = Radar(2,i,j,k,2)
      Ele1   = Radar(3,i,j,k,1)
      Ele2   = Radar(3,i,j,k,2)

      If (((Istyle .le. 2 .or. Istyle .eq. 12) .and.
     #     (Ang1 .lt. 0.0 .or. Ang2 .lt. 0.0)) .or. 
     #     (Istyle .eq. 40 .and. Ang1 .lt. 0.0))
     $     Then
         SynthDat(i,j,k,1) = Baddata
         SynthDat(i,j,k,6) = Baddata
         SynthDat(i,j,k,2) = Baddata
         SynthDat(i,j,k,7) = Baddata
         Return
      End If

      Cos_Ele1 = 0.0
      Sin_Ele1 = 0.0
      Cos_Ele2 = 0.0
      Sin_Ele2 = 0.0
      Cos_Ang1 = 0.0
      Sin_Ang1 = 0.0
      Cos_Ang2 = 0.0
      Sin_Ang2 = 0.0
  
      If (Istyle .le. 2 .or. Istyle .eq. 12) Then  ! Dual-Doppler

         Fact = 1.0 / Sin((Ang2 - Ang1)*Ctr)

         Call Calc_Vt (i,j,k,SynthDat)
         Term_Vel = SynthDat(i,j,k,9)

         Cos_Ele1 = Cos(Ele1*Ctr)
         Sin_Ele1 = Sin(Ele1*Ctr)
         Cos_Ele2 = Cos(Ele2*Ctr)
         Sin_Ele2 = Sin(Ele2*Ctr)
         Cos_Ang1 = Cos(Ang1*Ctr)
         Sin_Ang1 = Sin(Ang1*Ctr)
         Cos_Ang2 = Cos(Ang2*Ctr)
         Sin_Ang2 = Sin(Ang2*Ctr)
 
         Term1 = ((Velt1 + (w+Term_Vel)*Cos_Ele1) * Cos_Ang2)/Sin_Ele1
         Term2 = ((Velt2 + (w+Term_Vel)*Cos_Ele2) * Cos_Ang1)/Sin_Ele2
 
         u = Fact * (Term1 - Term2)
         ru = u - Su
 
         Term3 = ((Velt1 + (w+Term_Vel)*Cos_Ele1) * Sin_Ang2)/Sin_Ele1
         Term4 = ((Velt2 + (w+Term_Vel)*Cos_Ele2) * Sin_Ang1)/Sin_Ele2
 
         v = -Fact * (Term3 - Term4)
         rv = v - Sv

      Else                           ! Single-Doppler 

         Fact = 1.0 / Sin((Ang2 - Ang1)*Ctr)

         Call Calc_Vt(i,j,k,SynthDat)
         Term_Vel = SynthDat(i,j,k,9)

         CorRadVel = Velt1+(w+Term_Vel)*Abs(Cos(Ele1*Ctr))
  
         If (Ele1 .ne. 0.0 .and. Ele1 .ne. 180.0) Then
            Ws = Abs(CorRadVel/Sin(Ele1*Ctr))
  
            If (CorRadVel .le. 0.0) Then
               If (Ang1 .ge. 180.0) Then
                  Wd = Ang1 - 180.0
               Else
                  Wd = Ang1 + 180.0
               End If
            Else
               Wd = Ang1
            End If

            Ru = Ucomp(Ws,Wd)
            Rv = Vcomp(Ws,Wd)
         Else       ! No horizontal wind can be calculated at zenith/nadir
            Ru = Baddata
            Rv = Baddata
         End If
      End If

c  Store the computed relative wind components

      SynthDat(i,j,k,1) = Ru
      SynthDat(i,j,k,2) = Rv

c  Store the computed standard errors of the horizontal winds computed
c  from the geometry of the viewing angles

      Var = 1.0/(Sin((Ang2 - Ang1)*Ctr) * Sin((Ang2 - Ang1)*Ctr)) *
     #     (1.0/(Sin_Ele1*Sin_Ele1) + 1.0/(Sin_Ele2*Sin_Ele2))
      SE = Sqrt(Var)
      SynthDat(i,j,k,6) = SE
      SynthDat(i,j,k,7) = SE

c  Get rid of the bad data

      If (SE .gt. ThresH) Then
         SynthDat(i,j,k,1) = Baddata
         SynthDat(i,j,k,2) = Baddata
         SynthDat(i,j,k,6) = Baddata
         SynthDat(i,j,k,7) = Baddata
      End If

      Return
      End Subroutine P3_P3
      
c     ***********************************************

      Subroutine Grnd_Grnd (SynthDat, Radar, i, j, k, w)
  
c  Subroutine to compute winds from the two radial velocities
c  of the ground-based Doppler radars
c  Dual-Doppler synthesis will not be done if the beam crossing
c  angle is < 20 deg
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata

      Dimension SynthDat (Imax,Jmax,Kmax,Nflds)  
      Dimension Radar(6,Imax,Jmax,Kmax,Nrdrs)

      Azm1 = Radar(1,i,j,k,1)
      Vel1 = Radar(2,i,j,k,1)
      Ele1 = Radar(3,i,j,k,1)
      Azm2 = Radar(1,i,j,k,2)
      Vel2 = Radar(2,i,j,k,2)
      Ele2 = Radar(3,i,j,k,2)
  
c  Both radar estimates ok?
  
      If ((Istyle .eq. 20 .and. (Azm1 .lt. 0.0 .or. Azm2 .lt. 0.0))
     1    .or. (Istyle .eq. 30 .and. Azm1 .lt. 0.0)) Then
         SynthDat(i,j,k,1) = Baddata
         SynthDat(i,j,k,2) = Baddata
         SynthDat(i,j,k,6) = Baddata
         SynthDat(i,j,k,7) = Baddata
         Return
      End If

      Call Calc_Vt(i,j,k,SynthDat)  
      Term_Vel = SynthDat(i,j,k,9)

      Cos_Ele1 = 0.0
      Sin_Ele1 = 0.0
      Cos_Ele2 = 0.0
      Sin_Ele2 = 0.0
  
      If (Istyle .eq. 20) Then  ! Dual-Doppler
  
c  Beam crossing angle > 20 degrees?  If not, skip it
  
         If (Abs(Sin((Azm1-Azm2)*Ctr)) .lt. 0.342) Then
            SynthDat(i,j,k,1) = Baddata
            SynthDat(i,j,k,2) = Baddata
            SynthDat(i,j,k,6) = Baddata
            SynthDat(i,j,k,7) = Baddata
            Return
         End If

         If (i .eq. 56 .and. j .eq. 55) Then
            write (6,*) 'Azm1, Azm2', i,j,k, Azm1, Azm2, Ele1, Ele2, Sin((Azm1-Azm2)*Ctr)
         End If

         Cos_Ele1 = Cos(Ele1*Ctr)
         Sin_Ele1 = Sin(Ele1*Ctr)
         Cos_Ele2 = Cos(Ele2*Ctr)
         Sin_Ele2 = Sin(Ele2*Ctr)

         Cos_Azm1 = Cos(Azm1*Ctr)
         Cos_Azm2 = Cos(Azm2*Ctr)
         Sin_Azm1 = Sin(Azm1*Ctr)
         Sin_Azm2 = Sin(Azm2*Ctr)

         TermA = ((Vel2 - (w+Term_Vel)*Sin_Ele2) * Cos_Azm1)/Cos_Ele2
         TermB = ((Vel1 - (w+Term_Vel)*Sin_Ele1) * Cos_Azm2)/Cos_Ele1

c  Compute u component relative to the moving grid

         u = (1.0 / Sin((Azm2-Azm1)*Ctr)) * (TermA - TermB)
         u = u - Su

         TermA = ((Vel2 - (w+Term_Vel)*Sin_Ele2) * Sin_Azm1)/Cos_Ele2
         TermB = ((Vel1 - (w+Term_Vel)*Sin_Ele1) * Sin_Azm2)/Cos_Ele1

c  Compute v component relative to the moving grid

         v = (1.0 / Sin((Azm1-Azm2)*Ctr)) * (TermA - TermB)
         v = v - Sv
      Else  ! Single-Doppler
         CorRadVel = (Vel1-(w+Term_Vel)*Sin(Ele1*Ctr))
         Ws = Abs(CorRadVel)/Cos(Ele1*Ctr)
 
         If (CorRadVel .gt. 0.) Then
            If (Azm1 .ge. 180.) Then
               Wd = Azm1 - 180.
            Else
               Wd = Azm1 + 180.
            End If
         Else
            Wd = Azm1
         End If

         u = Ucomp(Ws,Wd)
         v = Vcomp(Ws,Wd)
      End If
 
      SynthDat(i,j,k,1) = u
      SynthDat(i,j,k,2) = v

      Var = 1.0/(Sin((Azm2 - Azm1)*Ctr) * Sin((Azm2 - Azm1)*Ctr)) *
     #     (1.0/(Cos_Ele1*Cos_Ele1) + 1.0/(Cos_Ele2*Cos_Ele2))
      SE = Sqrt(Var)

      SynthDat(i,j,k,6) = SE
      SynthDat(i,j,k,7) = SE

      Return
      End Subroutine Grnd_Grnd
      
c     ***********************************************

      Subroutine Grnd_P3 (SynthDat, Radar, i, j, k, w)
 
c  Subroutine to compute winds from the two radial velocities
c  from the P-3 leg and the ground-based radar
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata

      Dimension SynthDat (Imax,Jmax,Kmax,Nflds)  
      Dimension Radar(6,Imax,Jmax,Kmax,Nrdrs)

 
      Azm_P3 = Radar(1,i,j,k,1)
      Vel_P3 = Radar(2,i,j,k,1)
      Ele_P3 = Radar(3,i,j,k,1)
      Azm_Gr = Radar(1,i,j,k,2)
      Vel_Gr = Radar(2,i,j,k,2)
      Ele_Gr = Radar(3,i,j,k,2)
  
      If (Azm_P3 .lt. 0.0 .or. Azm_Gr .lt. 0.0) Then
         SynthDat(i,j,k,1) = Baddata
         SynthDat(i,j,k,2) = Baddata
         SynthDat(i,j,k,6) = Baddata
         SynthDat(i,j,k,7) = Baddata
         Return
      End If

      Call Calc_Vt(i,j,k,SynthDat)  
      Term_Vel = SynthDat(i,j,k,9)
  
c  Beam crossing angle > 20 degrees?  If not, skip it
  
      If (Abs(Sin((Azm_Gr-Azm_P3)*Ctr)) .lt. 0.342) Then
         SynthDat(i,j,k,1) = Baddata
         SynthDat(i,j,k,2) = Baddata
         SynthDat(i,j,k,6) = Baddata
         SynthDat(i,j,k,7) = Baddata
         Return
      End If
  
      Fact = 1.0 / Sin((Azm_Gr - Azm_P3)*Ctr)

      Cos_Gr_El = Cos(Ele_Gr*Ctr)
      Sin_Gr_El = Sin(Ele_Gr*Ctr)
      Sin_P3_El = Abs(Sin(Ele_P3*Ctr))
      Cos_P3_El = Cos(Ele_P3*Ctr)
      Cos_Gr_Az = Cos(Azm_Gr*Ctr)
      Sin_Gr_Az = Sin(Azm_Gr*Ctr)
      Sin_P3_Az = Sin(Azm_P3*Ctr)
      Cos_P3_Az = Cos(Azm_P3*Ctr)
 
      Term1 = ((Vel_Gr - (w+Term_Vel)*Sin_Gr_El) * Cos_P3_Az)/Cos_Gr_El
      Term2 = ((Vel_P3 + (w+Term_Vel)*Cos_P3_El) * Cos_Gr_Az)/Sin_P3_El
 
      u = Fact * (Term1 + Term2)
      u = u - Su
 
      Term3 = ((Vel_Gr - (w+Term_Vel)*Sin_Gr_El) * Sin_P3_Az)/Cos_Gr_El
      Term4 = ((Vel_P3 + (w+Term_Vel)*Cos_P3_El) * Sin_Gr_Az)/Sin_P3_El
 
      v = -Fact * (Term3 + Term4)
      v = v - Sv
 
      SynthDat(i,j,k,1) = u
      SynthDat(i,j,k,2) = v

      Var = 1.0/(Sin((Azm_Gr-Azm_P3)*Ctr) * Sin((Azm_Gr-Azm_P3)*Ctr)) *
     #     (1.0/(Cos_Gr_El*Cos_Gr_El) + 1.0/(Sin_P3_El*Sin_P3_El))
      SE = Sqrt(Var)
      SynthDat(i,j,k,6) = SE
      SynthDat(i,j,k,7) = SE

      Return
      End Subroutine Grnd_P3
      
c     ***********************************************

      Subroutine Calc_Divergence (SynthDat)
  
C   This Subroutine computes horizontal divergence using
C   a horizontal grid spacing of SX and SY (km)
c   Filtering (smoothing) is done at the same time in order to calculate
c   the correct standard errors of the filtered divergence
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
  
      Common /Print/ LuMes

      Dimension SynthDat (Imax,Jmax,Kmax,Nflds)

      Dimension Wt_X (-100:100), Wt_Y (-100:100)

      Delx = Sx * 1000.0           ! x grid spacing in m
      Dely = Sy * 1000.0           ! y grid spacing in m
      Max_Reach = 100

c  Set up the weights for the binomial and Leise filters

      If (IsmTyp .eq. 1) Then          ! the "Leise" filter
         Call Multi_Step_Leise_Filter_Weights(Nsmth, Max_Reach,
     #        Wt_X, Irx)
         Call Multi_Step_Leise_Filter_Weights(Nsmth, Max_Reach,
     #        Wt_Y, Iry)           
      Else If (IsmTyp .eq. 2) Then     ! the "binominal" filter
         Call Multi_Step_Binomial_Filter_Weights(Nsmth, Max_Reach,
     #        Wt_X, Irx)
         Call Multi_Step_Binomial_Filter_Weights(Nsmth, Max_Reach,
     #        Wt_Y, Iry)           
      Else
         Return
      End If

c  Compute the smoothed divergence  

      Do k = 1, Kmax
         Call Divergence_Filter_2D(
     #        SynthDat(1,1,k,1), SynthDat(1,1,k,6),
     #        SynthDat(1,1,k,2), SynthDat(1,1,k,7),
     #        Imax, Jmax, Imax, Jmax,
     #        Delx, Dely, Baddata,
     #        Wt_X, Max_Reach, Irx,
     #        Wt_Y, Iry, Imax, Jmax,
     #        SynthDat(1,1,k,5), SynthDat(1,1,k,10))
      End Do
  
      Return
      End Subroutine Calc_Divergence
      
c     ***********************************************

      Subroutine Synthesis_New (SynthDat, Radar, Type)
  
C  This Subroutine computes vertical motion (w) using the Matjeka
c  variational adjustment technique.

      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /w_Techniques/ Iw_at_top, I_w_0
      Common /Print/ LuMes
      Common /Debug/ Idebug, Ix1, Ix2, Jy1, Jy2
      Common /Sounding/ Z1, Tv1, Z2, Tv2, Z3, Rho3
      Common /Limits/ Zsfc, Iters, Kup, Klim, Top_Hgt

      Dimension SynthDat (Imax,Jmax,Kmax,Nflds)
      Dimension Radar(6,Imax,Jmax,Kmax,Nrdrs)
      Character*3 Type(Nrdrs)

      Logical Success, Use_Top_w

c  Iterate the calculation Iters times to converge w

      Do Loop = 1, Iters
         Write (LuMes,'(//"ITERATION Number:",i3,
     #        " Calculating Winds and Reflectivity")') Loop

c  Calculate the reflectivity

         Call ChooZ (SynthDat, Radar)

c  Create the wind vectors from the radial velocities

         Call Comp_Wind(SynthDat, Radar, Type)

c  Fill small holes in the horizonatal planes

         If (Ihole_Fill .eq. 1) Then
            Call Hole_Fill(SynthDat, Loop)
         End If

c  Fill small gaps in the vertical columns

         Nu_Gaps = 0
         Nv_Gaps = 0
         Nw_Gaps = 0
         Nz_Gaps = 0
         Nt_Gaps = 0

         If (Klim .gt. 0) Then
            Call Vfill (SynthDat,1,Nu_Gaps) ! u
            Call VFill (SynthDat,6,Nu_Gaps) ! SD u
            Call Vfill (SynthDat,2,Nv_Gaps) ! v
            Call Vfill (SynthDat,7,Nv_Gaps) ! SD v
            Call VFill (SynthDat,3,Nz_Gaps) ! dBz
            
            If (Istyle .ge. 5 .and. Istyle .le. 7) Then
               Call Vfill (SynthDat,4,Nw_Gaps) ! w
               Call VFill (SynthDat,8,Nw_Gaps) ! SD w

               If (Istyle .eq. 7) Then
                  Call VFill (SynthDat,9,Nt_Gaps)  ! Vt
                  Call VFill (SynthDat,11,Nt_Gaps) ! SD Vt
               End If
            End If
         End If

c  Calculate the horizontal divergence with a normal centered
c  difference operator with smoothing
      
         Call Calc_Divergence (SynthDat)
         
c  Interpolate over small gaps in the divergence

         If (Klim .gt. 0) Then
            Call VFill (SynthDat,5,Ndiv_Gaps) ! DIV
         Else
            Ndiv_Gaps = 0
         End If

c  Smooth the u,v,w,Vt fields (if desired)

         If (Nsmth .gt. 0 .and. Ismtyp .gt. 0)
     #        Call Smooth_It (SynthDat)
         
c  Initialize the statistics variables

         Ngood = 0
         Nbad = 0
         Nfailed = 0
         Nduds = 0
         Nowp = 0
         Nhy = 0
         Nw0 = 0
         Nobn = 0

c  Is this a "dual" solution
c-------------------------------------------------------------------------
         If (Istyle .eq. 4) Then         ! dual
c-------------------------------------------------------------------------
            Do j = 1, Jmax
               Do i = 1, Imax
                  Use_Top_w = .False. ! w=0 at the top boundary

c O'Brien adjustment required?

                  If (Iobr .eq. 1) Then
                     Call w_calc_OBrien (i,j, SynthDat, Success, Ireason,
     #                    Use_Top_w)
                     
                     If (Success) Then
                        Ngood = Ngood + 1
                        Cycle ! Jump to the end of the loop
                     Else
                        If (Ireason .eq. 1) Nduds = Nduds + 1
                        If (Ireason .eq. 2) Nbad = Nbad + 1
                        If (Ireason .eq. 3) Then   ! try no O'Brien
                           Nfailed = Nfailed + 1

                           Call w_Calc_Unadjusted (i,j, SynthDat, Success,
     #                          Ireason, Use_Top_w)

                           If (Success) Then
                              Nobn = Nobn + 1
                              Cycle
                           End If
                        End If
                     End If

c     A complete failure to find acceptable w's in the column? Erase the w's

                     If (.not. Success) Then
                        Do k = 1, Kmax
                           SynthDat(i,j,k,4) = Baddata
                           SynthDat(i,j,k,8) = Baddata
                        End Do
                     End If

c  No O'Brien adjustment

                  Else
                     Call w_Calc_Unadjusted (i,j, SynthDat, Success,
     #                    Ireason, Use_Top_w)
                     
                     If (Success) Then
                        Ngood = Ngood + 1
                        Cycle ! finished with this column
                     Else
                        If (Ireason .eq. 1) Nduds = Nduds + 1
                        If (Ireason .eq. 2) Nbad = Nbad + 1
                        If (Ireason .eq. 3) Nfailed = Nfailed + 1
                     End If

c     A complete failure to find acceptable w's in the column? Erase the w's

                     If (.not. Success) Then
                        Do k = 1, Kmax
                           SynthDat(i,j,k,4) = Baddata
                           SynthDat(i,j,k,8) = Baddata
                        End Do
                     End If
                  End If
               End Do
            End Do
c-------------------------------------------------------------------------
         Else If (Istyle .eq. 5 .or. Istyle .eq. 6) Then         ! Hybrid
c-------------------------------------------------------------------------
            Do j = 1, Jmax
               Do i = 1, Imax
                  Use_Top_w = .True. ! w=wp-f at the top boundary

c O'Brien adjustment required?

                  If (Iobr .eq. 1) Then
                     Call w_Calc_OBrien (i,j, SynthDat, Success, Ireason,
     #                    Use_Top_w)
                     
                     If (Success) Then
                        Ngood = Ngood + 1
                        Cycle ! finished with this column
                     Else
                        If (Ireason .eq. 1) Nduds = Nduds + 1       ! dud column
                        If (Ireason .eq. 2) Nbad = Nbad + 1         ! no good wp's
                        If (Ireason .eq. 3) Nfailed = Nfailed + 1   ! failure

c  Try the w=0 at the top condition?

                        If (I_w_0 .eq. 1) Then
                           Use_Top_w = .False. ! w=0 at the top boundary

                           If (Iobr .eq. 1)
     #                          Call W_Calc_OBrien (i,j, SynthDat,
     #                          Success,
     #                          Ireason, Use_Top_w)
                           If (Iobr .eq. 0)
     #                          Call w_Calc_Unadjusted (i,j, SynthDat,
     #                          Success,
     #                          Ireason, Use_Top_w)

                           If (Success) Then
                              Nw0 = Nw0 + 1
                              Cycle ! finished with this column
                           End If
                        End If
                     End If

c     A complete failure to find acceptable w's in the column? Erase the w's

                     If (.not. Success) Then
                        Do k = 1, Kmax
                           SynthDat(i,j,k,4) = Baddata
                           SynthDat(i,j,k,8) = Baddata
                        End Do
                     End If

c  No O'Brien adjustment

                  Else
                     Call w_Calc_Unadjusted (i,j, SynthDat, Success,
     #                    Ireason, Use_Top_w)
                     
                     If (Success) Then
                        Ngood = Ngood + 1
                        Cycle ! finished with this column
                     Else
                        If (Ireason .eq. 1) Nduds = Nduds + 1
                        If (Ireason .eq. 2) Nbad = Nbad + 1
                        If (Ireason .eq. 3) Nfailed = Nfailed + 1
                     End If

c     A complete failure to find acceptable w's in the column? Erase the w's

                     If (.not. Success) Then
                        Do k = 1, Kmax
                           SynthDat(i,j,k,4) = Baddata
                           SynthDat(i,j,k,8) = Baddata
                        End Do
                     End If
                  End If
               End Do
            End Do
c-------------------------------------------------------------------------
         Else If (Istyle .eq. 7) Then   ! The variational adjustment approach
c-------------------------------------------------------------------------
            Do j = 1, Jmax
               Do i = 1, Imax

                  Call Var_Adjust_Col (i,j, SynthDat, Success, Ireason)
                  
                  If (Success) Then
                     Ngood = Ngood + 1
                     Cycle             ! finished with this column
                  Else
                     If (Ireason .eq. 1) Nduds = Nduds + 1
                     If (Ireason .eq. 2) Nowp = Nowp + 1
                     If (Ireason .eq. 3) Nfailed = Nfailed + 1

c  For an inaccurate wp at the top, try w=0 as a boundary condition

                     If (Ireason .eq. 2 .and. I_w_0 .eq. 1) Then
                        Use_Top_w = .False. ! w=0 at the top boundary
                        If (Iobr .eq. 1)
     #                       Call W_Calc_OBrien (i,j, SynthDat, Success,
     #                       Ireason, Use_Top_w)
                        If (Iobr .eq. 0)
     #                       Call w_Calc_Unadjusted (i,j, SynthDat, Success,
     #                       Ireason, Use_Top_w)

                        If (Success) Then
                           Nw0 = Nw0 + 1
                           Cycle ! finished with this column
                        End If
                     End If

c     Try the Hybrid approach on the columns that failed

                     If (Ireason .eq. 3) Then
                        If (Iw_at_Top .eq. 1) Then
                           Use_Top_w = .True. ! w=0 at the top boundary
                           If (Iobr .eq. 1)
     #                          Call w_Calc_OBrien (i,j, SynthDat,
     #                          Success,
     #                          Ireason, Use_Top_w)
                           If (Iobr .eq. 0)
     #                          Call w_Calc_Unadjusted (i,j, SynthDat,
     #                          Success,
     #                          Ireason, Use_Top_w)
                           
                           If (Success) Then
                              Nhy = Nhy + 1
                              Cycle ! finished with this column
                           End If
                        End If

c  Try the w=0 at the top condition?

                        If (I_w_0 .eq. 1) Then
                           Use_Top_w = .False. ! w=0 at the top boundary
                           If (Iobr .eq. 1)
     #                          Call W_Calc_OBrien (i,j, SynthDat,
     #                          Success,
     #                          Ireason, Use_Top_w)
                           If (Iobr .eq. 0)
     #                          Call w_Calc_Unadjusted (i,j, SynthDat,
     #                          Success,
     #                          Ireason, Use_Top_w)

                           If (Success) Then
                              Nw0 = Nw0 + 1
                              Cycle ! finished with this column
                           End If
                        End If
                     End If

c     A complete failure to find acceptable w's in the column? Erase the w's
                     
                     If (.not. Success) Then
                        Do k = 1, Kmax
                           SynthDat(i,j,k,4) = Baddata
                           SynthDat(i,j,k,8) = Baddata
                        End Do
                     End If
                  End If
               End Do           ! j
            End Do            ! i
         End If               ! Style of synthesis
      End Do                  ! Loop

c  Print out the column statistics

      Ncolumns = Imax * Jmax
c-------------------------------------------------------------------------
      If (Istyle .eq. 7) Then
c-------------------------------------------------------------------------
         Nfailed = Nfailed + Nowp
         
         Write (LuMes,100) Ncolumns, Nduds, Ngood, Nfailed, 
     #        Nowp
         Write (LuMes,101) Nfailed, Nhy, Nw0
         Write (LuMes,102) Ndiv_Gaps, Nu_Gaps, Nv_Gaps,
     #        Nw_Gaps, Nz_Gaps, Nt_Gaps
         
 100     Format (/'Total number of columns:',
     #        i5/'# of columns of no data:'
     #        i5/'# of columns where variational adjustment worked:'
     #        i5/'# of columns where variational adjustment failed:'
     #        i5/'# of columns that had Wp_SD larger than thresh:'
     #        i5)

 101     Format(/'Of the',i5,' columns that failed,',i5,
     #        ' succeeded with the hybrid approach'/32x,i5
     #        ' succeeded with w=0 as a top boundary condition')

 102     Format(/'# of vertical column DIV gaps interpolated over:'
     #        i5/'# of vertical column u gaps interpolated over:'
     #        i5/'# of vertical column v gaps interpolated over:'
     #        i5/'# of vertical column w gaps interpolated over:'
     #        i5/'# of vertical column dBZ gaps interpolated over:'
     #        i5/'# of vertical column Vt gaps interpolated over:'
     #        i5)
c-------------------------------------------------------------------------
      Else If (Istyle .ge. 5 .and. Istyle .le. 6) Then
c-------------------------------------------------------------------------
         Ntot_failed = Nfailed + Nbad + Nduds
         Write (LuMes,200) Ncolumns, Nduds, Ngood,
     #        Nbad, Nfailed
         Write (LuMes,201) Ntot_failed, Nw0
         Write (LuMes,202) Ndiv_Gaps, Nu_Gaps, Nv_Gaps,
     #        Nw_Gaps, Nz_Gaps

 200     Format (/'# of columns:'
     #        i5/'# of columns of no data:'
     #        i5/'# of columns where hybrid approach worked:'
     #        i5/'# of columns that did not have a good w at top:'
     #        i5/'# of columns where vertical integration failed:'
     #        i5)

 201     Format(/'Of the',i5,' columns that failed,',i5,
     #        ' succeeded with w=0 at cloud top')

 202     Format (/'# of vertical column DIV gaps interpolated over:'
     #        i5/'# of vertical column u gaps interpolated over:'
     #        i5/'# of vertical column v gaps interpolated over:'
     #        i5/'# of vertical column w gaps interpolated over:'
     #        i5/'# of vertical column dBZ gaps interpolated over:'
     #        i5)
c-------------------------------------------------------------------------
      Else     ! Istyle = 4 should be the only one left
c-------------------------------------------------------------------------
         Write (LuMes,300) Ncolumns, Nduds, Ngood, Nfailed
         If (Iobr .eq. 1) Write (LuMes,301) Nfailed, Nobn
         Write (LuMes,302) Ndiv_Gaps, Nu_Gaps, Nv_Gaps, Nz_Gaps

 300     Format (/'# of columns:'
     #        i5/'# of columns of no data:'
     #        i5/'# of columns where vertical integration worked:'
     #        i5/'# of columns where vertical integration failed:'
     #        i5)

 301     Format(/'Of the',i5,' columns that failed,',i5,
     #        ' succeeded with no OBrien adjustment')

 302     Format (/'# of vertical column DIV gaps interpolated over:'
     #        i5/'# of vertical column u gaps interpolated over:'
     #        i5/'# of vertical column v gaps interpolated over:'
     #        i5/'# of vertical column dBZ gaps interpolated over:'
     #        i5)
      End If

c  Lastly, adjust the (u,v) fields for the adjusted divergences

      Call Adj_Winds_from_Div (SynthDat)

      Return
      End Subroutine Synthesis_New
      
c     ***********************************************

      Subroutine Var_Adjust_Col (i,j, SynthDat, Success, Ireason)

c  Subroutine to perform variational adjustment for the calculation
c  of vertical air velocity in a vertical column specified by (i,j)
c  Ireason is a flag to tell the calling program why the adjustment
c  failed.
c   Ireason = 1   No divergence or wp's in the column
c   Ireason = 2   No wp's in the column lower than the threshold
c   Ireason = 3   Success .false. returned by subroutine

      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata

      Common /Print/ LuMes
      Common /Debug/ Idebug, Ix1, Ix2, Jy1, Jy2
      Common /Sounding/ Z1, Tv1, Z2, Tv2, Z3, Rho3
      Common /Limits/ Zsfc, Iters, Kup, Klim, Top_Hgt

      Dimension SynthDat (Imax,Jmax,Kmax,Nflds)

      Dimension 
     #     Z(0:Kmax),           !  heights of levels in column
     #     W_Cnstrnt(0:Kmax),   !  strong w constraint
     #     W_Cnstrnt_SD(0:Kmax),!  strong w constraint standard dev

     #     Div(0:Kmax),         !  divergence
     #     Div_Adj(0:Kmax),     !  "adjusted" horizontal divergence
     #     Div_SD(0:Kmax),      !  standard error of hor divergence
     #     Div_Adj_SD(0:Kmax),  !  Stand error of "adjusted" divergence

     #     Wp(0:Kmax),          !  vertical hydrometeor velocity (cap w)
     #     Wp_Adj(0:Kmax),      !  "adjusted" vert hydromet. velocity
     #     Wp_SD(0:Kmax),       !  standard error of wp
     #     Wp_Adj_SD(0:Kmax),   !  Stand error of "adjusted" wp

     #     F(0:Kmax),           !  terminal fall speed
     #     F_Adj(0:Kmax),       !  "adjusted" terminal fall speed
     #     F_SD(0:Kmax),        !  standard error of terminal fall speed
     #     F_Adj_SD(0:Kmax),    !  Stand error of "adjusted" fallspeeds

     #     W(0:Kmax),           !  vertical air velocity 
     #     W_SD(0:Kmax)         !  Stand error of "adjusted" w

      Logical Success

c  Define the columns to perform the w calculation upon
c  An extra layer is added for the surface boundary conditions

      Ndiv = 0
      Nwp = 0
      SDWP_Min = 999.0
      Ireason = 0

c  Initialize the arrays

      Do k = 0, Kmax
         Wp(k)           = Baddata
         Wp_SD(k)        = Baddata
         Div(k)          = Baddata
         Div_SD(k)       = Baddata
         F(k)            = Baddata
         F_SD(k)         = Baddata
         W_Cnstrnt(k)    = Baddata
         W_Cnstrnt_SD(k) = Baddata
         Div_Adj(k)      = Baddata
         Div_Adj_SD(k)   = Baddata
         Wp_Adj(k)       = Baddata
         Wp_Adj_SD(k)    = Baddata
         F_Adj(k)        = Baddata
         F_Adj_SD(k)     = Baddata
         W(k)            = Baddata
         W_SD(k)         = Baddata
      End do

c  Set up the column data

      Do k = 1, Kmax
         Z(k)      = Hgt_of_Level(Z0,k,Sz) * 1000.0
         Div(k)    = SynthDat(i,j,k,5)
         Div_SD(k) = SynthDat(i,j,k,10)
         If (Div(k) .ne. Baddata) Ndiv = Ndiv + 1
         Vt        = SynthDat(i,j,k,9)
         wv        = SynthDat(i,j,k,4)
         SD_Vt     = SynthDat(i,j,k,11)

         If (Vt .ne. Baddata .and. SD_Vt .ne. Baddata) Then
            F(k)    = -Vt                ! the minus sign is required for
            F_SD(k) = SD_Vt              ! the Matejka routines
         End If
         
         If (Vt .ne. Baddata .and. wv .ne. Baddata) Then
            Wp(k)    = wv + Vt
            Wp_SD(k) = SynthDat(i,j,k,8)
         End If

c  Find the minimum standard error of Wp in the column

         If (Wp(k) .ne. Baddata) Then
            Nwp = Nwp + 1
            SD_w = Wp_SD(k)
            
            If (F_SD(k) .ne. Baddata) Then
               SD_w = Wp_SD(k)*Wp_SD(k) + F_SD(k)*F_SD(k)
               SD_w = Sqrt(SD_w)
            End If
            
            SDWP_Min = Amin1(SDWP_Min,SD_w)
         End If
      End Do                    ! end column definition
      
c  Find out if this is a dud column (no div, no wp's)
            
      If (Ndiv .eq. 0 .and. Nwp .eq. 0) Then
         Success = .False.
         Ireason = 1
         Return
      End IF
               
c  Are the Wp's in the column too inaccrate to use?

      If ( SDWP_Min .gt. ThresV .and. 
     #     SDWP_Min .ne. 999.0) Then
         Success = .False.
         Ireason = 2
         Return
      End If

c  First level missing? Then bring down upper levels up to Kup

      If (Div(1) .eq. Baddata .and. Kup .gt. 1) Then
         Do k = 2, Kup
            If (Div(k) .ne. Baddata) Then
               Do kk = k-1, 1, -1
                  Div(kk)    = Div(k)
                  Div_SD(kk) = Div_SD(k)
               End Do
               
               Go To 1
            End If
         End Do
      End If

c  Bogus in the surface values for a boundary condition

 1    Z(0)      = Zsfc
      Div(0)    = Div(1)
      Div_SD(0) = Div_SD(1)
      Wp(0)     = Wp(1)
      Wp_SD(0)  = Wp_SD(1)
      F(0)      = F(1)
      F_SD(0)   = F_SD(1)
      Div_Adj_SD(0) = Baddata
      Wp_Adj_SD(0) = Baddata
      F_Adj_Sd(0) = Baddata
      W_SD(0) = Baddata
      
c     Perform the equivalent of an O'Brien correction at the bottom?

      W_Cnstrnt(0)    = Baddata
      W_Cnstrnt_SD(k) = Baddata
      
      If (Iobr .eq. 1) Then
         W_Cnstrnt(0)    = 0.0
         W_Cnstrnt_SD(0) = 0.0
      End If

c     Perform the variational adjustment column by column

      Call Complete_Var_Col_Soln (
     #     Z1, Tv1, Z2, Tv2, Z3, Rho3, Kmax,
     #     Z, Div, Div_SD, Wp, Wp_SD, F, F_SD, W_Cnstrnt,
     #     Baddata, 
     #     Div_Adj, Wp_Adj, F_Adj, W, Success)
               
c     Print out debug column information?

      If ( i .ge. Ix1 .and. i .le. Ix2 .and.
     #     j .ge. Jy1 .and. j .le. Jy2) Then
         Write (LuMes,100) i,j
 100     Format (//2i3/
     #        ' n    Height     Div     Div_SD       WP      ',
     #        'Wp_SD       F       F_SD    W_constr')
         
         Do k = 0, Kmax
            Write (LuMes,'(i2,8e10.3)') k, Z(k), Div(k), Div_SD(k),
     #           Wp(k), Wp_SD(k), F(k), F_SD(k), W_Cnstrnt(k)
         End Do
         
         Write (LuMes,'(/2i3," Success:",L3)') i,j,Success
         Write (LuMes,'(" n   Height    Div_Adj    WP_ADJ    ",
     #        "F_ADJ      W")')

         Do k = 0, Kmax
            Write (LuMes,'(i2,5e10.3)') k, Z(k), Div_Adj(k), Wp_Adj(k),
     #           F_Adj(k), W(k)
         End Do
      End If

c  Was it successful?

      If (.not. Success) Then
         Ireason = 3
         Return
      End If

c  Calculate the standard deviation of the div, F, and w

      Call Complete_Var_Col_Soln_SD (
     #     Z1, Tv1, Z2, Tv2, Z3, Rho3, Kmax,
     #     Z, Div, Div_SD, Wp, Wp_SD, F, F_SD,
     #     W_Cnstrnt, W_Cnstrnt_SD, Baddata,
     #     Div_Adj_SD, Wp_Adj_SD, F_Adj_SD, W_SD)
      
c     Print out debug column information?

      If ( i .ge. Ix1 .and. i .le. Ix2 .and.
     #     j .ge. Jy1 .and. j .le. Jy2) Then
         Write (LuMes,104) i,j
 104     Format (//2i3,'         Adjusted standard errors'/
     #        ' n       Height     Div_Adj_SD   Wp_Adj_SD',
     #        '    F_Adj_SD        W_SD')
         
         Do k = 0, Kmax
            Write (LuMes,'(i2,1x,5E13.4)')
     #           k, Z(k), Div_Adj_SD(k),
     #           Wp_Adj_SD(k), F_Adj_Sd(k), W_SD(k)
         End Do
      End If
      
c     Store the adjusted variables
      
      Do k = 1, Kmax
         SynthDat(i,j,k,4) = W(k)
         SynthDat(i,j,k,5) = Div_Adj(k)
         SynthDat(i,j,k,9) = -F_Adj(k)
         
         If (W_SD(k) .gt. ThresV)  SynthDat(i,j,k,4) = Baddata
      End Do
         
      Return
      End Subroutine Var_Adjust_Col
      
c     ***********************************************

      Subroutine w_Calc_OBrien (i,j, SynthDat,
     #     Success, Ireason, Use_Top_w)

c  Routine to perform the "hybrid" vertical velocity solution where
c  w is defined at the column top by the Wp-Vt derived from the overdetermined
c  triple equation solution and the rest of the column's w's are derived
c  by continuity integration subject to a w=0 lower boundry (O'Brien-type
c  correction.

c  If this O'Brien solution does not work (usually because of no divergence
c  estimate at the lower boundary), then a straight downward integration is
c  used.

c  Ireason is a parameter to explain why the vertical integration failed:
c     Ireason = 1  No div estimates in the column
c     Ireason = 2  No good wp's exist at cloud top (valid for hybrid only)
c     Ireason = 3  Vertical integration failed (usually because the divergence
c                  column did not extend to the ground)

      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /Print/ LuMes
      Common /Limits/ Zsfc, Iters, Kup, Klim, Top_Hgt
      Common /Sounding/ Z1, Tv1, Z2, Tv2, Z3, Rho3
      Common /Debug/ Idebug, Ix1, Ix2, Jy1, Jy2

      Dimension SynthDat (Imax,Jmax,Kmax,Nflds)

      Dimension 
     #     Z(0:Kmax),           !  heights of levels in column
     #     Div(0:Kmax),         !  divergence
     #     SE_Div(0:Kmax),      !  standard error of hor divergence
     #     Div_Adj(0:Kmax),     !  "adjusted" horizontal divergence
     #     SE_Div_Adj(0:Kmax),  !  standard error of adjusted div
     #     SE_W(0:Kmax),        !  standard error of vertical velocity
     #     W(0:Kmax),           !  vertical air velocity 
     #     W_Adj(0:Kmax),       !  adj. w
     #     SE_W_Adj(0:Kmax)     !  stand error of adj w

      Logical Success, Use_Top_w

      W_BC1 = 0.0         ! lower boundary condition
      SE_W_BC1 = 0.0      ! guess at the standard deviation of lower w
      I_BC1 = 0           ! index of the lower boundary condition
      Ireason = 0         ! assume its going to work

c  Initialize all columns to initial values

 4    Do k = 0, Kmax
         Z(k)         = Baddata
         Div(k)       = Baddata               ! divergence column
         SE_Div(k)    = Baddata               ! standard deviation of div
         Div_Adj(k)   = Baddata               ! "adjusted" div
         SE_Div_Adj(k)= Baddata               ! "adjusted" div standard dev
         W(k)         = Baddata               ! vertical air motion
         SE_W(k)      = Baddata               ! standard deviation of w
         W_Adj(k)     = Baddata               ! "adjusted" w
         SE_W_Adj(k)  = Baddata               ! standard deviation of adj w
      End Do

c  Define the heights

      Do k = 1, Kmax
         Z(k) = Hgt_of_Level(Z0,k,Sz) * 1000.0
      End Do

      Z(0) = Zsfc

c  Define the columns to perform the w calculation upon
c  An extra layer is added for the surface boundary conditions
c  Set up the column data

      Do k = 1, Kmax
         W(k)         = SynthDat(i,j,k,4)           ! vertical air velocity (w)
         SE_W(k)      = SynthDat(i,j,k,8)           ! standard dev of w
         Div(k)       = SynthDat(i,j,k,5)           ! divergence
         SE_Div(k)    = SynthDat(i,j,k,10)          ! standard dev of div
      End Do                                    ! end column definition

c  Find out if this is a dud column
c  Search from the top down to find the top of the column

      Ktop = -1

      If (Use_Top_w) Then                       ! must find a good div and w
         Do k = Kmax, 1, -1                     ! to be a valid cloud top
            If ( W(k)   .ne. Baddata .and. 
     #           Div(k) .ne. Baddata) Then
               Ktop = k
               Go To 5
            End If
         End Do
      
         Success = .False.
         Ireason = 1
         Return
      Else                                      ! look only for a good div
         Do k = Kmax, 1, -1
            If (Div(k) .ne. Baddata) Then
               Ktop = k
               Go To 5
            End If
         End Do
         
         Success = .False.
         Ireason = 1
         Return
      End If

c  Valid column?

 5    If (Hgt_of_Level(Z0, Ktop, Sz) .lt. Top_Hgt .and. Use_Top_w) Then
         Ireason = 2
         Success = .False.
         Return
      End If

c  First level missing? Then bring down upper levels (up to Kup)

      If (Div(1) .eq. Baddata .and. Kup .gt. 1) Then
         Do k = 2, Kup
            If (Div(k) .ne. Baddata) Then
               Do kk = k-1, 1, -1
                  Div(kk)    = Div(k)
                  SE_Div(kk) = SE_Div(k)
               End Do
               
               Go To 6
            End If
         End Do
      End If

c  Bogus in the surface values for a boundary condition

 6    Div(0)     = Div(1)
      SE_Div(0)  = SE_Div(1)
      
c  Define the boundary conditions at the column top

      If (Use_Top_w) Then
         W_BC2    = W(Ktop)           ! vertical air velocity
         SE_W_BC2 = SE_W(Ktop)        ! standard dev of w

c  Average the top 2 levels to get a more representative w

         If (W(Ktop-1) .ne. Baddata) Then
            W_BC2    = (   W(Ktop) +    W(Ktop-1))/2.0
            SE_W_BC2 = (SE_W(Ktop) + SE_W(Ktop-1))/2.0
         End If

      Else
         W_BC2    = 0.0
         SE_W_BC2 = 0.0
      End If

c  If no good w exists at cloud top, then abandon the column

      If (W_BC2 .eq. Baddata .or. SE_W_BC2 .gt. ThresV) Then
         Ireason = 2
         Success = .False.
         Return
      End If
      
c  Print out diagnostic column information if requested

      If ( i .ge. Ix1 .and. i .le. Ix2 .and.
     #     j .ge. Jy1 .and. j .le. Jy2) Then
         Write (LuMes,'(/"i,j I_BC1, IBC_2, W_BC2, SE_W_BC2:"
     #        2i3,1x,2i3,1x,2E11.4)') i,j, I_BC1, Ktop,
     #        W_BC2, SE_W_BC2
         write (LuMes,'(/" n    Height     Div     SE_DIV      W",
     #        "        SE_W")')
         
         Do k = 0, Kmax
            write (LuMes,'(i2,5e10.3)') k, Z(k), 
     #           Div(k), SE_DIV(k), W(k), SE_W(k)
         End Do
      End If

c  Perform the O'Brien adjustment column by column

      Call W_Adjust_Levels (
     #     Z1, Tv1, Z2, Tv2, Z3, Rho3, Kmax,
     #     Z, Div, SE_Div, I_BC1, W_BC1, SE_W_BC1,
     #     Ktop, W_BC2, SE_W_BC2, Baddata,
     #     Div_Adj, SE_Div_Adj, W_Adj, SE_W_Adj, Success)
      
c     Was it successful?  Keep track of how many good and bad columns
      
      If (.not. Success) Then
         Ireason = 3
         Return
      End If

c     Store the adjusted variables if the w standard error is less than
c     the prescribed threshold
         
      Do k = 1, Kmax
         If ( SE_W_Adj(k) .lt. ThresV .and.
     #        SE_W_Adj(k) .ne. Baddata) Then
            SynthDat(i,j,k,4) = W_Adj(k)
            SynthDat(i,j,k,5) = Div_Adj(k)
            SynthDat(i,j,k,8) = SE_W_Adj(k)
            SynthDat(i,j,k,10)= SE_Div_Adj(k)
         Else
            SynthDat(i,j,k,4) = Baddata
            SynthDat(i,j,k,5) = Baddata
            SynthDat(i,j,k,8) = Baddata
            SynthDat(i,j,k,10)= Baddata
         End If
      End Do
      
c     Diagnostic printout of the columns?
      
      If ( i .ge. Ix1 .and. i .le. Ix2 .and.
     #     j .ge. Jy1 .and. j .le. Jy2) Then
         Write (LuMes,'("Success:",L1/,
     #        /" n     Height    Div_Adj   SE_DIV_ADJ",
     #        "   W_Adj     SE_W_Adj")') Success
         Do k = 0, Kmax
            write (LuMes,'(i3,5E11.4)') k, Z(k), Div_Adj(k),
     #           SE_Div_Adj(k), W_Adj(k), SE_W_Adj(k)
         End Do
      End If
      
c  In case there are still good divergences above the level where w defined
c  the column top, use a straight upward integration starting at
c  column top with the good w

      If (Ktop .lt. Kmax .and. Use_Top_w) Then
         Call W_Upward_Downward_Levels (
     #        Z1, Tv1, Z2, Tv2, Z3, Rho3, Kmax,
     #        Z, Div, SE_Div, Ktop, W_BC2, SE_W_BC2, 
     #        Baddata, W_Adj, SE_W_Adj)
         
c     Store the w's above the column top

         Do k = Ktop+1, Kmax
            If ( SE_W_Adj(k) .lt. ThresV .and.
     #           SE_W_Adj(k) .ne. Baddata) Then
               SynthDat(i,j,k,4) = W_Adj(k)
               SynthDat(i,j,k,8) = SE_W_Adj(k)
               SynthDat(i,j,k,5) = Div(k)
               SynthDat(i,j,k,10)= SE_Div(k)
            Else
               SynthDat(i,j,k,4) = Baddata
               SynthDat(i,j,k,8) = Baddata
               SynthDat(i,j,k,5) = Baddata
               SynthDat(i,j,k,10)= Baddata
            End If
         End Do
      
c     Diagnostic printout of the columns?
      
         If ( i .ge. Ix1 .and. i .le. Ix2 .and.
     #        j .ge. Jy1 .and. j .le. Jy2) Then
            Write (LuMes,'(/"Column after upward integration to top"/,
     #           /" n     Height    Div_Adj   SE_DIV_ADJ",
     #           "   W_Adj     SE_W_Adj")') 
            Do k = 0, Kmax
               write (LuMes,'(i3,5E11.4)') k, Z(k), SynthDat(i,j,k,5),
     #               SynthDat(i,j,k,10), SynthDat(i,j,k,4), SynthDat(i,j,k,8)
            End Do
         End If
      End If

      Return
      End Subroutine w_Calc_OBrien
      
c     ***********************************************

      Subroutine w_Calc_Unadjusted (i,j, SynthDat, Success, Ireason,
     #     Use_Top_w)
  
C  This Subroutine performs upward or downward vertical
c  integration of horizontal divergence to generate estimates
c  of vertical air velocity.
c  The boundary condition is the top w from the
c  overdetermined dual or triple equation solution technique (for
c  the hybrid approach) or w=0 (for the dual approach).

c  Ireason is a parameter to explain why the vertical integration failed:
c     Ireason = 1  No div estimates in the column
c     Ireason = 2  No good wp's exist at cloud top (valid for hybrid only)
c     Ireason = 3  Vertical integration failed (usually because the divergence
c                  column did not extend to the ground)
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /Print/ LuMes
      Common /Sounding/ Z1, Tv1, Z2, Tv2, Z3, Rho3
      Common /Limits/ Zsfc, Iters, Kup, Klim, Top_Hgt
      Common /Debug/ Idebug, Ix1, Ix2, Jy1, Jy2

      Dimension SynthDat (Imax,Jmax,Kmax,Nflds)

      Logical Success, Use_Top_w

      Dimension 
     #     Z(0:Kmax),           !  heights of levels in column
     #     Div(0:Kmax),         !  divergence
     #     SE_Div(0:Kmax),      !  standard error of hor divergence
     #     SE_W(0:Kmax),        !  standard error of vertical velocity
     #     W(0:Kmax)            !  vertical air velocity 

c  Initialize all columns to initial values

 4    Do k = 0, Kmax
         Div(k)       = Baddata               ! divergence column
         SE_Div(k)    = Baddata               ! standard deviation of div
         W(k)         = Baddata               ! vertical air motion
         SE_W(k)      = Baddata               ! standard deviation of w
      End Do

c  Define the heights

      Do k = 1, Kmax
         Z(k) = Hgt_of_Level(Z0,k,Sz) * 1000.0
      End Do

      Z(0) = Zsfc
      Ireason = 0

c  Define the columns to perform the w calculation upon
c  An extra layer is added for the surface boundary conditions
c  Set up the column data

      Do k = 1, Kmax
         W(k)         = SynthDat(i,j,k,4)
         SE_W(k)      = SynthDat(i,j,k,8)
         Div(k)       = SynthDat(i,j,k,5)
         SE_Div(k)    = SynthDat(i,j,k,10)
      End Do

c  Find out if this is a dud column
c  Search from the top down to find the first good level

      If (Use_Top_w) Then
         Do k = Kmax, 1, -1
            If ( W(k)   .ne. Baddata .and. 
     #           Div(k) .ne. Baddata) Then
               Ktop = k
               Go To 5
            End If
         End Do
      
         Success = .False.
         Ireason = 1
         Return
      Else
         Do k = Kmax, 1, -1
            If (Div(k) .ne. Baddata) Then
               Ktop = k
               Go To 5
            End If
         End Do
         
         Success = .False.
         Ireason = 1
         Return
      End If
               
c  Valid column?

 5    If (Hgt_of_Level(Z0, Ktop, Sz) .lt. Top_Hgt .and. Use_Top_w) Then
         Ireason = 2
         Success = .False.
         Return
      End If

c  Bogus in the surface values for a boundary condition

      Div(0)     = Div(1)
      SE_Div(0)  = SE_Div(1)

      If (Idir_Int .eq. 1) Then ! Upward integration
         I_BC = 0
         W_BC = 0.0
         SE_W_BC = 0.0          ! Downward hybrid integration
      Else                      ! starting from Z=Ktop with
         If (Use_Top_w) Then
            W_BC    = SynthDat(i,j,Ktop,4) ! boundary conditions at top
            SE_W_BC = SynthDat(i,j,Ktop,8)
            I_BC = Ktop
         Else
            W_BC = 0.0
            SE_W_BC = 0.0
            I_BC = Ktop
         End If

         If (W_BC .eq. Baddata .or. SE_W_BC .gt. ThresV) Then
            W_BC = 0.0
            SE_W_BC = 1.0
            Ireason = 2
            Success = .False.
            Return
         End If
      End If

c     Perform the variational adjustment column by column

      Call W_Upward_Downward_Levels (
     #     Z1, Tv1, Z2, Tv2, Z3, Rho3, Kmax,
     #     Z, Div, SE_Div, I_BC, W_BC, SE_W_BC, 
     #     Baddata, W, SE_W)

      Success = .True.

c  Store the adjusted variables if the standard deviation is acceptable
               
      Do k = 1, Kmax
         If (SE_W(k) .lt. ThresV .and. SE_W(k)
     #        .ne. Baddata) Then
            SynthDat(i,j,k,4) = W(k)
            SynthDat(i,j,k,8) = SE_W(k)
         Else
            SynthDat(i,j,k,4) = Baddata
            SynthDat(i,j,k,8) = Baddata
         End If
      End Do
      
c     Print out columns if debug info is requested

      If ( i .ge. Ix1 .and. i .le. Ix2 .and.
     #     j .ge. Jy1 .and. j .le. Jy2) Then
         Write (LuMes,100) i,j, Itop
 100     Format (//2i3,' Top Index:',i3//
     #        ' n    Height      Div       SE_DIV')
         
         Do k = 0, Kmax
            Write (LuMes,'(i2,3E11.4)') k, Z(k), Div(k), SE_DIV(k)
         End Do
         
         Write (LuMes,'(/" I_BC, W_BC, SE_W_BC:",
     #        i4,2E11.4)') I_BC, W_BC, SE_W_BC
         Write (LuMes,'(" n    Height       W         SE_W")')
         
         Do k = 0, Kmax
            Write (LuMes,'(I2,3E11.4)') k, Z(k), W(k), SE_W(k)
         End Do
      End If

      Return
      End Subroutine w_Calc_Unadjusted
      
c     ***********************************************

      Subroutine Synthesis (SynthDat, Radar, Type)
  
C  This Subroutine computes vertical motion (w) using horizontal
C  divergence derived from dual-Doppler measurements of the
C  horizontal wind through application of the equation of anelastic
C  continuity.
C  Integration can be either from the top to the bottom or
C  from the bottom to the top of any layer comprising two or more
C  levels lying within the vertical grid.  Upon completion of the
C  integration, an O'Brien-type adjustment may be applied to force
C  satisfaction of boundary conditions by distributing the diagnosed
C  divergence "error" (associated with deviation of w from 0 at the Last
C  level of integration) uniformly throughout the column.
C  The column is assumed to encompass the layer extending from
C  one-half grid level (i.e. a distance of Sz/2) above the top-most
C  level at which divergence is known to below the bottom-most non-missing
C  level.  If the bottom-most level is at k=1, the column is extended
C  to the surface with divergence assumed to be constant in the
C  layer from z=0 to z=Z0, and the lower b.c. is applied at z=0; other-
C  wise the column is assumed to extend one-half grid level below the
C  bottom-most level and the lower b.c. is applied there.
C  Written 9/87 by Brad Smull.
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /Vert/ IFirst, Last, Szm, Z0m
      Common /Print/ LuMes
      Common /Limits/ Zsfc, Main_loop_Max, Kup, Klim, Top_Hgt

      Dimension SynthDat(Imax,Jmax,Kmax,Nflds)
      Dimension Radar(6,Imax,Jmax,Kmax,Nrdrs)

      Character*3 Type(Nrdrs)

      Real*8 Alpha, Beta, Kappa, Dens
  
C  Statement function to define the mean density variation with
C  altitude (km).  First constant is mean surface density (kg/m**3)
C  value across analysis domain.  Second constant is scale height.
 
      Dens(z) = 1.1904 / Exp(z/9.58125)
  
c   Perform O''Brien''s Adjustment?

 15   If (Iobr .eq. 1) Then
         Write (LuMes,'(/"O''Brien Adjustment will be performed")')
      Else
         Write (LuMes,'(/"O''Brien Adjustment will NOT be performed")')
         Iobr = 0
      End If
  
      Write (LuMes,'(/"The scheme will be iterated ",i2," times")')
     1         Main_Loop_Max
  
      Szm = Sz * 1000.0
      Z0m = Z0 * 1000.0
  
      Do Main_Loop = 1, Main_Loop_Max
         Write (LuMes,'(/"Vertical velocity scheme: Iteration ",i2,
     #                   " of ",i2)')  Main_Loop, Main_Loop_Max

c  Calculate the reflectivity
         
         Call ChooZ (SynthDat, Radar)
  
c  Create the wind vectors from the radial velocities

         Call Comp_Wind (SynthDat, Radar, Type)

c  Fill small holes in the horizonatal planes

         If (Ihole_Fill .eq. 1) Then
            Call Hole_Fill (SynthDat, Main_Loop)
         End If

c  Fill small gaps in the vertical columns

         If (Klim .gt. 0) Then
            Call Vfill (SynthDat,1,Nu_Gaps) ! u
            Call VFill (SynthDat,6,Nu_Gaps) ! SE u
            Call Vfill (SynthDat,2,Nv_Gaps) ! v
            Call Vfill (SynthDat,7,Nv_Gaps) ! SE v
            Call Vfill (SynthDat,4,Nw_Gaps) ! w
            Call VFill (SynthDat,8,Nw_Gaps) ! SE w
            Call VFill (SynthDat,3,Nz_Gaps) ! dBz
         Else
            Nu_Gaps = 0
            Nv_Gaps = 0
            Nw_Gaps = 0
            Nz_Gaps = 0
         End If

c  Calculate the divergence with a normal centered-difference operator

         Call Calc_Divergence (SynthDat)

c  Interpolate over small gaps in the divergence

         If (Klim .gt. 0) Then
            Call VFill (SynthDat,5,Ndiv_Gaps) ! DIV
         Else
            Ndiv_Gaps = 0
         End If

c  Smooth the u,v,w,Vt fields if desired

         If (Nsmth .gt. 0 .and. Ismtyp .gt. 0)
     #        Call Smooth_It (SynthDat)

         Iflagupr = 0
         Iflaglwr = 0
  
         Do j = 1, Jmax
            Do i = 1, Imax

c  Wl is the b.c. to be applied at the appropriate level adjacent to
c  k=First
c  If "quad" solution, then use first good point from the top
  
               Wl = 0.0
               Wtop = 0.0

c Initialize lower & upper column index limits assuming all Data missing

               If (Idir_Int .eq. +1) Then ! Upward integration
                  Last = -99
                  IFirst = 99
                  Kstrt = 1
                  Kstop = Kmax
               Else             ! Downward integration
                  Last = 99
                  IFirst = -99
                  Kstrt = Kmax
                  Kstop = 1
               End If
 
c  Vertical integration to obtain w

               Do K = Kstrt, Kstop, Idir_Int
                  If (SynthDat(i,j,k,5) .eq. Baddata) Then ! missing data
 
c  Allow no integration across data gaps (missing divergence) in column;
c  if good divergence data have already been found, terminate the column
c  for purposes of vertical velocity calculation.

                     If (((Idir_Int .eq. +1) .and. (Last .ne. -99)) .or.
     #               ((Idir_Int .eq. -1) .and. (Last .ne. 99))) Then

                        Do kk = K, Kstop, Idir_Int
                           SynthDat(i,j,kk,4) = Baddata
                        End Do

                        Go To 50
                     End If
 
c  If no good data has yet been found in this column, set w(k)=missing
c  and advance to next level if not at bottom level.

                     SynthDat(i,j,k,4) = Baddata
                     Cycle
                  End If
 
c  Since a good data point has been found, re-define column index limits
c  to be used in Subroutine W_Calc and later by O'Brien scheme.
 
                  If (Idir_Int .eq. +1) Then ! Upward integration
                     Last = Max0(k,Last)
                     IFirst = Min0(k,IFirst)
                  Else          ! Downward integration
                     If ((Istyle .eq. 8) .and. IFirst .eq. -99
     #                    .and. SynthDat(i,j,k,4) .ne. Baddata) Then
                        Wtop = SynthDat(i,j,k,4)
                        Wl = Wtop
c                        Write (7,'(" Wtop:",f6.1," at i,j,k:",3i3)')
c     #                       Wtop, i,j,k
                     End If

                     Last = Min0(k,Last)
                     IFirst = Max0(k,IFirst)
                  End If
 
                  Call W_Calc (SynthDat,i,j,k,wl,w)

c  Save vertical velocity at current level for use at next level.

                  Wl = w
                  SynthDat(i,j,k,4) = w

                  If (Abs(w) .gt. 40.0 .and. W .ne. Baddata)
     #                 Write (7,'(3i4," W>40 (before OBrien):",
     #                 E11.4," wtop:",E11.4," First:",i2," Last:",
     #                 i2)') i,j,k, w, wtop, IFirst, Last
               End Do           ! K
  
c  Diagnostics for b.c. applied at level k=First
  
 50            If (Idir_Int .eq. +1) Then
                  If ((Kstrt .eq. 1) .and. (IFirst .ne. 1) .and.
     #                 (IFirst .ne. 99) .and. (IFirst .ne. Last)) Then
c                     dBZ = SynthDat(i,j,IFirst,3)
c                     Write (7,'("lower b.c. applied at (",2(i2,","),
c     #                    i2,") where dBZ= ",f4.0)')i,j,IFirst,dBZ
                     Iflaglwr = Iflaglwr + 1
                  End If
               Else             ! Idir_Int = -1
                  If ((IFirst .ne. -99) .and. (IFirst .ne. Last)) Then
                     dBZ =  SynthDat(i,j,IFirst,3)
 
                     If (dBZ .gt. 20.) Then
c                        Write (7,'(" UPPER b.c. applied at (",
c     #                       2(i2,","),i2,") where dBZ= ",f4.0)')
c     #                       i,j,IFirst,dBZ
                        Iflagupr = Iflagupr + 1
                     End If
                  End If
               End If
  
               If (Iobr .eq. 1) Then ! invoke O'Brien scheme
  
c  Adjust w using technique described in Sec. 2 of O'Brien (1970) and
c  Ray et al. (1980), i.e. assuming divergence error is constant
c  with height.  Re-calculate iteratively until  b.c. at appropriate
c  level adjacent to k=Last is satisfied.
c  Perform adjustment starting from First non-missing level, and
c  only if there is more than one non-missing level.
  
                  If (((Idir_Int .eq. +1) .and. (IFirst .eq.  99)) .or.
     #                 ((Idir_Int .eq. -1) .and. (IFirst .eq. -99)) .or.
     #                 (IFirst .eq. Last)) Go to 70 ! no data to adjust
  
cbfs Only adjust columns extending below 1.5 km
  
                  If (Idir_Int .eq. -1 .and. Last .gt. 2) Go to 70
  
c  Diagnostics for b.c. applied at level k=Last via O'Brien scheme
  
                  If (Idir_Int .eq. +1) Then
                     dBZ =  SynthDat(i,j,Last,3)
  
                     If (dBZ .gt. 20.) Then
c                        Write (7,'(" UPPER b.c. applied at (",
c     #                       2(i2,","),i2,") where dBZ= ",f4.0)')
c     #                       i,j,Last,dBZ
                        Iflagupr = Iflagupr + 1
                     End If
                  Else          ! Idir_Int = -1
                     If ((Kstop .eq. 1) .and. (Last. ne. 1)) Then
                        dBZ = SynthDat(i,j,Last,3)
c                        Write (7,'(" lower b.c. applied at (",2(i2,
c     #                       ","),i2,") where dBZ= ",f4.0)')
c     #                       i,j,Last,dBZ
                        Iflaglwr = Iflaglwr + 1
                     End If
                  End If
  
c  Calculate depth of column over which divergence error accrued so that
c  adjustment may be made.
c  Extend column upward one-half grid level above topmost data point.
  
                  dk = Abs(Last - IFirst) + 1.5
  
                  If (Last .eq. 1 .or. IFirst .eq. 1) Then
                     Depth = dk*Szm + Z0m
                  Else
                     Depth = dk*Szm
                  End If
 
c  Wt is b.c. to be applied at appropriate level adjacent to k=Last
 
                  wt = 0.0
                                             ! or sfc if Last = 1
c  Progressively under-relax to obtain convergence; if diagnostic write
c  statements indicate some points have not converged, First try
c  increasing the upper limit of IRepeat, then if necessary decrease
c  the step (third Do-loop parameter) used in the iteration of Relax.
c  Currently choice of parameters limits total number of iterations to
c  100 per column, but typically convergence requires < 10.
 
                  Do iRelax = 1, 10   !2.5, 7.0, 0.5
                     Relax = 2.0 + Float(iRelax)*0.5

                     Do IRepeat = 1, 10
                        If (Idir_Int .eq. -1) Then ! downward integration
                           If (Last .eq. 1) Then
 
c  Impose lower b.c. at surface, assuming constant divergence in layer
c  between surface and Z0.
 
                              Kappa = (Dens(0.) - Dens(Z0))/
     #                             Dens(Z0/2.0)/2.0
                              Beta  = 1.0 - Kappa
                              Alpha = 1.0 + Kappa
                              Wz = SynthDat(i,j,1,4) * Beta/Alpha +
     1                             SynthDat(i,j,1,5) * Z0m/Alpha
                              wdif = wz - wt
                              If (Abs(wdif) .lt. 0.1) Go to 70
                              Adj_Div = -wdif/depth/Relax
                           Else
 
c     Impose b.c. one-half grid level below bottom-most data point.
 
                              Kappa = (Dens(Z0 + (Float(Last)-1.5)*Sz) -
     1                             Dens(Z0 + (Last-1)*Sz)) /
     2                             Dens(Z0 + (Float(Last)-1.25)*Sz)
                              Beta  = 1.0 - Kappa
                              Alpha = 1.0 + Kappa
                              Wz = SynthDat(i,j,Last,4) * Beta/Alpha +
     1                             SynthDat(i,j,Last,5) * Szm/2.0/Alpha
                              wdif = wz - wt
                              If (Abs(wdif) .lt. 0.1) Go to 70
                              Adj_Div = -wdif/depth/Relax
                           End If
                        Else    ! upward integration
  
c  Impose upper b.c. one-half grid level above top-most data point
  
                           Kappa = (Dens(Z0+(Last-1)*Sz) -
     #                          Dens(Z0+(Float(Last)-0.50)*Sz)) /
     #                          Dens(Z0+(Float(Last)-0.25)*Sz) / 2.0
                           Beta  = 1.0 + Kappa
                           Alpha = 1.0 - Kappa
                           Wz = SynthDat(i,j,Last,4) * Beta/Alpha -
     #                          SynthDat(i,j,Last,5) * Szm/2.0/Alpha
                           wdif = wz - wt
                           If (Abs(wdif) .lt. 0.1) Go to 70
                           Adj_Div = wdif/depth/Relax
                        End If
  
                        Wl = 0.0 ! assumed b.c. at k=First

                        If (Istyle .eq. 8 .and. Idir_Int .eq. -1) 
     #                       Wl = Wtop
  
                        Do k = IFirst, Last, Idir_Int ! adjust div & recalc w
  
C  *** assume constant divergence error with height ***
  
                           SynthDat(i,j,k,5) = SynthDat(i,j,k,5) + Adj_Div
                           Call W_Calc(SynthDat,i,j,k,Wl,Adj_w)
                           Wl = Adj_w
                           SynthDat(i,j,k,4) = Adj_W
                        End Do
                     End Do     ! IRepeat
                  End Do        ! Relax
  
                  Write (7,'(" O''Brien failed to converge"
     #                 " at boundary"
     #                 " near (",2(i2,","),i2,") where w= ",
     #                 e12.5)')i,j,Last,Wz
               End If           ! End of O'Brien scheme

 70            Continue
            End Do              ! I
         End Do                 ! J
         
         Write (LuMes,'("# of points flagged due to questionable ",
     1        "upper b.c.: ",i4,/," Number of points flagged due to ",
     2        "questionable lower b.c.: ",i4,/)') Iflagupr, Iflaglwr
         Write (LuMes,'("Number of single-level divergence gaps",1x,
     #        "interpolated over:",i6)') Ndiv_gaps
         Write (LuMes,'("Number of single-level u gaps",1x,
     #        "interpolated over:",i6)') Nu_gaps
         Write (LuMes,'("Number of single-level v gaps",1x,
     #        "interpolated over:",i6)') Nv_gaps
      End Do                                 !  Main_loop

c  Lastly, adjust the (u,v) fields for the adjusted divergences

      Call Adj_Winds_from_Div (SynthDat)

      Return
      End Subroutine Synthesis
      
c     ***********************************************
  
      Subroutine W_Calc (SynthDat, i, j, k, wl, w)
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /Vert/ IFirst, Last, Szm, Z0m
  
      Dimension SynthDat(Imax,Jmax,Kmax,Nflds)
  
      Real*8 Alpha, Beta, Kappa, Dens
  
C   Statement function to define the mean density variation with
C   altitude (km).  First constant is mean surface density (kg/m**3)
C   value across analysis domain.  Second constant is scale height.
 
      Dens(z) = 1.1904 / Exp(z/9.58125)
  
      If (Idir_Int .eq. +1) Then             ! upward integration
         If (k .eq. 1) Then                  ! extend column down to sfc
            Kappa = (Dens(0.) - Dens(Z0))/(Dens(Z0/2.))/2.0
            Beta  = 1.0 + Kappa
            Alpha = 1.0 - Kappa
            w = Wl*Beta/Alpha - SynthDat(i,j,1,5) * Z0m/Alpha
         Else
            If (k .eq. IFirst) Then
               Kappa=(Dens(Z0+(Float(k)-1.5)*Sz)-Dens(Z0+(k-1)*Sz))/
     1              Dens(Z0 + (Float(k)-1.25)*Sz)/2.0
               Beta  = 1.0 + Kappa
               Alpha = 1.0 - Kappa
               w = Wl*Beta/Alpha - SynthDat(i,j,k,5) * Szm/2.0/Alpha
            Else
               Kappa = (Dens(Z0 + (k-2)*Sz) - Dens(Z0 + (k-1)*Sz))/
     1              Dens(Z0 + (Float(k)-1.5)*Sz)/2.0
               Beta  = 1.0 + Kappa
               Alpha = 1.0 - Kappa
               w = Wl*Beta/Alpha - (SynthDat(i,j,k,5) + SynthDat(i,j,k-1,5))*
     1             Szm/Alpha/2.0
            End If
         End If
      Else                                   ! downward integration
         If (k .eq. IFirst) Then
  
c  Apply upper b.c. one-half grid level above top-most data point
 
            Kappa=(Dens(Z0+Float(k-1)*Sz)-Dens(Z0+(Float(k)-0.5)*Sz))/
     1           Dens(Z0 + (Float(k)-0.75)*Sz)/2.0
            Beta  = 1.0 - Kappa
            Alpha = 1.0 + Kappa
            w = Wl*Beta/Alpha + SynthDat(i,j,k,5) * Szm/2.0/Alpha
         Else
            Kappa = (Dens(Z0 + Float(k-1)*Sz) - Dens(Z0 + k*Sz))/
     1               Dens(Z0 + (Float(k)-0.5)*Sz)/2.0
            Beta  = 1.0 - Kappa
            Alpha = 1.0 + Kappa
            w = Wl*Beta/Alpha + (SynthDat(i,j,k+1,5) + SynthDat(i,j,k,5)) *
     1          Szm/Alpha/2.0
         End If
      End If
  
      Return
      End Subroutine W_Calc
      
c     ***********************************************

      Subroutine Triple (SynthDat, Radar, Type)
  
C  This Subroutine computes (u,v,w) using the Matjeka
c  overdetermined triple equation solution technique.
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /Print/ LuMes
      Common /Limits/ Zsfc, Iters, Kup, Klim, Top_Hgt

      Dimension SynthDat(Imax,Jmax,Kmax,Nflds)
      Dimension Radar(6,Imax,Jmax,Kmax,Nrdrs)
      Character*3 Type(Nrdrs)

c  Calculate the reflectivity

      Call ChooZ (SynthDat, Radar)

c  Create the wind vectors from the radial velocities

      Call Comp_Wind (SynthDat, Radar, Type)

c  Fill small holes in the horizonatal planes

      If (Ihole_Fill .eq. 1) Then
         Call Hole_Fill(SynthDat, 1)
      End If

c  Fill small gaps in the vertical columns

      If (Klim .gt. 0) Then
         Call Vfill (SynthDat,1,Nu_Gaps) ! u
         Call VFill (SynthDat,6,Nu_Gaps) ! SE u
         Call Vfill (SynthDat,2,Nv_Gaps) ! v
         Call Vfill (SynthDat,7,Nv_Gaps) ! SE v
         Call Vfill (SynthDat,4,Nw_Gaps) ! w
         Call VFill (SynthDat,8,Nw_Gaps) ! SE w
         Call VFill (SynthDat,3,Nz_Gaps) ! dBz
      Else
         Nu_Gaps = 0
         Nv_Gaps = 0
         Nw_Gaps = 0
         Nz_Gaps = 0
      End If

c  Smooth the u,v,w,Vt fields (if desired)

      If (Nsmth .gt. 0 .and. Ismtyp .gt. 0)
     #     Call Smooth_It(SynthDat)

c  Print out the column statistics

      Write (LuMes,200) Nu_Gaps, Nv_Gaps, Nw_Gaps, Nz_Gaps
 200  Format (/'Number of vertical column u gaps interpolated over:'
     #      i5/'Number of vertical column v gaps interpolated over:'
     #      i5/'Number of vertical column w gaps interpolated over:'
     #      i5/'Number of vertical column dBZ gaps interpolated over:'
     #      i5)

      Return

      End Subroutine Triple
      
c     ***********************************************
  
      Subroutine Init_Field (SynthDat)
  
c  This routine initializes data to no data flag
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata

      Dimension SynthDat(Imax,Jmax,Kmax,Nflds)  
  
      Do k = 1, Kmax
         Do j = 1, Jmax
            Do i = 1, Imax
               Do n = 1, Nflds
                  SynthDat (i,j,k,n) = Baddata
               End Do
            End Do
         End Do
      End Do

      Return
      End Subroutine Init_Field
      
c     ***********************************************

      Subroutine Examine_Column (X, Itop, Itop_Miss, Ibot_Miss)
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata

      Common /Print/ LuMes

      Dimension X(*)

c  Find the height of the echo top

      Itop = -99
      
      Do k = Kmax, 1, -1
         If (X(k) .ne. Baddata) Then
            Itop = k
            Go To 10
         End If
      End Do
      
c     Find the highest missing below the top level
      
 10   Itop_Miss = -99
      
      Do k = Itop, 1, -1
         If (X(k) .eq. Baddata) Then
            Itop_Miss = k
            Go To 20
         End If
      End Do
      
c     Find the lowest missing level above the surface
      
 20   Ibot_Miss = -99
      
      Do k = 1, Itop
         If (X(k) .eq. Baddata) Then
            Ibot_Miss = k
            Go To 30
         End If
      End Do
      
 30   Return
      End Subroutine Examine_Column
      
c     ***********************************************

      Subroutine Hole_Fill(SynthDat, Iprint)

c  Routine to "hole fill" the u,v,w,dBZ, and Vt arrays
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /Print/ LuMes
      Common /HoleFillParms/ Max_Cons_Unoccupied_Octants,
     #                       Max_Search_Radius,
     #                       Min_Values_Per_Octant, 
     #                       Max_Values_Per_Octant

      Dimension SynthDat(Imax,Jmax,Kmax,Nflds)
      Dimension Dum(Imax,Jmax)

      Logical Candidate(Imax,Jmax)

c  Initialize the candidate array

      Do j = 1, Jmax
         Do i = 1, Imax
            Candidate(i,j) = .True.
         End Do
      End Do

c  Define a few constants

      Xmin = 0.0
      Ymin = 0.0
      Delx = Sx
      Dely = Sy

c  Fill holes a plane at a time

      Do k = 1, Kmax

c  Start with the U field

         Call Hole_Fill_2D(
     #        SynthDat(1,1,k,1), SynthDat(1,1,k,6), Candidate, 
     #        Imax, Jmax, Delx, Dely,
     #        Imax, Jmax, Baddata,
     #        Max_Cons_Unoccupied_Octants,
     #        Max_Search_Radius,
     #        Min_Values_Per_Octant,
     #        Max_Values_Per_Octant,
     #        Imax, Jmax,
     #        SynthDat(1,1,k,1), SynthDat(1,1,k,6), 
     #        Nmissing, Nfilled)

         If (Iprint .eq. 1)
     #        Write (LuMes,'(/"u  array.  On level:",i3,
     #        " filled ",i5," values out of ",i5,
     #        " missing data points")') k, Nfilled, Nmissing

c  v field

         Call Hole_Fill_2D(
     #        SynthDat(1,1,k,2), SynthDat(1,1,k,7), Candidate, 
     #        Imax, Jmax, Delx, Dely,
     #        Imax, Jmax, Baddata,
     #        Max_Cons_Unoccupied_Octants,
     #        Max_Search_Radius,
     #        Min_Values_Per_Octant,
     #        Max_Values_Per_Octant,
     #        Imax, Jmax,
     #        SynthDat(1,1,k,2), SynthDat(1,1,k,7), 
     #        Nmissing, Nfilled)

         If (Iprint .eq. 1)
     #        Write (LuMes,'("v  array.  On level:",i3,
     #        " filled ",i5," values out of ",i5,
     #        " missing data points")') k, Nfilled, Nmissing

c  w array

         Call Hole_Fill_2D(
     #        SynthDat(1,1,k,4), SynthDat(1,1,k,8), Candidate, 
     #        Imax, Jmax, Delx, Dely,
     #        Imax, Jmax, Baddata,
     #        Max_Cons_Unoccupied_Octants,
     #        Max_Search_Radius,
     #        Min_Values_Per_Octant,
     #        Max_Values_Per_Octant,
     #        Imax, Jmax,
     #        SynthDat(1,1,k,4), SynthDat(1,1,k,8), 
     #        Nmissing, Nfilled)

         If (Iprint .eq. 1)
     #        Write (LuMes,'("w  array.  On level:",i3,
     #        " filled ",i5," values out of ",i5,
     #        " missing data points")') k, Nfilled, Nmissing

c  Vt array

         Call Hole_Fill_2D(
     #        SynthDat(1,1,k,9), SynthDat(1,1,k,11), Candidate, 
     #        Imax, Jmax, Delx, Dely,
     #        Imax, Jmax, Baddata,
     #        Max_Cons_Unoccupied_Octants,
     #        Max_Search_Radius,
     #        Min_Values_Per_Octant,
     #        Max_Values_Per_Octant,
     #        Imax, Jmax,
     #        SynthDat(1,1,k,9), SynthDat(1,1,k,11), 
     #        Nmissing, Nfilled)

         If (Iprint .eq. 1)
     #        Write (LuMes,'("Vt array.  On level:",i3,
     #        " filled ",i5," values out of ",i5,
     #        " missing data points")') k, Nfilled, Nmissing

c  dBZ array, first dummy up the SD of dBZ

         Do i = 1, Imax
            Do j = 1, Jmax
               Dum(i,j) = Baddata
               If (SynthDat(i,j,k,3) .ne. Baddata) Dum(i,j) = 0.5
            End Do
         End Do

         Call Hole_Fill_2D(
     #        SynthDat(1,1,k,3), Dum(1,1), Candidate, 
     #        Imax, Jmax, Delx, Dely,
     #        Imax, Jmax, Baddata,
     #        Max_Cons_Unoccupied_Octants,
     #        Max_Search_Radius,
     #        Min_Values_Per_Octant,
     #        Max_Values_Per_Octant,
     #        Imax, Jmax,
     #        SynthDat(1,1,k,3), Dum, 
     #        Nmissing, Nfilled)

         If (Iprint .eq. 1)
     #        Write (LuMes,'("Z  array.  On level:",i3,
     #        " filled ",i5," values out of ",i5,
     #        " missing data points")') k, Nfilled, Nmissing
      End Do

      Return
      End Subroutine Hole_Fill
      
c     ***********************************************

      Subroutine VFill (SynthDat, n, Ngaps)

c  Routine to interpolate over small gaps in the array (1 missing level max)
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /Limits/ Zsfc, Iters, Kup, Klim, Top_Hgt

      Dimension SynthDat (Imax,Jmax,Kmax,Nflds)
  
      Ngaps = 0

      Do j = 1, Jmax
         Do i = 1, Imax
            Do k = 2, Kmax
               Val  = SynthDat(i,j,k,n)
               Val1 = SynthDat(i,j,k-1,n)

               If (Val .eq. Baddata .and. Val1 .ne. Baddata) Then
                  Do kk = k, Klim+k
                     If (kk .gt. Kmax) Cycle

                     If (SynthDat(i,j,kk,n) .ne. Baddata) Then
                        Grad = SynthDat(i,j,kk,n) - SynthDat(i,j,k-1,n)
                        Ngood = kk - k + 1
                        Knt = 0

                        Do kkk = k, kk-1
                           Knt = Knt + 1
                           Fac = Float(Knt)/Float(Ngood)
                           SynthDat(i,j,kkk,n) = SynthDat(i,j,k-1,n) + 
     #                          Fac * Grad
                           Ngaps = Ngaps + 1
                        End Do

                        Cycle
                     End If
                  End Do
               End If
            End Do
         End Do
      End Do

      Return
      End Subroutine VFill
      
c     ***********************************************

      Subroutine Smooth_It (SynthDat)
  
C  Smooth the (u,v,w,Vt) fields

      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /Print/ LuMes

      Dimension SynthDat(Imax,Jmax,Kmax,Nflds)

      Dimension Wt_X (-100:100), Wt_Y (-100:100)

      Logical Uncorrelated

      Max_Reach = 100

c  Set up the weights for the binomial and Leise filters

      If (IsmTyp .eq. 1) Then          ! the "Leise" filter
         Call Multi_Step_Leise_Filter_Weights (Nsmth, Max_Reach,
     #        Wt_X, Irx)
         Call Multi_Step_Leise_Filter_Weights (Nsmth, Max_Reach,
     #        Wt_Y, Iry)           
      Else If (IsmTyp .eq. 2) Then     ! the "binominal" filter
         Call Multi_Step_Binomial_Filter_Weights(Nsmth, Max_Reach,
     #        Wt_X, Irx)
         Call Multi_Step_Binomial_Filter_Weights(Nsmth, Max_Reach,
     #        Wt_Y, Iry)           
      Else
         Return
      End If

c  Load up the arrays and perform the filtering

      Do k = 1, Kmax

c  First, the u field

         Uncorrelated = .True.

         Call Filter_2D(
     #        SynthDat(1,1,k,1),             ! original u array
     #        SynthDat(1,1,k,6),             ! standary dev of u array
     #        Imax, Jmax,       ! dimensions of array
     #        Imax, Jmax,       ! size of grid
     #        Baddata,          ! bad flag
     #        Uncorrelated,     ! correlated data?
     #        Wt_X,             ! x filter weights
     #        Max_Reach,        ! number of weights
     #        Irx,              ! reach of x
     #        WT_Y,             ! y filter weights
     #        Iry,              ! reach of y
     #        Imax, Jmax,       ! size of the output arry
     #        SynthDat(1,1,k,1),            ! filtered u field
     #        SynthDat(1,1,k,6))            ! filtered standard dev of u

c  Next the v field

         Uncorrelated = .True.

         Call Filter_2D(
     #        SynthDat(1,1,k,2),             ! original v array
     #        SynthDat(1,1,k,7),             ! standary dev of v array
     #        Imax, Jmax,       ! dimensions of array
     #        Imax, Jmax,       ! size of grid
     #        Baddata,          ! bad flag
     #        Uncorrelated,     ! correlated data?
     #        Wt_X,             ! x filter weights
     #        Max_Reach,        ! number of weights
     #        Irx,              ! reach of x
     #        WT_Y,             ! y filter weights
     #        Iry,              ! reach of y
     #        Imax, Jmax,       ! size of the output arry
     #        SynthDat(1,1,k,2),            ! filtered v field
     #        SynthDat(1,1,k,7))            ! filtered standard dev of v

c  The w field

         Uncorrelated = .True.

         Call Filter_2D(
     #        SynthDat(1,1,k,4),             ! original w array
     #        SynthDat(1,1,k,8),             ! standary dev of w array
     #        Imax, Jmax,       ! dimensions of array
     #        Imax, Jmax,       ! size of grid
     #        Baddata,          ! bad flag
     #        Uncorrelated,     ! correlated data?
     #        Wt_X,             ! x filter weights
     #        Max_Reach,        ! number of weights
     #        Irx,              ! reach of x
     #        WT_Y,             ! y filter weights
     #        Iry,              ! reach of y
     #        Imax, Jmax,       ! size of the output arry
     #        SynthDat(1,1,k,4),            ! filtered w field
     #        SynthDat(1,1,k,8))            ! filtered standard dev of w

c  Lastly, the Vt field

         Uncorrelated = .False.

         Call Filter_2D(
     #        SynthDat(1,1,k,9),             ! original Vt array
     #        SynthDat(1,1,k,11),            ! standary dev of Vt array
     #        Imax, Jmax,       ! dimensions of array
     #        Imax, Jmax,       ! size of grid
     #        Baddata,          ! bad flag
     #        Uncorrelated,     ! correlated data?
     #        Wt_X,             ! x filter weights
     #        Max_Reach,        ! number of weights
     #        Irx,              ! reach of x
     #        WT_Y,             ! y filter weights
     #        Iry,              ! reach of y
     #        Imax, Jmax,       ! size of the output arry
     #        SynthDat(1,1,k,9),            ! filtered Vt field
     #        SynthDat(1,1,k,11))           ! filtered standard dev of Vt
      End Do

      Return
      End Subroutine Smooth_It
      
c     ***********************************************

      Subroutine Calc_Vt (i,j,k,SynthDat)

c  Routine to estimate reflectivity-weighted terminal fallspeeds from
c  estimates of reflectivity
c   Based on empirical relationships between drop size and terminal
c   fallspeeds

      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /Print/ LuMes
      Common /Sounding/ Z1, Tv1, Z2, Tv2, Z3, Rho3

      Dimension SynthDat(Imax,Jmax,Kmax,Nflds)

c  Find the terminal velocity at the grid point

      Hgt = Hgt_of_Level(Z0,k,SZ) * 1000.0

      Ref = SynthDat(i,j,k,3)

c  Term_Vel returns a value for terminal fallspeed based on dbZ

      If (Ref .ne. Baddata) Then
         HRain = Vt_Rain * 1000.0
         HSnow = Vt_Snow * 1000.0
         Call F_From_Dbz (Z1, Tv1, Z2, Tv2, Z3, Rho3, HRain,
     #        HSnow, Hgt, Ref, Term_Vel, Term_Vel_SD)
      Else
         Term_Vel    = Baddata
         Term_Vel_SD = Baddata
      End If

      SynthDat(i,j,k,9)  = -Term_Vel
      SynthDat(i,j,k,11) =  Term_Vel_SD

      Return
      End Subroutine Calc_Vt
      
c     ***********************************************

      Function Hgt_of_Level (Z0, k, Sz)

c  Calculate the height (km) of any given level
c    Z0 = height of first level (km)
c    k  = index number of the level
c    Sz = Spacing between levels (km)

      Hgt_of_Level = Z0 + Float(k-1)*Sz

      Return
      End Function Hgt_of_Level
      
c     ***********************************************

      Subroutine Adj_Winds_from_Div (SynthDat)
c  Routine to adjust the (u,v) winds for the adjusted divergence
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /Print/ LuMes
      Common /HoleFillParms/ Max_Cons_Unoccupied_Octants,
     #                       Max_Search_Radius,
     #                       Min_Values_Per_Octant, 
     #                       Max_Values_Per_Octant

      Dimension SynthDat(Imax,Jmax,Kmax,Nflds)
      Dimension U_Out(Imax,Jmax), V_Out(Imax,Jmax)

      Parameter (Ind_Max = 30)

      Dimension Ihist(Ind_Max)

      Data Steps /0.20/

      Do l = 1, Ind_Max
         Ihist(l) = 0
      End Do

      Delx = Sx * 1000.0           ! x grid spacing in m
      Dely = Sy * 1000.0           ! y grid spacing in m
      Vol_Total = 0.0
      N_Lev = 0

      Write (LuMes,'(/"u,v adjuster for divergence:"/
     #     " n   ave u dif Npts   Largest u   ave v dif   ",
     #     "Largest v  wind diff  wind dif sd")')

c  Do the adjustment a plane at a time

      Do k = 1, Kmax
         Call Adjust_Winds_to_Div_2D (
     #        SynthDat(1,1,k,1), SynthDat(1,1,k,6),
     #        SynthDat(1,1,k,2), SynthDat(1,1,k,7), SynthDat(1,1,k,5),
     #        Imax, Jmax, Delx, Dely, Imax, Jmax, 
     #        Baddata, Imax, Jmax,
     #        U_Out, V_Out)

         Ularge = 0.0
         Vlarge = 0.0
         Utot = 0.0
         Vtot = 0.0
         Npts = 0
         Grid_Total = 0.0
         UVtot = 0
         UVss = 0.0

c  Calculate some statistics of the amount of change in the wind
c    field due to the divergence adjustment

         Do j = 1, Jmax
            Do i = 1, Imax
               If ( SynthDat(i,j,k,1) .ne. Baddata .and. 
     #                 U_Out(i,j) .ne. Baddata .and.
     #              SynthDat(i,j,k,2) .ne. Baddata .and.
     #                 V_Out(i,j) .ne. Baddata) Then
                  Udif = Abs(SynthDat(i,j,k,1) - U_Out(i,j))   ! u difference
                  If (Udif .gt. Ularge) Ularge = Udif      ! biggest u change
                  Usq = Udif*Udif                          ! u diff squared
                  Utot = Utot + Udif                       ! u total diff
                  Npts = Npts + 1                          ! count the point
                  Vdif = Abs(SynthDat(i,j,k,2) - V_Out(i,j))   ! v difference
                  If (Vdif .gt. Vlarge) Vlarge = Vdif      ! biggest v change
                  Vsq = Vdif*Vdif                          ! v diff squared
                  Vtot = Vtot + Vdif                       ! v total diff

                  Var = Sqrt(Usq + Vsq)                    ! total wind change
                  UVtot = UVtot + Var                      ! total wind diff
                  UVss = UVss + Usq + Vsq
                  Ind = Var/Steps + 1.5                    ! index for histagram
                  If (Ind .gt. Ind_Max) Ind = Ind_Max      ! cut off the upper limit
                  Ihist(Ind) = Ihist(Ind) + 1              ! count the point
               End If

c  Put the (u,v) data back into the master array

               SynthDat(i,j,k,1) = U_Out(i,j)
               SynthDat(i,j,k,2) = V_Out(i,j)
            End Do
         End Do

         If (Npts .gt. 0) Then
            Uave = Utot/Float(Npts)
            Vave = Vtot/Float(Npts)
            UVave= UVtot/Float(Npts)
            UVsd = Sqrt(UVss/Float(Npts) - UVave*UVave)
            Vol_Total = Vol_Total + UVave
            N_Lev = N_Lev + 1
            Write (LuMes,'(i2,E12.5,i6,5E12.5)')
     #           k, Uave, Npts, Ularge, Vave, Vlarge, UVave, UVsd
         End If
      End Do

      If (N_Lev .gt. 0) Then
         Vol_Ave = Vol_Total/Float(N_Lev)
         Write (LuMes,'(/"Volume averaged change:",E12.5,i5)')
     #        Vol_Ave, N_Lev
         Write (LuMes,'(/"Distribution of changes:"/
     #        "Abs Dif (m/s)    Number")')

         Do l=1, Ind_Max
            H_low = Float(l-1) * Steps
            H_high = Float(l) * Steps
            Write (LuMes,'(F4.1,"-",F4.1,7x,i6)')
     #           H_low, H_high, Ihist(l)
         End Do
      End If

      Return
      End Subroutine Adj_Winds_from_Div
      
c     ***********************************************

      Block Data
  
      Common /Size/ Imax, Jmax, Kmax, Nrdrs, Istyle, Ctr, ThresH, Nflds,
     #     ThresV, Sx, Sy, Sz, Z0, Nsmth, IsmTyp, Iobr, Ihole_Fill,
     #     Su, Sv, Idbz_Parm, Idir_Int, VT_Snow, VT_Rain, Baddata
      Common /Vert/ IFirst, Last, Szm, Z0m
      Common /Print/ LuMes
      Common /Sounding/ Z1, Tv1, Z2, Tv2, Z3, Rho3
      Common /Limits/ Zsfc, Iters, Kup, Klim, Top_Hgt
      Common /HoleFillParms/ Max_Cons_Unoccupied_Octants,
     #                       Max_Search_Radius,
     #                       Min_Values_Per_Octant, 
     #                       Max_Values_Per_Octant
      Common /Debug/ Idebug, Ix1, Ix2, Jy1, Jy2
      Common /w_Techniques/ Iw_at_top, I_w_0
 
      Data LuMes /60/,
     #     Nsmth /0/,
     #     IsmTyp /0/,
     #     Klim /1/,
     #     Iobr /1/, Kup /1/,
     #     Idebug /0/,
     #     Idir_Int /-99/,
     #     Idbz_Parm /-1/,
     #     Baddata /-999.0/
     #     Ctr /0.01745329/,
     #     Iters /1/,
     #     Iw_at_top /1/, I_w_0 /1/,
     #     Top_Hgt /0.6/,
     #     Ix1 /-1/, Ix2 /-1/, Jy1 /-1/, Jy2 /-1/
      End
