      Function VCOMP(WS,WD)

C  Returns THE V COMPONENT (NORTHERLY) OF THE WIND
C  WS  = WIND SPEED                               (M/S)
C  WD  = WIND DIRECTION                           (DEG)

       A = ANGLE(WD) + 180.0
       If (A.gt.360.0) A = A - 360.0
       VCOMP = WS * SIN(A*0.01745329)
       Return
       End
