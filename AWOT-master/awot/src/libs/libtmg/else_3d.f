      SUBROUTINE ELSE_3D(A_IN,B_IN,MAXX_IN,MAXY_IN,MAXZ_IN,NX,NY,NZ,
     $BADDATA,MAXX_OUT,MAXY_OUT,MAXZ_OUT,A_OUT)

C  Thomas Matejka NOAA/NSSL 6 April 1994

C  This subroutine creates a three-dimensional data field that consists
C  at each point of the value, if it exists, at that point in a first
C  input field, or, if it does not, the value at that point in a second
C  input field.

C  Input:

C  A_IN is a three-dimensional real array.  A_IN(I,J,K) specifies the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, and the Kth grid point in the third
C  dimension, of the data field to be copied to the output field
C  wherever it exists.  If it is missing, it should equal BADDATA.

C  B_IN is a three-dimensional real array.  B_IN(I,J,K) specifies the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, and the Kth grid point in the third
C  dimension, of the data field to be copied to the output field
C  wherever A_IN does not exist.  If it is missing, it should equal
C  BADDATA.

C  MAXX_IN is an integer variable that specifies the first dimension of
C  A_IN and B_IN in the calling program.

C  MAXY_IN is an integer variable that specifies the second dimension of
C  A_IN and B_IN in the calling program.

C  MAXZ_IN is an integer variable that specifies the third dimension of
C  A_IN and B_IN in the calling program.

C  NX is an integer variable that specifies the number of grid points in
C  the first dimension for A_IN, B_IN, and A_OUT.

C  NY is an integer variable that specifies the number of grid points in
C  the second dimension for A_IN, B_IN, and A_OUT.

C  NZ is an integer variable that specifies the number of grid points in
C  the third dimension for A_IN, B_IN, and A_OUT.

C  BADDATA is a real variable that indicates a missing value as
C  described.

C  MAXX_OUT is an integer variable that specifies the first dimension of
C  A_OUT in the calling program.

C  MAXY_OUT is an integer variable that specifies the second dimension
C  of A_OUT in the calling program.

C  MAXZ_OUT is an integer variable that specifies the third dimension of
C  A_OUT and the calling program.

C  Output:

C  A_OUT is a three-dimensional real array.  A_OUT(I,J,K) returns the
C  value, at the Ith grid point in the first dimension, the Jth grid
C  point in the second dimension, and the Kth grid point in the third
C  dimension, of A_IN(I,J,K) wherever it exists.  Otherwise, it returs
C  the value of B_IN(I,J,K).  If it is missing, it is returned as
C  BADDATA.  A_OUT may be identical to A_IN or B_IN, in which case A_IN
C  or B_IN is overwritten.

      IMPLICIT NONE
      INTEGER MAXX_IN,MAXY_IN,MAXZ_IN,MAXX_OUT,MAXY_OUT,MAXZ_OUT,NX,NY,
     $NZ,IX,IY,IZ
      REAL BADDATA
      REAL A_IN(MAXX_IN,MAXY_IN,MAXZ_IN),B_IN(MAXX_IN,MAXY_IN,MAXZ_IN)
      REAL A_OUT(MAXX_OUT,MAXY_OUT,MAXZ_OUT)

C  Loop through all the points.
      DO 1 IZ=1,NZ
         DO 2 IY=1,NY
            DO 3 IX=1,NX

C  Find the value of A_IN else B_IN.
               IF(A_IN(IX,IY,IZ).NE.BADDATA)THEN
                  A_OUT(IX,IY,IZ)=A_IN(IX,IY,IZ)
               ELSE
                  A_OUT(IX,IY,IZ)=B_IN(IX,IY,IZ)
               ENDIF
3           CONTINUE
2        CONTINUE
1     CONTINUE

C  Done.
      RETURN
      END
