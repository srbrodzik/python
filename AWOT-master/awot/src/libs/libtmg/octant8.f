      SUBROUTINE OCTANT8(A_IN,SD_A_IN,MAXX_IN,MAXY_IN,NX,NY,BADDATA,
     $MAX_SEARCH_RADIUS,MAX_VALUES_PER_OCTANT,IX,IY,POINTS_FOUND,
     $IX_FOUND,IY_FOUND)

      IMPLICIT NONE
      INTEGER MAXX_IN,MAXY_IN,NX,NY,IX,IY,IX_SEARCH,IY_SEARCH,
     $POINTS_FOUND,MAX_SEARCH_RADIUS,MAX_VALUES_PER_OCTANT
      INTEGER IX_FOUND(1),IY_FOUND(1)
      REAL BADDATA
      REAL A_IN(MAXX_IN,MAXY_IN),SD_A_IN(MAXX_IN,MAXY_IN)

      POINTS_FOUND=0
      DO 1 IY_SEARCH=IY+1,IY+MAX_SEARCH_RADIUS
         IF(IY_SEARCH.GT.NY)THEN
             RETURN
         ENDIF
         DO 2 IX_SEARCH=IX-1,IX-(IY_SEARCH-IY),-1
            IF(IX_SEARCH.LT.1)THEN
               GOTO 3
            ENDIF
            IF(A_IN(IX_SEARCH,IY_SEARCH).NE.BADDATA.AND.
     $      SD_A_IN(IX_SEARCH,IY_SEARCH).NE.BADDATA)THEN
               POINTS_FOUND=POINTS_FOUND+1
               IX_FOUND(POINTS_FOUND)=IX_SEARCH
               IY_FOUND(POINTS_FOUND)=IY_SEARCH
               IF(POINTS_FOUND.GE.MAX_VALUES_PER_OCTANT)THEN
                  RETURN
               ENDIF
            ENDIF
2        CONTINUE
3        CONTINUE
1     CONTINUE
      RETURN
      END
