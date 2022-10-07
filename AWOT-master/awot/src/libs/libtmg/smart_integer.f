      SUBROUTINE SMART_INTEGER(I,STRING_OUT)

C  Thomas Matejka NOAA/NSSL 22 May 2000

      IMPLICIT NONE
      INCLUDE 'tmmlib.inc'
      INTEGER,EXTERNAL::S_L
      CHARACTER(LEN=MAX_STRING)::STRING
      CHARACTER(LEN=*)::STRING_OUT
      INTEGER::I,L,IEND

C  Initialize.
      STRING=''
      STRING_OUT=''
      L=LEN(STRING_OUT)

      WRITE(STRING,*)I
      IF(S_L(STRING).GE.MAX_STRING)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'SMART_INTEGER: MEMORY EXCEEDED. ',
     $   'INCREASE MAX_STRING.'
         STOP
      ENDIF
      CALL LEFT_JUSTIFY(STRING)
      IEND=S_L(STRING)
      IF(IEND.GT.L)THEN
         WRITE(TMMLIB_MESSAGE_UNIT,*)'SMART_INTEGER: STRING_OUT IS ',
     $   'TOO SHORT.'
         STOP
      ENDIF
      STRING_OUT(1:IEND)=STRING(1:IEND)

      END SUBROUTINE SMART_INTEGER
