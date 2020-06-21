      double precision FUNCTION GEN(ANOISE)
c
C  GENERATE NUMBERS FOR CONSTRUCTION OF TEST CASES.
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1972 DEC 15, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
c     ------------------------------------------------------------------
      integer I, J, MI, MJ
      double precision AI, AJ, ANOISE, ZERO
      parameter(ZERO = 0.0d0)
      SAVE
c     ------------------------------------------------------------------
      IF (ANOISE) 10,30,20
   10 MI=891
      MJ=457
      I=5
      J=7
      AJ= ZERO
      GEN= ZERO
      RETURN
C
C     THE SEQUENCE OF VALUES OF J  IS BOUNDED BETWEEN 1 AND 996
C     IF INITIAL J = 1,2,3,4,5,6,7,8, OR 9, THE PERIOD IS 332
   20 J=J*MJ
      J=J-997*(J/997)
      AJ=J-498
C     THE SEQUENCE OF VALUES OF I  IS BOUNDED BETWEEN 1 AND 999
C     IF INITIAL I = 1,2,3,6,7, OR 9,  THE PERIOD WILL BE 50
C     IF INITIAL I = 4 OR 8   THE PERIOD WILL BE 25
C     IF INITIAL I = 5        THE PERIOD WILL BE 10
   30 I=I*MI
      I=I-1000*(I/1000)
      AI=I-500
      GEN=AI+AJ*ANOISE
      RETURN
      END