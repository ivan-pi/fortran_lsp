      program PROG4
c
C  DEMONSTRATE SINGULAR VALUE ANALYSIS.
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------
      integer MDA, MX
      parameter(MDA = 15, MX = 5)
      integer I, J, KPVEC(4), M, N
      double precision A(MDA,MX), B(MDA), D(MX), SING(MX), WORK(2*MX)
      character NAMES(MX)*8
      data KPVEC / 1, 111111, -1, 76 /
      data NAMES/'  Fire',' Water',' Earth','   Air','Cosmos'/
C     ------------------------------------------------------------------
      M = MDA
      N = MX
      write (*,'(a)')
     * ' PROG4.  DEMONSTRATE SINGULAR VALUE ANALYSIS',
     * ' PROG4 will read data from file: data4.dat'
      open(unit= 10, file= 'data4.dat', status= 'old')
      read (10,'(6F12.0)') ((A(I,J),J=1,N),B(I),I=1,M)
      write (*,'(/a)')
     * ' LISTING OF INPUT MATRIX, A, AND VECTOR, B, FOLLOWS..'
      write (*,'(/(1x,f11.7,4F12.7,F13.4))') ((A(I,J),J=1,N),B(I),I=1,M)
      write (*,'(/)')
C
      call SVA (A,MDA,M,N,MDA,B,SING,KPVEC, NAMES,1,D, WORK)
C
      end