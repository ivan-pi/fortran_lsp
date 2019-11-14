      program PROG2
c
C  DEMONSTRATE ALGORITHM  HFTI  FOR SOLVING LEAST SQUARES PROBLEMS
C  AND ALGORITHM  COV  FOR COMPUTING THE ASSOCIATED UNSCALED
C  COVARIANCE MATRIX.
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------
      integer MDA
      parameter(MDA = 8)
      integer I, II, IP(MDA), IP1, J, JM1, K, KM1, KP1, KRANK, L
      integer M, MN1, MN2, N, NM1, NOISE
      double precision A(MDA,MDA), ANOISE, ANORM, B(MDA)
      double precision DUMMY, FIVE00, G(MDA), GEN
      double precision H(MDA), HALF
      double precision ONE, SM, SMALL, SRSMSQ, TAU, TMP, ZERO
      parameter(FIVE00 = 500.0d0, HALF = 0.5d0)
      parameter(ONE = 1.0d0, SMALL = 1.0d-4, ZERO = 0.0d0)
C     ------------------------------------------------------------------
  190 format (/1x,8X,'RESIDUAL LENGTH = ',E12.4)
  200 format (/1x,8X,
     * 'ESTIMATED PARAMETERS,  X=A**(+)*B, COMPUTED BY ''HFTI'''//
     * (9X,I6,E16.8,I6,E16.8,I6,E16.8,I6,E16.8,I6,E16.8))
  210 format (/1x,8X,
     * 'COVARIANCE MATRIX (UNSCALED) OF ESTIMATED PARAMETERS'/
     * 9x,'COMPUTED BY ''COV''.'/1X)
  220 format (9X,2I3,E16.8,2I3,E16.8,2I3,E16.8,2I3,E16.8,2I3,E16.8)
  230 format (/' PROG2.  THIS PROGRAM DEMONSTATES THE ALGORITHMS'/
     * 9x,'HFTI  AND  COV.')
  240 format (/
     * ' THE RELATIVE NOISE LEVEL OF THE GENERATED DATA WILL BE ',
     * g11.3/
     * ' THE MATRIX NORM IS APPROXIMATELY ',E12.4/
     * ' THE ABSOLUTE PSEUDORANK TOLERANCE, TAU, IS ',E12.4)
  250 format (/////'    M   N'/1X,2I4)
  260 format (/1x,8X,'PSEUDORANK = ',I4)
C     ------------------------------------------------------------------
          DO 180 NOISE=1,2
          ANORM= FIVE00
          ANOISE= ZERO
          TAU= HALF
          IF (NOISE.EQ.1) GO TO 10
          ANOISE= SMALL
          TAU=ANORM*ANOISE*10.
   10     CONTINUE
C     INITIALIZE THE DATA GENERATION FUNCTION
C    ..
          DUMMY=GEN(-ONE)
          write (*,230)
          write (*,240) ANOISE,ANORM,TAU
C
              DO 180 MN1=1,6,5
              MN2=MN1+2
                  DO 180 M=MN1,MN2
                      DO 180 N=MN1,MN2
                      write (*,250) M,N
C     GENERATE DATA
C    ..
                          DO 20 I=1,M
                              DO 20 J=1,N
   20                         A(I,J)=GEN(ANOISE)
                          DO 30 I=1,M
   30                     B(I)=GEN(ANOISE)
C
C     ****** CALL HFTI   ******
C
                      CALL HFTI(A,MDA,M,N,B,1,1,TAU,KRANK,SRSMSQ,H,G,IP)
C
                      write (*,260) KRANK
                      write (*,200) (I,B(I),I=1,N)
                      write (*,190) SRSMSQ
                      IF (KRANK.LT.N) GO TO 180
C     ******  ALGORITHM COV BEGINS HERE  ******
C    ..
                          DO 40 J=1,N
   40                     A(J,J)= ONE/A(J,J)
                      IF (N.EQ.1) GO TO 70
                      NM1=N-1
                          DO 60 I=1,NM1
                          IP1=I+1
                              DO 60 J=IP1,N
                              JM1=J-1
                              SM= ZERO
                                  DO 50 L=I,JM1
   50                             SM=SM+A(I,L)*A(L,J)
   60                         A(I,J)=-SM*A(J,J)
C    ..
C     THE UPPER TRIANGLE OF A HAS BEEN INVERTED
C     UPON ITSELF.
   70                     DO 90 I=1,N
                              DO 90 J=I,N
                              SM= ZERO
                                  DO 80 L=J,N
   80                             SM=SM+A(I,L)*DBLE(A(J,L))
   90                         A(I,J)=SM
                      IF (N.LT.2) GO TO 160
                          DO 150 II=2,N
                          I=N+1-II
                          IF (IP(I).EQ.I) GO TO 150
                          K=IP(I)
                          TMP=A(I,I)
                          A(I,I)=A(K,K)
                          A(K,K)=TMP
                          IF (I.EQ.1) GO TO 110
                              DO 100 L=2,I
                              TMP=A(L-1,I)
                              A(L-1,I)=A(L-1,K)
  100                         A(L-1,K)=TMP
  110                     IP1=I+1
                          KM1=K-1
                          IF (IP1.GT.KM1) GO TO 130
                              DO 120 L=IP1,KM1
                              TMP=A(I,L)
                              A(I,L)=A(L,K)
  120                         A(L,K)=TMP
  130                     IF (K.EQ.N) GO TO 150
                          KP1=K+1
                              DO 140 L=KP1,N
                              TMP=A(I,L)
                              A(I,L)=A(K,L)
  140                         A(K,L)=TMP
  150                     CONTINUE
  160                 CONTINUE
C    ..
C     COVARIANCE HAS BEEN COMPUTED AND REPERMUTED.
C     THE UPPER TRIANGULAR PART OF THE
C     SYMMETRIC MATRIX (A**T*A)**(-1) HAS
C     REPLACED THE UPPER TRIANGULAR PART OF
C     THE A ARRAY.
                      write (*,210)
                          DO 170 I=1,N
  170                     write (*,220) (I,J,A(I,J),J=I,N)
  180                 CONTINUE
      STOP
      END