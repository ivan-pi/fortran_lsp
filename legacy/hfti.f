      SUBROUTINE HFTI (A,MDA,M,N,B,MDB,NB,TAU,KRANK,RNORM,H,G,IP)
c
C  SOLVE LEAST SQUARES PROBLEM USING ALGORITHM, HFTI.
c  Householder Forward Triangulation with column Interchanges.
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------
      integer I, II, IP1, J, JB, JJ, K, KP1, KRANK
      integer L, LDIAG, LMAX, M, MDA, MDB, N, NB
c     integer IP(N)
c     double precision A(MDA,N),B(MDB,NB),H(N),G(N),RNORM(NB)
      integer IP(*)
      double precision A(MDA,*),B(MDB, *),H(*),G(*),RNORM( *)
      double precision DIFF, FACTOR, HMAX, SM, TAU, TMP, ZERO
      parameter(FACTOR = 0.001d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------
C
      K=0
      LDIAG=min(M,N)
      IF (LDIAG.LE.0) GO TO 270
          DO 80 J=1,LDIAG
          IF (J.EQ.1) GO TO 20
C
C     UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX
C    ..
          LMAX=J
              DO 10 L=J,N
              H(L)=H(L)-A(J-1,L)**2
              IF (H(L).GT.H(LMAX)) LMAX=L
   10         CONTINUE
          IF(DIFF(HMAX+FACTOR*H(LMAX),HMAX)) 20,20,50
C
C     COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX
C    ..
   20     LMAX=J
              DO 40 L=J,N
              H(L)=0.
                  DO 30 I=J,M
   30             H(L)=H(L)+A(I,L)**2
              IF (H(L).GT.H(LMAX)) LMAX=L
   40         CONTINUE
          HMAX=H(LMAX)
C    ..
C     LMAX HAS BEEN DETERMINED
C
C     DO COLUMN INTERCHANGES IF NEEDED.
C    ..
   50     CONTINUE
          IP(J)=LMAX
          IF (IP(J).EQ.J) GO TO 70
              DO 60 I=1,M
              TMP=A(I,J)
              A(I,J)=A(I,LMAX)
   60         A(I,LMAX)=TMP
          H(LMAX)=H(J)
C
C     COMPUTE THE J-TH TRANSFORMATION AND APPLY IT TO A AND B.
C    ..
   70     CALL H12 (1,J,J+1,M,A(1,J),1,H(J),A(1,J+1),1,MDA,N-J)
   80     CALL H12 (2,J,J+1,M,A(1,J),1,H(J),B,1,MDB,NB)
C
C     DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE, TAU.
C    ..
          DO 90 J=1,LDIAG
          IF (ABS(A(J,J)).LE.TAU) GO TO 100
   90     CONTINUE
      K=LDIAG
      GO TO 110
  100 K=J-1
  110 KP1=K+1
C
C     COMPUTE THE NORMS OF THE RESIDUAL VECTORS.
C
      IF (NB.LE.0) GO TO 140
          DO 130 JB=1,NB
          TMP=ZERO
          IF (KP1.GT.M) GO TO 130
              DO 120 I=KP1,M
  120         TMP=TMP+B(I,JB)**2
  130     RNORM(JB)=SQRT(TMP)
  140 CONTINUE
C                                           SPECIAL FOR PSEUDORANK = 0
      IF (K.GT.0) GO TO 160
      IF (NB.LE.0) GO TO 270
          DO 150 JB=1,NB
              DO 150 I=1,N
  150         B(I,JB)=ZERO
      GO TO 270
C
C     IF THE PSEUDORANK IS LESS THAN N COMPUTE HOUSEHOLDER
C     DECOMPOSITION OF FIRST K ROWS.
C    ..
  160 IF (K.EQ.N) GO TO 180
          DO 170 II=1,K
          I=KP1-II
  170     CALL H12 (1,I,KP1,N,A(I,1),MDA,G(I),A,MDA,1,I-1)
  180 CONTINUE
C
C
      IF (NB.LE.0) GO TO 270
          DO 260 JB=1,NB
C
C     SOLVE THE K BY K TRIANGULAR SYSTEM.
C    ..
              DO 210 L=1,K
              SM=ZERO
              I=KP1-L
              IF (I.EQ.K) GO TO 200
              IP1=I+1
                  DO 190 J=IP1,K
  190             SM=SM+A(I,J)*B(J,JB)
  200         continue
  210         B(I,JB)=(B(I,JB)-SM)/A(I,I)
C
C     COMPLETE COMPUTATION OF SOLUTION VECTOR.
C    ..
          IF (K.EQ.N) GO TO 240
              DO 220 J=KP1,N
  220         B(J,JB)=ZERO
              DO 230 I=1,K
  230         CALL H12 (2,I,KP1,N,A(I,1),MDA,G(I),B(1,JB),1,MDB,1)
C
C      RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE
C      COLUMN INTERCHANGES.
C    ..
  240         DO 250 JJ=1,LDIAG
              J=LDIAG+1-JJ
              IF (IP(J).EQ.J) GO TO 250
              L=IP(J)
              TMP=B(L,JB)
              B(L,JB)=B(J,JB)
              B(J,JB)=TMP
  250         CONTINUE
  260     CONTINUE
C    ..
C     THE SOLUTION VECTORS, X, ARE NOW
C     IN THE FIRST  N  ROWS OF THE ARRAY B(,).
C
  270 KRANK=K
      RETURN
      END