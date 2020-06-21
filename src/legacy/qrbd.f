C     SUBROUTINE QRBD (IPASS,Q,E,NN,V,MDV,NRV,C,MDC,NCC)
c
C  QR ALGORITHM FOR SINGULAR VALUES OF A BIDIAGONAL MATRIX.
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------
C     THE BIDIAGONAL MATRIX
C
C                       (Q1,E2,0...    )
C                       (   Q2,E3,0... )
C                D=     (       .      )
C                       (         .   0)
C                       (           .EN)
C                       (          0,QN)
C
C                 IS PRE AND POST MULTIPLIED BY
C                 ELEMENTARY ROTATION MATRICES
C                 RI AND PI SO THAT
C
C                 RK...R1*D*P1**(T)...PK**(T) = DIAG(S1,...,SN)
C
C                 TO WITHIN WORKING ACCURACY.
C
C  1. EI AND QI OCCUPY E(I) AND Q(I) AS INPUT.
C
C  2. RM...R1*C REPLACES 'C' IN STORAGE AS OUTPUT.
C
C  3. V*P1**(T)...PM**(T) REPLACES 'V' IN STORAGE AS OUTPUT.
C
C  4. SI OCCUPIES Q(I) AS OUTPUT.
C
C  5. THE SI'S ARE NONINCREASING AND NONNEGATIVE.
C
C     THIS CODE IS BASED ON THE PAPER AND 'ALGOL' CODE..
C REF..
C  1. REINSCH,C.H. AND GOLUB,G.H. 'SINGULAR VALUE DECOMPOSITION
C     AND LEAST SQUARES SOLUTIONS' (NUMER. MATH.), VOL. 14,(1970).
C
C     ------------------------------------------------------------------
      SUBROUTINE QRBD (IPASS,Q,E,NN,V,MDV,NRV,C,MDC,NCC)
C     ------------------------------------------------------------------
      integer MDC, MDV, NCC, NN, NRV
c     double precision C(MDC,NCC), E(NN), Q(NN),V(MDV,NN)
      double precision C(MDC,*  ), E(* ), Q(* ),V(MDV,* )
      integer I, II, IPASS, J, K, KK, L, LL, LP1, N, N10, NQRS
      double precision CS, DIFF, DNORM, F, G, H, SMALL
      double precision ONE, SN, T, TEMP, TWO, X, Y, Z, ZERO

      logical WNTV ,HAVERS,FAIL
      parameter(ONE = 1.0d0, TWO = 2.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------
      N=NN
      IPASS=1
      IF (N.LE.0) RETURN
      N10=10*N
      WNTV=NRV.GT.0
      HAVERS=NCC.GT.0
      FAIL=.FALSE.
      NQRS=0
      E(1)=ZERO
      DNORM=ZERO
           DO 10 J=1,N
   10      DNORM=max(abs(Q(J))+abs(E(J)),DNORM)
           DO 200 KK=1,N
           K=N+1-KK
C
C     TEST FOR SPLITTING OR RANK DEFICIENCIES..
C         FIRST MAKE TEST FOR LAST DIAGONAL TERM, Q(K), BEING SMALL.
   20       IF(K.EQ.1) GO TO 50
            IF(DIFF(DNORM+Q(K),DNORM) .ne. ZERO) go to 50
C
C     SINCE Q(K) IS SMALL WE WILL MAKE A SPECIAL PASS TO
C     TRANSFORM E(K) TO ZERO.
C
           CS=ZERO
           SN=-ONE
                DO 40 II=2,K
                I=K+1-II
                F=-SN*E(I+1)
                E(I+1)=CS*E(I+1)
                CALL G1 (Q(I),F,CS,SN,Q(I))
C         TRANSFORMATION CONSTRUCTED TO ZERO POSITION (I,K).
C
                IF (.NOT.WNTV) GO TO 40
                     DO 30 J=1,NRV
c
c                          Apply procedure G2 (CS,SN,V(J,I),V(J,K))
c
                        TEMP = V(J,I)
                        V(J,I) = CS*TEMP + SN*V(J,K)
                        V(J,K) =-SN*TEMP + CS*V(J,K)
   30                continue
C              ACCUMULATE RT. TRANSFORMATIONS IN V.
C
   40           CONTINUE
C
C         THE MATRIX IS NOW BIDIAGONAL, AND OF LOWER ORDER
C         SINCE E(K) .EQ. ZERO..
C
   50           DO 60 LL=1,K
                  L=K+1-LL
                  IF(DIFF(DNORM+E(L),DNORM) .eq. ZERO) go to 100
                  IF(DIFF(DNORM+Q(L-1),DNORM) .eq. ZERO) go to 70
   60           CONTINUE
C     THIS LOOP CAN'T COMPLETE SINCE E(1) = ZERO.
C
           GO TO 100
C
C         CANCELLATION OF E(L), L.GT.1.
   70      CS=ZERO
           SN=-ONE
                DO 90 I=L,K
                F=-SN*E(I)
                E(I)=CS*E(I)
                IF(DIFF(DNORM+F,DNORM) .eq. ZERO) go to 100
                CALL G1 (Q(I),F,CS,SN,Q(I))
                IF (HAVERS) then
                     DO 80 J=1,NCC
c
c                          Apply procedure G2 ( CS, SN, C(I,J), C(L-1,J)
c
                        TEMP = C(I,J)
                        C(I,J)   = CS*TEMP + SN*C(L-1,J)
                        C(L-1,J) =-SN*TEMP + CS*C(L-1,J)
   80                continue
                endif
   90           CONTINUE
C
C         TEST FOR CONVERGENCE..
  100      Z=Q(K)
           IF (L.EQ.K) GO TO 170
C
C         SHIFT FROM BOTTOM 2 BY 2 MINOR OF B**(T)*B.
           X=Q(L)
           Y=Q(K-1)
           G=E(K-1)
           H=E(K)
           F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(TWO*H*Y)
           G=sqrt(ONE+F**2)
           IF (F .ge. ZERO) then
              T=F+G
           else
              T=F-G
           endif
           F=((X-Z)*(X+Z)+H*(Y/T-H))/X
C
C         NEXT QR SWEEP..
           CS=ONE
           SN=ONE
           LP1=L+1
                DO 160 I=LP1,K
                G=E(I)
                Y=Q(I)
                H=SN*G
                G=CS*G
                CALL G1 (F,H,CS,SN,E(I-1))
                F=X*CS+G*SN
                G=-X*SN+G*CS
                H=Y*SN
                Y=Y*CS
                IF (WNTV) then
C
C              ACCUMULATE ROTATIONS (FROM THE RIGHT) IN 'V'
c
                     DO 130 J=1,NRV
c
c                          Apply procedure G2 (CS,SN,V(J,I-1),V(J,I))
c
                        TEMP = V(J,I-1)
                        V(J,I-1) = CS*TEMP + SN*V(J,I)
                        V(J,I)   =-SN*TEMP + CS*V(J,I)
  130                continue
                endif
                CALL G1 (F,H,CS,SN,Q(I-1))
                F=CS*G+SN*Y
                X=-SN*G+CS*Y
                IF (HAVERS) then
                     DO 150 J=1,NCC
c
c                          Apply procedure G2 (CS,SN,C(I-1,J),C(I,J))
c
                        TEMP = C(I-1,J)
                        C(I-1,J) = CS*TEMP + SN*C(I,J)
                        C(I,J)   =-SN*TEMP + CS*C(I,J)
  150                continue
                endif
c
C              APPLY ROTATIONS FROM THE LEFT TO
C              RIGHT HAND SIDES IN 'C'..
C
  160           CONTINUE
           E(L)=ZERO
           E(K)=F
           Q(K)=X
           NQRS=NQRS+1
           IF (NQRS.LE.N10) GO TO 20
C          RETURN TO 'TEST FOR SPLITTING'.
C
           SMALL=ABS(E(K))
           I=K
C          IF FAILURE TO CONVERGE SET SMALLEST MAGNITUDE
C          TERM IN OFF-DIAGONAL TO ZERO.  CONTINUE ON.
C      ..
                DO 165 J=L,K
                TEMP=ABS(E(J))
                IF(TEMP .EQ. ZERO) GO TO 165
                IF(TEMP .LT. SMALL) THEN
                     SMALL=TEMP
                     I=J
                end if
  165           CONTINUE
           E(I)=ZERO
           NQRS=0
           FAIL=.TRUE.
           GO TO 20
C     ..
C     CUTOFF FOR CONVERGENCE FAILURE. 'NQRS' WILL BE 2*N USUALLY.
  170      IF (Z.GE.ZERO) GO TO 190
           Q(K)=-Z
           IF (WNTV) then
                DO 180 J=1,NRV
  180           V(J,K)=-V(J,K)
           endif
  190      CONTINUE
C         CONVERGENCE. Q(K) IS MADE NONNEGATIVE..
C
  200      CONTINUE
      IF (N.EQ.1) RETURN
           DO 210 I=2,N
           IF (Q(I).GT.Q(I-1)) GO TO 220
  210      CONTINUE
      IF (FAIL) IPASS=2
      RETURN
C     ..
C     EVERY SINGULAR VALUE IS IN ORDER..
  220      DO 270 I=2,N
           T=Q(I-1)
           K=I-1
                DO 230 J=I,N
                IF (T.GE.Q(J)) GO TO 230
                T=Q(J)
                K=J
  230           CONTINUE
           IF (K.EQ.I-1) GO TO 270
           Q(K)=Q(I-1)
           Q(I-1)=T
           IF (HAVERS) then
                DO 240 J=1,NCC
                T=C(I-1,J)
                C(I-1,J)=C(K,J)
  240           C(K,J)=T
           endif

  250      IF (WNTV) then
                DO 260 J=1,NRV
                T=V(J,I-1)
                V(J,I-1)=V(J,K)
  260           V(J,K)=T
           endif
  270      CONTINUE
C         END OF ORDERING ALGORITHM.
C
      IF (FAIL) IPASS=2
      RETURN
      END