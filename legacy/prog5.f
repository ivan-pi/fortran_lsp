      program PROG5
c
C  EXAMPLE OF THE USE OF SUBROUTINES BNDACC AND BNDSOL TO SOLVE
C  SEQUENTIALLY THE BANDED LEAST SQUARES PROBLEM THAT ARISES IN
C  SPLINE CURVE FITTING.
C
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------
      integer MDG, MXY
      parameter(MDG = 12, MXY = 12)
      integer I, IC, IG, IP, IR, J, JT, L, M, MT, NBAND, NBP, NC
      double precision B(10), C(MXY), COV(MXY), G(MDG,5), H
      double precision ONE, P1, P2, Q(4), R, RDUMMY, RNORM
      double precision SIGFAC, SIGSQ, U, X(MXY), Y(MXY), YFIT, ZERO
      parameter(ONE = 1.0d0, ZERO = 0.0d0)
      data Y /2.2d0, 4.0d0, 5.0d0, 4.6d0, 2.8d0, 2.7d0,
     *        3.8d0, 5.1d0, 6.1d0, 6.3d0, 5.0d0, 2.0d0/
C     ------------------------------------------------------------------
  160 format (/'   I',8X,'X',10X,'Y',6X,'YFIT',4X,'R=Y-YFIT/1X')
  170 format (1X,I3,4X,F6.0,4X,F6.2,4X,F6.2,4X,F8.4)
  180 format (/' C =',6F10.5/(4X,6F10.5))
  190 format (3(2X,2I4,g15.7))
  200 format (/' COVARIANCE MATRIX OF THE SPLINE COEFFICIENTS.')
  210 format (/' RNORM  =',g15.8)
  220 format (/' SIGFAC =',g15.8)
  230 format (/' PROG5.  EXECUTE A SEQUENCE OF CUBIC SPLINE FITS'/
     * 9x,'TO A DISCRETE SET OF DATA.')
  240 format (////' NEW CASE..'/
     * ' NUMBER OF BREAKPOINTS, INCLUDING ENDPOINTS, IS ',I5/
     * ' BREAKPOINTS =',6f10.5,/(14X,6f10.5))
C     ------------------------------------------------------------------
      write (*,230)
      NBAND=4
      M=MDG
C                                  SET ABCISSAS OF DATA
           DO 10 I=1,M
   10      X(I)=2*I
C
C     BEGIN LOOP THRU CASES USING INCREASING NOS OF BREAKPOINTS.
C
           DO 150 NBP=5,10
           NC=NBP+2
C                                  SET BREAKPOINTS
           B(1)=X(1)
           B(NBP)=X(M)
           H=(B(NBP)-B(1))/(NBP-1)
           IF (NBP.LE.2) GO TO 30
                DO 20 I=3,NBP
   20           B(I-1)=B(I-2)+H
   30      CONTINUE
           write (*,240) NBP,(B(I),I=1,NBP)
C
C     INITIALIZE IR AND IP BEFORE FIRST CALL TO BNDACC.
C
           IR=1
           IP=1
           I=1
           JT=1
   40      MT=0
   50      CONTINUE
           IF (X(I).GT.B(JT+1)) GO TO 60
C
C                        SET ROW  FOR ITH DATA POINT
C
           U=(X(I)-B(JT))/H
           IG=IR+MT
           G(IG,1)=P1(ONE - U)
           G(IG,2)=P2(ONE - U)
           G(IG,3)=P2(U)
           G(IG,4)=P1(U)
           G(IG,5)=Y(I)
           MT=MT+1
           IF (I.EQ.M) GO TO 60
           I=I+1
           GO TO 50
C
C                   SEND BLOCK OF DATA TO PROCESSOR
C
   60      CONTINUE
           CALL BNDACC (G,MDG,NBAND,IP,IR,MT,JT)
           IF (I.EQ.M) GO TO 70
           JT=JT+1
           GO TO 40
C                   COMPUTE SOLUTION C()
   70      CONTINUE
           CALL BNDSOL (1,G,MDG,NBAND,IP,IR,C,NC,RNORM)
C                  WRITE SOLUTION COEFFICIENTS
           write (*,180) (C(L),L=1,NC)
           write (*,210) RNORM
C
C              COMPUTE AND PRINT X,Y,YFIT,R=Y-YFIT
C
           write (*,160)
           JT=1
                DO 110 I=1,M
   80           IF (X(I).LE.B(JT+1)) GO TO 90
                JT=JT+1
                GO TO 80
C
   90           U=(X(I)-B(JT))/H
                Q(1)=P1(ONE - U)
                Q(2)=P2(ONE - U)
                Q(3)=P2(U)
                Q(4)=P1(U)
                YFIT=ZERO
                     DO 100 L=1,4
                     IC=JT-1+L
  100                YFIT=YFIT+C(IC)*Q(L)
                R=Y(I)-YFIT
                write (*,170) I,X(I),Y(I),YFIT,R
  110           CONTINUE
C
C     COMPUTE RESIDUAL VECTOR NORM.
C
           IF (M.LE.NC) GO TO 150
           SIGSQ=(RNORM**2)/(M-NC)
           SIGFAC=sqrt(SIGSQ)
           write (*,220) SIGFAC
           write (*,200)
C
C     COMPUTE AND PRINT COLS. OF COVARIANCE.
C
                DO 140 J=1,NC
                     DO 120 I=1,NC
  120                COV(I)=ZERO
                COV(J)= ONE
                CALL BNDSOL (2,G,MDG,NBAND,IP,IR,COV,NC,RDUMMY)
                CALL BNDSOL (3,G,MDG,NBAND,IP,IR,COV,NC,RDUMMY)
C
C    COMPUTE THE JTH COL. OF THE COVARIANCE MATRIX.
                     DO 130 I=1,NC
  130                COV(I)=COV(I)*SIGSQ
  140           write (*,190) (L,J,COV(L),L=1,NC)
  150      CONTINUE
      STOP
      END
c     ==================================================================
C
C        DEFINE FUNCTIONS P1 AND P2 TO BE USED IN CONSTRUCTING
C        CUBIC SPLINES OVER UNIFORMLY SPACED BREAKPOINTS.
C
c     ==================================================================
      double precision function P1(T)
      double precision T
      P1 = 0.25d0 * T**2 * T
      return
      end
c     ==================================================================
      double precision function P2(T)
      double precision T
      P2 = -(1.0d0 - T)**2 * (1.0d0 + T) * 0.75d0 + 1.0d0
      return
      end