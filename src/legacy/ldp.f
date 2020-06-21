      SUBROUTINE LDP (G,MDG,M,N,H,X,XNORM,W,INDEX,MODE)
c
C  Algorithm LDP: LEAST DISTANCE PROGRAMMING
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1974 MAR 1, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------
      integer I, IW, IWDUAL, IY, IZ, J, JF, M, MDG, MODE, N, NP1
c     integer INDEX(M)
c     double precision G(MDG,N), H(M), X(N), W(*)
      integer INDEX(*)
      double precision G(MDG,*), H(*), X(*), W(*)
      double precision DIFF, FAC, ONE, RNORM, XNORM, ZERO
      parameter(ONE = 1.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------
      IF (N.LE.0) GO TO 120
          DO 10 J=1,N
   10     X(J)=ZERO
      XNORM=ZERO
      IF (M.LE.0) GO TO 110
C
C     THE DECLARED DIMENSION OF W() MUST BE AT LEAST (N+1)*(M+2)+2*M.
C
C      FIRST (N+1)*M LOCS OF W()   =  MATRIX E FOR PROBLEM NNLS.
C       NEXT     N+1 LOCS OF W()   =  VECTOR F FOR PROBLEM NNLS.
C       NEXT     N+1 LOCS OF W()   =  VECTOR Z FOR PROBLEM NNLS.
C       NEXT       M LOCS OF W()   =  VECTOR Y FOR PROBLEM NNLS.
C       NEXT       M LOCS OF W()   =  VECTOR WDUAL FOR PROBLEM NNLS.
C     COPY G**T INTO FIRST N ROWS AND M COLUMNS OF E.
C     COPY H**T INTO ROW N+1 OF E.
C
      IW=0
          DO 30 J=1,M
              DO 20 I=1,N
              IW=IW+1
   20         W(IW)=G(J,I)
          IW=IW+1
   30     W(IW)=H(J)
      JF=IW+1
C                                STORE N ZEROS FOLLOWED BY A ONE INTO F.
          DO 40 I=1,N
          IW=IW+1
   40     W(IW)=ZERO
      W(IW+1)=ONE
C
      NP1=N+1
      IZ=IW+2
      IY=IZ+NP1
      IWDUAL=IY+M
C
      CALL NNLS (W,NP1,NP1,M,W(JF),W(IY),RNORM,W(IWDUAL),W(IZ),INDEX,
     *           MODE)
C                      USE THE FOLLOWING RETURN IF UNSUCCESSFUL IN NNLS.
      IF (MODE.NE.1) RETURN
      IF (RNORM) 130,130,50
   50 FAC=ONE
      IW=IY-1
          DO 60 I=1,M
          IW=IW+1
C                               HERE WE ARE USING THE SOLUTION VECTOR Y.
   60     FAC=FAC-H(I)*W(IW)
C
      IF (DIFF(ONE+FAC,ONE)) 130,130,70
   70 FAC=ONE/FAC
          DO 90 J=1,N
          IW=IY-1
              DO 80 I=1,M
              IW=IW+1
C                               HERE WE ARE USING THE SOLUTION VECTOR Y.
   80         X(J)=X(J)+G(I,J)*W(IW)
   90     X(J)=X(J)*FAC
          DO 100 J=1,N
  100     XNORM=XNORM+X(J)**2
      XNORM=sqrt(XNORM)
C                             SUCCESSFUL RETURN.
  110 MODE=1
      RETURN
C                             ERROR RETURN.       N .LE. 0.
  120 MODE=2
      RETURN
C                             RETURNING WITH CONSTRAINTS NOT COMPATIBLE.
  130 MODE=4
      RETURN
      END