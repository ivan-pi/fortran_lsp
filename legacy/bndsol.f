      SUBROUTINE BNDSOL (MODE,G,MDG,NB,IP,IR,X,N,RNORM)     
c
C  SEQUENTIAL SOLUTION OF A BANDED LEAST SQUARES PROBLEM..  
C  SOLUTION PHASE.   FOR THE ACCUMULATION PHASE USE BNDACC.     
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C   
C     MODE = 1     SOLVE R*X=Y   WHERE R AND Y ARE IN THE G( ) ARRAY    
C                  AND X WILL BE STORED IN THE X( ) ARRAY.  
C            2     SOLVE (R**T)*X=Y   WHERE R IS IN G( ),   
C                  Y IS INITIALLY IN X( ), AND X REPLACES Y IN X( ),    
C            3     SOLVE R*X=Y   WHERE R IS IN G( ).
C                  Y IS INITIALLY IN X( ), AND X REPLACES Y IN X( ).    
C   
C     THE SECOND SUBSCRIPT OF G( ) MUST BE DIMENSIONED AT LEAST 
C     NB+1 IN THE CALLING PROGRAM.  
      integer I, I1, I2, IE, II, IP, IR, IX, J, JG, L
      integer MDG, MODE, N, NB, NP1, IRM1
      double precision G(MDG,*), RNORM, RSQ, S, X(N), ZERO
      parameter(ZERO = 0.0d0)
C   
      RNORM=ZERO
      GO TO (10,90,50), MODE
C                                   ********************* MODE = 1  
C                                   ALG. STEP 26
   10      DO 20 J=1,N  
   20      X(J)=G(J,NB+1)   
      RSQ=ZERO  
      NP1=N+1   
      IRM1=IR-1 
      IF (NP1.GT.IRM1) GO TO 40     
           DO 30 J=NP1,IRM1 
   30      RSQ=RSQ+G(J,NB+1)**2     
      RNORM=SQRT(RSQ)   
   40 CONTINUE  
C                                   ********************* MODE = 3  
C                                   ALG. STEP 27
   50      DO 80 II=1,N 
           I=N+1-II     
C                                   ALG. STEP 28
           S=ZERO   
           L=max(0,I-IP)   
C                                   ALG. STEP 29
           IF (I.EQ.N) GO TO 70     
C                                   ALG. STEP 30
           IE=min(N+1-I,NB)
                DO 60 J=2,IE
                JG=J+L  
                IX=I-1+J
   60           S=S+G(I,JG)*X(IX)   
C                                   ALG. STEP 31
   70      continue
           IF (G(I,L+1) .eq. ZERO) go to 130
   80      X(I)=(X(I)-S)/G(I,L+1)   
C                                   ALG. STEP 32
      RETURN
C                                   ********************* MODE = 2  
   90      DO 120 J=1,N 
           S=ZERO   
           IF (J.EQ.1) GO TO 110    
           I1=max(1,J-NB+1)
           I2=J-1   
                DO 100 I=I1,I2  
                L=J-I+1+max(0,I-IP)
  100           S=S+X(I)*G(I,L)     
  110      L=max(0,J-IP)   
           IF (G(J,L+1) .eq. ZERO) go to 130
  120      X(J)=(X(J)-S)/G(J,L+1)   
      RETURN
C   
  130 write (*,'(/a/a,4i6)')' ZERO DIAGONAL TERM IN BNDSOL.',
     *      ' MODE,I,J,L = ',MODE,I,J,L  
      STOP  
      END