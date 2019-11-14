      SUBROUTINE BNDACC (G,MDG,NB,IP,IR,MT,JT)  
c
C  SEQUENTIAL ALGORITHM FOR BANDED LEAST SQUARES PROBLEM..  
C  ACCUMULATION PHASE.      FOR SOLUTION PHASE USE BNDSOL.  
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C   
C     THE CALLING PROGRAM MUST SET IR=1 AND IP=1 BEFORE THE FIRST CALL  
C     TO BNDACC FOR A NEW CASE.     
C   
C     THE SECOND SUBSCRIPT OF G( ) MUST BE DIMENSIONED AT LEAST 
C     NB+1 IN THE CALLING PROGRAM.  
c     ------------------------------------------------------------------
      integer I, J, IE, IG, IG1, IG2, IP, IR, JG, JT, K, KH, L, LP1
      integer MDG, MH, MT, MU, NB, NBP1
c     double precision G(MDG,NB+1)
      double precision G(MDG,*)
      double precision RHO, ZERO
      parameter(ZERO = 0.0d0)
c     ------------------------------------------------------------------
C   
C              ALG. STEPS 1-4 ARE PERFORMED EXTERNAL TO THIS SUBROUTINE.
C   
      NBP1=NB+1 
      IF (MT.LE.0) RETURN   
C                                             ALG. STEP 5   
      IF (JT.EQ.IP) GO TO 70
C                                             ALG. STEPS 6-7
      IF (JT.LE.IR) GO TO 30
C                                             ALG. STEPS 8-9
      DO 10 I=1,MT  
        IG1=JT+MT-I     
        IG2=IR+MT-I     
        DO 10 J=1,NBP1  
   10   G(IG1,J)=G(IG2,J)   
C                                             ALG. STEP 10  
      IE=JT-IR  
      DO 20 I=1,IE  
        IG=IR+I-1   
        DO 20 J=1,NBP1  
   20   G(IG,J)=ZERO    
C                                             ALG. STEP 11  
      IR=JT 
C                                             ALG. STEP 12  
   30 MU=min(NB-1,IR-IP-1) 
      IF (MU.EQ.0) GO TO 60 
C                                             ALG. STEP 13  
      DO 50 L=1,MU  
C                                             ALG. STEP 14  
        K=min(L,JT-IP) 
C                                             ALG. STEP 15  
        LP1=L+1 
        IG=IP+L 
        DO 40 I=LP1,NB  
          JG=I-K
   40     G(IG,JG)=G(IG,I)  
C                                             ALG. STEP 16  
        DO 50 I=1,K     
        JG=NBP1-I   
   50   G(IG,JG)=ZERO   
C                                             ALG. STEP 17  
   60 IP=JT 
C                                             ALG. STEPS 18-19  
   70 MH=IR+MT-IP   
      KH=min(NBP1,MH)  
C                                             ALG. STEP 20  
      DO 80 I=1,KH  
        CALL H12 (1,I,max(I+1,IR-IP+1),MH,G(IP,I),1,RHO,   
     *            G(IP,I+1),1,MDG,NBP1-I)   
   80 continue
C                                             ALG. STEP 21  
      IR=IP+KH  
C                                             ALG. STEP 22  
      IF (KH.LT.NBP1) GO TO 100     
C                                             ALG. STEP 23  
      DO 90 I=1,NB  
   90   G(IR-1,I)=ZERO  
C                                             ALG. STEP 24  
  100 CONTINUE  
C                                             ALG. STEP 25  
      RETURN
      END   