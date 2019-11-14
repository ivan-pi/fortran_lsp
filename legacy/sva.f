      subroutine SVA(A,MDA,M,N,MDATA,B,SING,KPVEC,NAMES,ISCALE,D,WORK)
c
C  SINGULAR VALUE ANALYSIS.  COMPUTES THE SINGULAR VALUE
C  DECOMPOSITION OF THE MATRIX OF A LEAST SQUARES PROBLEM, AND
C  PRODUCES A PRINTED REPORT.
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------
c  This 1995 version differs from the original 1973 version by the
c  addition of the arguments KPVEC() and WORK(), and by allowing user to
c  choose the length of names in NAMES().
c  KPVEC() allows the user to exercise options regarding printing.
c  WORK() provides 2*N locations of work space.  Originally SING() was
c  required to have 3*N elements, of which the last 2*N were used for
c  work space.  Now SING() only needs N elements.
c     ------------------------------------------------------------------
c                         Subroutine Arguments
c     A(,)     [inout]  On entry, contains the M x N matrix of the least
c              squares problem to be analyzed.  This could be a matrix
c              obtained by preliminary orthogonal transformations
c              applied to the actual problem matrix which may have had
c              more rows (See MDATA below.)
c
c     MDA   [in]  First dimensioning parameter for A(,).  Require
c              MDA .ge. max(M, N).
c
c     M,N   [in]  No. of rows and columns, respectively, in the
c              matrix, A.  Either M > N or M .le. N is permitted.
c              Require M > 0 and N > 0.
c
c     MDATA [in]  No. of rows in actual least squares problem.
c              Generally MDATA .ge. M.  MDATA is used only in computing
c              statistics for the report and is not used as a loop
c              count or array dimension.
c
c     B()   [inout]  On entry, contains the right-side vector, b, of the
c              least squares problem.  This vector is of length, M.
c              On return, contains the vector, g = (U**t)*b, where U
c              comes from the singular value decomposition of A.  The
c              vector , g, is also of length M.
c
c     SING()   [out] On return, contains the singular values of A, in
c              descending order, in locations indexed 1 thru min(M,N).
c              If M < N, locations indexed from M+1 through N will be
c              set to zero.
c
c     KPVEC()  [integer array, in]  Array of integers to select print
c              options.  KPVEC(1) determines whether the rest of
c              the array is to be used or ignored.
c              If KPVEC(1) = 1, the contents of (KPVEC(I), I=2,4)
c              will be used to set internal variables as follows:
c                       PRBLK = KPVEC(2)
c                       UNIT  = KPVEC(3)
c                       WIDTH = KPVEC(4)
c              If KPVEC(1) = 0 default settings will be used.  The user
c              need not dimension KPVEC() greater than 1.  The subr will
c              set PRBLK = 111111, UNIT = -1, and WIDTH = 69.
c
c              The internal variables PRBLK, UNIT, and WIDTH are
c              interpreted as follows:
c
c              PRBLK    The decimal representation of PRBLK must be
c              representable as at most 6 digits, each being 0 or 1.
c              The decimal digits will be interpreted as independant
c              on/off flags for the 6 possible blocks of printed output.
c              Examples:  111111 selects all blocks, 0 suppresses all
c              printing,  101010 selects the 1st, 3rd, and 5th blocks,
c              etc.
c              The six blocks are:
c              1. Header, with size and scaling option parameters.
c              2. V-matrix.  Amount of output depends on M and N.
c              3. Singular values and related quantities.  Amount of
c                 output depends on N.
c              4. Listing of YNORM and RNORM and their logarithms.
c                 Amount of output depends on N.
c              5. Levenberg-Marquart analysis.
c              6. Candidate solutions.  Amount of output depends on
c                 M and N.
c
c              UNIT     Selects the output unit.  If UNIT .ge. 0,
c              UNIT will be used as the output unit number.
c              If UNIT = -1, output will be written to the "*" output
c              unit, i.e., the standard system output unit.
c              The calling program unit is responsible for opening
c              and/or closing the selected output unit if the host
c              system requires these actions.
c
c        WIDTH    Default value is 79.  Determines the width of
c        blocks 2, 3, and 6 of the output report.
c        Block 3 will use 95(+1) cols if WIDTH .ge. 95, and otherwise
c        69(+1) cols.
c        Blocks 2 and 6 are printed by subroutine MFEOUT.  These blocks
c        generally use at most WIDTH(+1) cols, but will use more if
c        the names are so long that more space is needed to print one
c        name and one numeric column.  The (+1)'s above are reminders
c        that in all cases there is one extra initial column for Fortran
c        "carriage control".  The carriage control character will always
c        be a blank.
c        Blocks 1, 4, and 5 have fixed widths of 63(+1), 66(+1) and
c        66(+1), respectively.
c
c     NAMES()  [in]  NAMES(j), for j = 1, ..., N, may contain a
c              name for the jth component of the solution
c              vector.  The declared length of the elements of the
c              NAMES() array is not specifically limited, but a
c              greater length reduces the space available for columns
c              of the matris to be printed.
c              If NAMES(1) contains only blank characters,
c              it will be assumed that no names have been provided,
c              and this subr will not access the NAMES() array beyond
c              the first element.
C
C     ISCALE   [in]  Set by the user to 1, 2, or 3 to select the column
c              scaling option.
C                  1   SUBR WILL USE IDENTITY SCALING AND IGNORE THE D()
C                      ARRAY.
C                  2   SUBR WILL SCALE NONZERO COLS TO HAVE UNIT EUCLID-
C                      EAN LENGTH AND WILL STORE RECIPROCAL LENGTHS OF
C                      ORIGINAL NONZERO COLS IN D().
C                  3   USER SUPPLIES COL SCALE FACTORS IN D(). SUBR
C                      WILL MULT COL J BY D(J) AND REMOVE THE SCALING
C                      FROM THE SOLN AT THE END.
c
c     D()   [ignored or out or in]  Usage of D() depends on ISCALE as
c              described above.  When used, its length must be
c              at least N.
c
c     WORK()   [scratch]   Work space of length at least 2*N.  Used
c              directly in this subr and also in _SVDRS.
c     ------------------------------------------------------------------
      integer I, IE, IPASS, ISCALE, J, K, KPVEC(4), M, MDA, MDATA
      integer MINMN, MINMN1, MPASS, N, NSOL
      integer PRBLK, UNIT, WIDTH
      double precision A(MDA,N), A1, A2, A3, A4, ALAMB, ALN10, B(M)
      double precision D(N), DEL, EL, EL2
      double precision ONE, PCOEF, RL, RNORM, RS
      double precision SB, SING(N), SL, TEN, THOU, TWENTY
      double precision WORK(2*N), YL, YNORM, YS, YSQ, ZERO
      character*(*) NAMES(N)
      logical BLK(6), NARROW, STAR
      parameter( ZERO = 0.0d0, ONE = 1.0d0)
      parameter( TEN = 10.0d0, TWENTY = 20.0d0, THOU = 1000.0d0)
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  220 format (1X/' INDEX  SING. VAL.    P COEF   ',
     *        '   RECIPROCAL    G COEF         G**2   ',
     *        '   CUMULATIVE  SCALED SQRT'/
     *    31x,'   SING. VAL.',26x,
     *        '  SUM of SQRS  of CUM.S.S.')
  221 format (1X/' INDEX  SING. VAL.    P COEF   ',
     *        '   RECIPROCAL    G COEF     SCALED SQRT'/
     *    31x,'   SING. VAL.',13x,'  of CUM.S.S.')
  222 format (1X/' INDEX  SING. VAL.    G COEF         G**2   ',
     *        '   CUMULATIVE  SCALED SQRT'/
     *    44x,'  SUM of SQRS  of CUM.S.S.')

  230 format (' ',4X,'0',64X,2g13.4)
  231 format (' ',4X,'0',51X, g13.4)
  232 format (' ',4X,'0',38X,2g13.4)

  240 format (' ',i5,g12.4,6g13.4)

  260 format (1X,' M  = ',I6,',   N  =',I4,',   MDATA  =',I8)
  270 format (1X/' Singular Value Analysis of the least squares',
     *        ' problem,  A*X = B,'/
     *        ' scaled as (A*D)*Y = B.')
  280 format (1X/' Scaling option No.',I2,'.  D is a diagonal',
     * ' matrix with the following diagonal elements..'/(5X,10E12.4))
  290 format (1X/' Scaling option No. 1.   D is the identity matrix.'/
     * 1X)
  300 format (1X/' INDEX',12X,'YNORM      RNORM',11X,
     * '      LOG10      LOG10'/
     * 45X,'      YNORM      RNORM'/1X)
  310 format (' ',I5,6X,2E11.3,11X,2F11.3)
  320 format (1X/
     *' Norms of solution and residual vectors for a range of values'/
     *' of the Levenberg-Marquardt parameter, LAMBDA.'//
     *        '      LAMBDA      YNORM      RNORM',
     *         '      LOG10      LOG10      LOG10'/
     *     34X,'     LAMBDA      YNORM      RNORM')
  330 format (1X, 3E11.3, 3F11.3)
c     ------------------------------------------------------------------
      IF (M.LE.0 .OR. N.LE.0) RETURN
      MINMN = min(M,N)
      MINMN1 = MINMN + 1
      if(KPVEC(1) .eq. 0) then
         PRBLK = 111111
         UNIT = -1
         WIDTH = 79
      else
         PRBLK = KPVEC(2)
         UNIT = KPVEC(3)
         WIDTH = KPVEC(4)
      endif
      STAR = UNIT .lt. 0
c                           Build logical array BLK() by testing
c                           decimal digits of PRBLK.
      do 20 I=6, 1, -1
         J = PRBLK/10
         BLK(I) = (PRBLK - 10*J) .gt. 0
         PRBLK = J
   20 continue
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                         Optionally print header and M, N, MDATA
      if(BLK(1)) then
         if(STAR) then
            write (*,270)
            write (*,260) M,N,MDATA
         else
            write (UNIT,270)
            write (UNIT,260) M,N,MDATA
         endif
      endif
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                                  Handle scaling as selected by ISCALE.
      if( ISCALE .eq. 1) then
         if(BLK(1)) then
            if(STAR) then
               write (*,290)
            else
               write (UNIT,290)
            endif
         endif
      else
C
C                                  Apply column scaling to A.
C
         DO 52 J = 1,N
            A1 = D(J)
            if( ISCALE .le. 2) then
               SB = ZERO
               DO 30 I = 1,M
   30             SB = SB + A(I,J)**2
               A1 = sqrt(SB)
               IF (A1.EQ.ZERO) A1 = ONE
               A1 = ONE/A1
               D(J) = A1
            endif
            DO 50 I = 1,M
               A(I,J) = A(I,J)*A1
   50       continue
   52    continue
         if(BLK(1)) then
            if(STAR) then
               write (*,280) ISCALE,(D(J),J = 1,N)
            else
               write (UNIT,280) ISCALE,(D(J),J = 1,N)
            endif
         endif
      endif
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C         Compute the Singular Value Decomposition of the scaled matrix.
C
      call SVDRS (A,MDA,M,N,B,M,1,SING,WORK)
c
c                                               Determine NSOL.
      NSOL = MINMN
      do 60 J = 1,MINMN
         if(SING(J) .eq. ZERO) then
            NSOL = J-1
            go to 65
         endif
   60 continue
   65 continue
C
c              The array B() contains the vector G.
C              Compute cumulative sums of squares of components of
C              G and store them in WORK(I), I = 1,...,MINMN+1
C
      SB = ZERO
      DO 70 I = MINMN1,M
         SB = SB + B(I)**2
   70 CONTINUE
      WORK(MINMN+1) = SB
      DO  75 J = MINMN, 1, -1
         SB = SB + B(J)**2
         WORK(J) = SB
   75 CONTINUE
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                                          PRINT THE V MATRIX.
C
      if(BLK(2)) CALL MFEOUT (A,MDA,N,N,NAMES,1, UNIT, WIDTH)
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                                REPLACE V BY D*V IN THE ARRAY A()
      if (ISCALE .gt.1) then
         do 82 I = 1,N
            do 80 J = 1,N
               A(I,J) = D(I)*A(I,J)
   80       continue
   82    continue
      endif
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(BLK(3)) then
c                     Print singular values and other summary results.
c
c     Output will be done using one of two layouts.  The narrow
c     layout uses 69 cols + 1 for carriage control, and makes two passes
c     through the computation.
c     The wide layout uses 95 cols + 1 for carriage control, and makes
c     only one pass through the computation.
c
C                       G  NOW IN  B() ARRAY.  V NOW IN A(,) ARRAY.
C
      NARROW = WIDTH .lt. 95
      MPASS = 1
      if(NARROW) MPASS = 2
      do 170 IPASS = 1, MPASS
         if(STAR) then
            if(NARROW) then
               if(IPASS .eq. 1) then
                  write(*,221)
               else
                  write(*,222)
               endif
            else
               write (*,220)
            endif
         else
            if(NARROW) then
               if(IPASS .eq. 1) then
                  write(UNIT,221)
               else
                  write(UNIT,222)
               endif
            else
               write (UNIT,220)
            endif
         endif
c                               The following stmt converts from
c                               integer to floating-point.
      A3 = WORK(1)
      A4 = sqrt(A3/ max(1,MDATA))
      if(STAR) then
         if(NARROW) then
            if(IPASS .eq. 1) then
               write(*,231) A4
            else
               write(*,232) A3, A4
            endif
         else
            write (*,230) A3,A4
         endif
      else
         if(NARROW) then
            if(IPASS .eq. 1) then
               write(UNIT,231) A4
            else
               write(UNIT,232) A3, A4
            endif
         else
            write (UNIT,230) A3,A4
         endif
      endif
C
      DO 160 K = 1,MINMN
         if (SING(K).EQ.ZERO) then
            PCOEF = ZERO
            if(STAR) then
               write (*,240) K,SING(K)
            else
               write (UNIT,240) K,SING(K)
            endif
         else
            PCOEF = B(K) / SING(K)
            A1 = ONE / SING(K)
            A2 = B(K)**2
            A3 = WORK(K+1)
            A4 = sqrt(A3/max(1,MDATA-K))
            if(STAR) then
               if(NARROW) then
                  if(IPASS .eq. 1) then
                     write(*,240) K,SING(K),PCOEF,A1,B(K),      A4
                  else
                     write(*,240) K,SING(K),         B(K),A2,A3,A4
                  endif
               else
                  write (*,240) K,SING(K),PCOEF,A1,B(K),A2,A3,A4
               endif
            else
               if(NARROW) then
                  if(IPASS .eq. 1) then
                     write(UNIT,240) K,SING(K),PCOEF,A1,B(K),      A4
                  else
                     write(UNIT,240) K,SING(K),         B(K),A2,A3,A4
                  endif
               else
                  write (UNIT,240) K,SING(K),PCOEF,A1,B(K),A2,A3,A4
               endif

            endif
         endif
  160 continue
  170 continue
      endif
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if( BLK(4) ) then
C
C         Compute and print values of YNORM, RNORM and their logarithms.
C
      if(STAR) then
         write (*,300)
      else
         write (UNIT,300)
      endif
      YSQ = ZERO
      do 180 J = 0, NSOL
         if(J .ne. 0) YSQ = YSQ + (B(J) / SING(J))**2
         YNORM = sqrt(YSQ)
         RNORM = sqrt(WORK(J+1))
         YL = -THOU
         IF (YNORM .GT. ZERO) YL = log10(YNORM)
         RL = -THOU
         IF (RNORM .GT. ZERO) RL = log10(RNORM)
         if(STAR) then
            write (*,310) J,YNORM,RNORM,YL,RL
         else
            write (UNIT,310) J,YNORM,RNORM,YL,RL
         endif
  180 continue
      endif
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if( BLK(5) .and. SING(1) .ne. ZERO ) then
C
C     COMPUTE VALUES OF XNORM AND RNORM FOR A SEQUENCE OF VALUES OF
C     THE LEVENBERG-MARQUARDT PARAMETER.
C
         EL = log10(SING(1)) + ONE
         EL2 = log10(SING(NSOL)) - ONE
         DEL = (EL2-EL) / TWENTY
         ALN10 = log(TEN)
         if(STAR) then
            write (*,320)
         else
            write (UNIT,320)
         endif
         DO 200 IE = 1,21
C                                    COMPUTE        ALAMB = 10.0**EL
            ALAMB = EXP(ALN10*EL)
            YS = ZERO
            RS = WORK(NSOL+1)
            DO 190 I = 1,MINMN
               SL = SING(I)**2 + ALAMB**2
               YS = YS + (B(I)*SING(I)/SL)**2
               RS = RS + (B(I)*(ALAMB**2)/SL)**2
  190       CONTINUE
            YNORM = sqrt(YS)
            RNORM = sqrt(RS)
            RL = -THOU
            IF (RNORM.GT.ZERO) RL = log10(RNORM)
            YL = -THOU
            IF (YNORM.GT.ZERO) YL = log10(YNORM)
            if(STAR) then
               write (*,330) ALAMB,YNORM,RNORM,EL,YL,RL
            else
               write (UNIT,330) ALAMB,YNORM,RNORM,EL,YL,RL
            endif
            EL = EL + DEL
  200    CONTINUE
      endif
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C                   Compute and optionally print candidate solutions.
C
      do 215 K = 1,NSOL
         PCOEF = B(K) / SING(K)
         DO 210 I = 1,N
                A(I,K) = A(I,K) * PCOEF
  210           IF (K.GT.1) A(I,K) = A(I,K) + A(I,K-1)
  215 continue
      if (BLK(6) .and. NSOL.GE.1)
     *       CALL MFEOUT (A,MDA,N,NSOL,NAMES,2,UNIT,WIDTH)
      return
      END