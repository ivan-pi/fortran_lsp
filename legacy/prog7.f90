      program PROG7
!
!  DEMONSTRATE THE USE OF THE SUBROUTINE BVLS  FOR LEAST
!  SQUARES SOLVING WITH BOUNDS ON THE VARIABLES.
!
!  The original version of this code was developed by
!  Charles L. Lawson and Richard J. Hanson and published in the book
!  "SOLVING LEAST SQUARES PROBLEMS." REVISED APRIL, 1995 to accompany
!  reprinting of the book by SIAM.

IMPLICIT NONE

INTERFACE

  SUBROUTINE BVLS (A, B, BND, X, RNORM, NSETP, W, INDEX, IERR)
    REAL(KIND(1E0)) A(:,:), B(:), BND(:,:), X(:), RNORM, W(:)
    INTEGER NSETP, INDEX(:), IERR
  END SUBROUTINE

END INTERFACE

!     Test driver for BVLS.  Bounded Variables Least Squares.
!     C.L.Lawson, & R.J.Hanson, Jet Propulsion Laboratory, 1973 June 12
!     Changes made in September 1982, June 1986, October 1987.
!     Conversion made to Fortran 90 by R. J. Hanson, April 1995.
!     ------------------------------------------------------------------
!     Subprograms referenced: RANDOM_NUMBER, BVLS
!     ------------------------------------------------------------------

      integer, parameter :: MM=10, NN=10, MXCASE = 6, JSTEP=5
      INTEGER I, J, ICASE, IERR, m, n, nsetp, j1, j2
      integer     MTAB(MXCASE), NTAB(MXCASE)
      real(kind(1E0))   UNBTAB(MXCASE), BNDTAB(2,NN,MXCASE)
      real(kind(1E0))   A(MM,NN),B(MM),X(NN),W(NN)
      real(kind(1E0))   A2(MM,NN),B2(MM), R(MM), D(NN),BND(2,NN)
      real(kind(1E0)) RNORM, RNORM2, UNBND

      integer     INDEX(NN)
      data MTAB  / 2, 2, 4,  5, 10, 6 /
      data NTAB  / 2, 4, 2, 10,  5, 4 /
      data UNBTAB / 5 * 1.0E6,  999.0E0 /
      data ((BNDTAB(I,J,1),I=1,2),J=1,2)/ 1.,2.,    3.,4.  /
      data ((BNDTAB(I,J,2),I=1,2),J=1,4)/ 0,10,  0,10,  0,10,  0,10/
      data ((BNDTAB(I,J,3),I=1,2),J=1,2)/ 0,100,   -100,100/
      data ((BNDTAB(I,J,4),I=1,2),J=1,10)/&
             0,0,   -.3994E0,-.3994E0,  -1,1,     -.3E0,-.2E0,    21,22,&
                    -4,-3,  45,46,          100,101,  1.E6,1.E6,  -1,1/
      data ((BNDTAB(I,J,5),I=1,2),J=1,5)/&
                     0,1,  -1,0,  0,1,  .3E0,.4E0,  .048E0,.049E0/
      data ((BNDTAB(I,J,6),I=1,2),J=1,4)/&
               -100.,100.,  999.,999.,   999.,999.,   999.,999. /
!     ------------------------------------------------------------------
      write(*,1002)
      DO ICASE = 1,MXCASE
      M = MTAB(ICASE)
      N = NTAB(ICASE)
      UNBND = UNBTAB(ICASE)
      DO J = 1,N
         BND(1,J) = BNDTAB(1,J,ICASE)
         BND(2,J) = BNDTAB(2,J,ICASE)
      END DO
    where(bnd(1,1:N) == UNBND) bnd(1,1:n)=-huge(1e0)
    where(bnd(2,1:N) == UNBND) bnd(2,1:n)= huge(1e0)

      write(*,'(1x/////1x,a,i3/1x)') 'Case No.',ICASE
      write(*,'(1x,a,i5,a,i5,a,g17.5)')&
             'M =',M,',   N =',N,',   UNBND =',UNBND
      write(*,'(1X/'' Bounds ='')')
      DO J1 = 1, N, JSTEP
         J2 = MIN(J1 - 1 + JSTEP, N)
         write(*,*) ' '
         write(*,1001) (BND(1,J),J=J1,J2)
         write(*,1001) (BND(2,J),J=J1,J2)
      END DO
      call random_number (b(1:m))
       DO J=1,N
          call random_number (A(1:M,J))
      END DO

      write(*,'(1X/'' A(,) ='')')
      DO J1 = 1, N, JSTEP
         J2 = MIN(J1 - 1 + JSTEP, N)
         write(*,*) ' '
         DO I = 1,M
            write(*,1001) (A(I,J),J=J1,J2)
         END DO
      END DO
!
      B2(1:M)=B(1:M)
      A2(1:M,1:N)=A(1:M,1:N)

      write(*,'(1X/'' B() =''/1X)')
      write(*,1001) (B(I),I=1,M)
!
      call BVLS  (A2(1:M,1:N), B2,BND,X,RNORM,NSETP,W,INDEX, IERR)
!
       IF(IERR > 0) THEN
           WRITE(*,'(1X, "ABNORMAL ERROR FLAG, IERR = ")') IERR
            STOP
       END IF
      write(*,'(1x/1x,a,i4)') &
        'After BVLS:  No. of components not at constraints =',NSETP
      write(*,'(1X/'' Solution vector, X() =''/1X)')
      write(*,1001) (X(J),J=1,N)

      R(1:M)=B(1:M)-matmul(A(1:M,1:N),X(1:N))
!
      RNORM2=sqrt(dot_product(R(1:M),R(1:M)))
!
      D(1:N)=matmul(R(1:M),A(1:M,1:N))
      write(*,'(1X/'' R = B - A*X Computed by the driver:''/1X)')
      write(*,1001) (R(I),I=1,M)
      write(*,'(1x,a,g17.5)') 'RNORM2 computed by the driver =',RNORM2
      write(*,'(1x,a,g17.5)') 'RNORM computed by BVLS       = ',RNORM
!
      write(*,'(1X/'' W = (A**T)*R Computed by the driver:''/1X)')
      write(*,1001) (D(J),J=1,N)
      write(*,'(1X/'' Dual vector from BVLS, W() =''/1X)')
      write(*,1001) (W(J),J=1,N)
!
      END DO ! ICASE
!

 1001 FORMAT(1X,5G14.5)
 1002 FORMAT(&
       ' TESTBVLS.. Test driver for BVLS,',&
       ' Bounded Variables Least Squares.'&
       /'          If the algorithm succeeds the solution vector, X(),'/&
       '           and the dual vector, W(),'/&
       '           should be related as follows:'/&
       '        X(i) not at a bound       =>  W(i) = 0'/&
       '        X(i) at its lower bound  =>  W(i) .le. 0'/&
       '        X(i) at its upper bound  =>  W(i) .ge. 0'/&
       ' except that if an upper bound and lower bound are equal then'/&
       ' the corresponding X(i) must take that value and W(i) may have'/&
       ' any value.'/1X)

      end program PROG7