module ls_legacy_mod

    implicit none
    private

    public :: bndacc, bndsol
    public :: g1, g2
    public :: nnls, bvls
    public :: diff, gen
    public :: sva, svdrs, mfeout, ldp, qrbd
    public :: h12, hfti

    interface
!>
!   Sequential algorithm for banded least squares problem. Accumulation phase.
!
!   For the solution phase use `bndsol`.
!
        subroutine bndacc(g,mdg,nb,ip,ir,mt,jt)
            integer(kind=4), intent(in) :: mdg
            real(kind=8), intent(inout) :: g(mdg,*)
            integer(kind=4), intent(in) :: nb
            integer(kind=4), intent(inout) :: ip
            integer(kind=4), intent(inout) :: ir
            integer(kind=4), intent(in) :: mt
            integer(kind=4), intent(in) :: jt
        end subroutine bndacc
!>
!   Sequential solution of a banded least squares problem. Solution phase.
!
!   For the accumulation phase use `bndacc`.
!
        subroutine bndsol(mode,g,mdg,nb,ip,ir,x,n,rnorm)
            integer(kind=4), intent(in) :: mode
            integer(kind=4), intent(in) :: mdg
            integer(kind=4), intent(in) :: nb
            real(kind=8), intent(inout) :: g(mdg,*)
            integer(kind=4), intent(inout) :: ip
            integer(kind=4), intent(inout) :: ir
            integer(kind=4), intent(in) :: n
            real(kind=8), intent(inout) :: x(n)
            real(kind=8), intent(out) :: rnorm
        end subroutine bndsol
!>
!   Bounded variables least squares
!
        subroutine bvls(a,b,bnd,x,rnorm,nsetp,w,index,ierr)
            real(kind=8), intent(inout) :: a(:,:)
            real(kind=8), intent(inout) :: b(:)
            real(kind=8), intent(in) :: bnd(:,:)
            real(kind=8), intent(out) :: x(:)
            real(kind=8), intent(out) :: rnorm
            integer(kind=4), intent(out) :: nsetp
            real(kind=8), intent(out) :: w(:)
            integer(kind=4), intent(out) :: index(:)
            integer(kind=4), intent(out) :: ierr
        end subroutine bvls
!>
!   Function used in tests that depend on machine precision
!
        function diff(x,y)
            real(kind=8), intent(in) :: x
            real(kind=8), intent(in) :: y
            real(kind=8) :: diff
        end function diff
!>
!   Compute orthogonal rotation matrix
!
        subroutine g1(a,b,cterm,sterm,sig)
            real(kind=8), intent(in) :: a
            real(kind=8), intent(in) :: b
            real(kind=8), intent(out) :: cterm
            real(kind=8), intent(out) :: sterm
            real(kind=8), intent(out) :: sig
        end subroutine g1
!>
!   Apply the rotation computed by `g1` to (x,y)
!
        subroutine g2(cterm,sterm,x,y)
            real(kind=8), intent(in) :: cterm
            real(kind=8), intent(in) :: sterm
            real(kind=8), intent(inout) :: x
            real(kind=8), intent(inout) :: y
        end subroutine g2
!>
!   Generate numbers for construction of test cases
!
        function gen(anoise)
            real(kind=8), intent(in) :: anoise
            real(kind=8) :: gen
        end function gen
!>
!   Construction and/or application of a single Householder transformation
!
        subroutine h12(mode,lpivot,l1,m,u,iue,up,c,ice,icv,ncv)
            integer(kind=4), intent(in) :: mode
            integer(kind=4), intent(in) :: lpivot
            integer(kind=4), intent(in) :: l1
            integer(kind=4), intent(in) :: m
            integer(kind=4), intent(in) :: iue
            real(kind=8), intent(inout) :: up
            real(kind=8), intent(inout) :: u(iue,*)
            real(kind=8), intent(inout) :: c(*)
            integer(kind=4), intent(in) :: ice
            integer(kind=4), intent(in) :: icv
            integer(kind=4), intent(in) :: ncv
        end subroutine h12
!>
!   Solution of the least squares problem by Householder transformations
!
        subroutine hfti(a,mda,m,n,b,mdb,nb,tau,krank,rnorm,h,g,ip)
            integer(kind=4), intent(in) :: mda
            integer(kind=4), intent(in) :: m
            integer(kind=4), intent(in) :: n
            real(kind=8), intent(inout) :: a(mda,*)
            integer(kind=4), intent(in) :: mdb
            integer(kind=4), intent(in) :: nb
            real(kind=8), intent(inout) :: b(mdb,*)
            real(kind=8), intent(in) :: tau
            integer(kind=4), intent(out) :: krank
            real(kind=8), intent(out) :: rnorm(*)
            real(kind=8), intent(inout) :: h(*)
            real(kind=8), intent(inout) :: g(*)
            integer(kind=4), intent(inout) :: ip(*)
        end subroutine hfti
!>
!   Least distance programming
!
!   Computes the solution of the constrained least squares problem
!
        subroutine ldp(g,mdg,m,n,h,x,xnorm,w,index,mode)
            integer(kind=4), intent(in) :: mdg
            integer(kind=4), intent(in) :: m
            integer(kind=4), intent(in) :: n
            real(kind=8), intent(in) :: g(mdg,*)
            real(kind=8), intent(in) :: h(*)
            real(kind=8), intent(inout) :: x(*)
            real(kind=8), intent(out) :: xnorm
            real(kind=8), intent(inout) :: w(*)
            integer(kind=4), intent(inout) :: index(*)
            integer(kind=4), intent(out) :: mode
        end subroutine ldp
!>
!   Labeled matrix output for use with singular value analysis
!
        subroutine mfeout(a,mda,m,n,names,mode,unit,width)
            integer(kind=4), intent(in) :: mda
            integer(kind=4), intent(in) :: m
            integer(kind=4), intent(in) :: n
            real(kind=8), intent(in) :: a(mda,n)
            character(*), intent(in) :: names(m)
            integer(kind=4), intent(in) :: mode
            integer(kind=4), intent(in) :: unit
            integer(kind=4), intent(in) :: width
        end subroutine mfeout
!>
!   Nonnegative linear least squares
!
        subroutine nnls(a,mda,m,n,b,x,rnorm,w,zz,index,mode)
            integer(kind=4), intent(in) :: mda
            integer(kind=4), intent(in) :: m
            integer(kind=4), intent(in) :: n
            real(kind=8), intent(inout) :: a(mda,*)
            real(kind=8), intent(inout) :: b(*)
            real(kind=8), intent(out) :: x(*)
            real(kind=8), intent(out) :: rnorm
            real(kind=8), intent(inout):: w(*)
            real(kind=8), intent(inout) :: zz(*)
            integer(kind=4), intent(inout) :: index(*)
            integer(kind=4), intent(out) :: mode
        end subroutine nnls
!>
!   QR algorithm for singular value decomposition of a bidiagonal matrix
!
        subroutine qrbd(ipass,q,e,nn,v,mdv,nrv,c,mdc,ncc)
            integer(kind=4), intent(out) :: ipass
            integer(kind=4), intent(in) :: nn
            real(kind=8), intent(inout) :: q(*)
            real(kind=8), intent(inout) :: e(*)
            integer(kind=4), intent(in) :: mdv
            integer(kind=4), intent(in) :: nrv
            real(kind=8), intent(inout) :: v(mdv,*)
            integer(kind=4), intent(in) :: mdc
            integer(kind=4), intent(in) :: ncc
            real(kind=8), intent(inout) :: c(mdc,*)
        end subroutine qrbd
!>
!   Singular value analysis
!
!   Computes the singular value decomposition of the matrix
!   of a least squares problem, and produces a printed report.
!
        subroutine sva(a,mda,m,n,mdata,b,sing,kpvec,names,iscale,d,work)
            integer(kind=4), intent(in) :: mda
            integer(kind=4), intent(in) :: m
            integer(kind=4), intent(in) :: n
            real(kind=8), intent(inout) :: a(mda,n)
            integer(kind=4), intent(in) :: mdata
            real(kind=8), intent(inout) :: b(m)
            real(kind=8), intent(out) :: sing(n)
            integer(kind=4), intent(in) :: kpvec(4)
            character(*), intent(in) :: names(n)
            integer(kind=4), intent(in) :: iscale
            real(kind=8), intent(inout) :: d(n)
            real(kind=8), intent(inout) :: work(2*n)
        end subroutine sva
!>
!   Singular value decomposition also treating right side vector
!
        subroutine svdrs(a,mda,m1,n1,b,mdb,nb,s,work)
            integer(kind=4), intent(in) :: mda
            integer(kind=4), intent(in) :: m1
            integer(kind=4), intent(in) :: n1
            real(kind=8), intent(inout) :: a(mda,*)
            integer(kind=4), intent(in) :: mdb
            integer(kind=4), intent(in) :: nb
            real(kind=8), intent(inout) :: b(mdb,*)
            real(kind=8), intent(out) :: s(*)
            real(kind=8), intent(inout) :: work(n1,2)
        end subroutine svdrs
    end interface


end module

module least_squares_mod

    use ls_legacy_mod, only: nnls_ => nnls, bvls_ => bvls, ldp_ => ldp, &
        & hfti_ => hfti

    implicit none
    private

    public :: nnls, bvls, ldp, hfti

    interface hfti
        module procedure hfti_b  ! single rhs
        module procedure hfti_bb ! multiple rhs
    end interface

contains

    subroutine hfti_b(A,b,x,tau,krank,rnorm)
        real(8), intent(inout) :: A(:,:)
        real(8), intent(in) :: b(size(A,1))
        real(8), intent(out) :: x(size(A,2))
        real(8), intent(in) :: tau
        integer, intent(out), optional :: krank
        real(8), intent(out), optional :: rnorm

        real(8), allocatable :: h(:), g(:), x_(:)
        integer, allocatable :: ip(:)
        real(8) :: rnorm_(1)
        integer :: mda, m, n, krank_

        mda = size(A,1)
        m = mda
        n = size(A,2)

        allocate(h(n),g(n),ip(n),x_(max(m,n)))

        x_(1:m) = b ! copy data into x
        call hfti_(A,mda,m,n,x_,m,1,tau,krank_,rnorm_,h,g,ip)
        x = x(1:n) ! copy data out

        if (present(krank)) krank = krank_
        if (present(rnorm)) rnorm = rnorm_(1)
    end subroutine

    subroutine hfti_bb(A,b,x,tau,krank,rnorm)
        real(8), intent(inout) :: A(:,:)
        real(8), intent(in) :: b(:,:)
        real(8), intent(out) :: x(size(A,2),size(b,2))
        real(8), intent(in) :: tau
        integer, intent(out), optional :: krank
        real(8), intent(out), optional :: rnorm(size(b,2))

        real(8), allocatable :: h(:), g(:), b_(:,:)
        integer, allocatable :: ip(:)
        real(8) :: rnorm_(size(b,2))
        integer :: mda, m, n, mdb, nb, krank_

        mda = size(A,1)
        m = mda
        n = size(A,2)

        if (size(b,1) < m) error stop "[hfti] wrong dimension for dummy variable b"
        mdb = m ! only use the first m columns of b
        nb = size(b,2)

        allocate(h(n),g(n),ip(n),b_(max(m,n),nb))
        b_(1:mdb,:) = b(1:mdb,:)
        call hfti_(A,mda,m,n,b_(1:m,:),mdb,nb,tau,krank_,rnorm_,h,g,ip)
        x = b_(1:n,1:nb)
        if (present(krank)) krank = krank_
        if (present(rnorm)) rnorm = rnorm_
    end subroutine

    subroutine nnls(A,b,x,rnorm,ierr)
        real(8), intent(inout) :: A(:,:)
        real(8), intent(inout) :: b(size(A,dim=1))
        real(8), intent(out) :: x(size(A,dim=2))
        real(8), intent(out), optional :: rnorm
        integer, intent(out), optional :: ierr

        real(8), allocatable :: w(:), zz(:)
        integer, allocatable :: index(:)
        real(8) :: rnorm_
        integer :: mda, m, n, mode

        ! Retrieve dimensions
        mda = size(A,dim=1)
        m = mda
        n = size(A,dim=2)

        ! Allocate work space
        allocate(w(n),zz(m),index(n))

        ! Call legacy routine
        call nnls_(a,mda,m,n,b,x,rnorm_,w,zz,index,mode)

        ! Return optional variables
        if (present(rnorm)) rnorm = rnorm_
        if (present(ierr)) ierr = mode

    end subroutine

    subroutine bvls(A,b,bnd,x,rnorm,ierr)
        real(8), intent(inout) :: A(:,:)
        real(8), intent(inout) :: b(size(A,1))
        real(8), intent(in) :: bnd(2,size(A,2))
        real(8), intent(out) :: x(size(A,2))
        real(8), intent(out), optional :: rnorm
        integer, intent(out), optional :: ierr

        real(8) :: rnorm_
        integer :: ierr_, nsetp
        integer, allocatable :: index(:)
        real(8), allocatable :: w(:)

        allocate(index(size(A,2)),w(size(A,2)))

        call bvls_(A,b,bnd,x,rnorm_,nsetp,w,index,ierr_)

        if (present(rnorm)) rnorm = rnorm_
        if (present(ierr)) ierr = ierr_

    end subroutine

    subroutine ldp(G,h,x,xnorm,ierr)
        real(8), intent(in) :: G(:,:)
        real(8), intent(in) :: h(size(G,1))
        real(8), intent(out) :: x(size(G,2))
        real(8), intent(out), optional :: xnorm
        integer, intent(out), optional :: ierr

        real(8) :: xnorm_
        real(8), allocatable :: w(:)
        integer :: mdg, m, n, mode
        integer, allocatable :: index(:)

        mdg = size(G,1)
        m = mdg
        n = size(G,2)

        allocate(w((n+1)*(m+2)+2*m),index(m))

        call ldp_(G,mdg,m,n,h,x,xnorm_,w,index,mode)

        if (present(xnorm)) xnorm = xnorm_
        if (present(ierr)) ierr = mode
    end subroutine

end module

program test_nnls

    use least_squares_mod, only: nnls
    implicit none

    call test1
    call test2
    call test3
    call test4
    call test5
    call test6
    call test7
    call test8
    call test9

contains

    ! Trivial unit tests with cases taken from
    ! http://www.turkupetcentre.net/reports/tpcmod0020_app_a.pdf

    ! 4x2 problem, unconstrained solution positive
    subroutine test1()
        real(8) :: A(4,2), b(4), x(2), rnorm
        integer :: ierr

        A = reshape([1, 2, 3, 4, 1, 4, 9, 16],[4,2])
        b = [real(8) :: 0.6, 2.2, 4.8, 8.4]

        call nnls(A,b,x,rnorm=rnorm,ierr=ierr)
        print *, "x = ", x
        print *, "rnorm = ", rnorm
        print *, "ierr = ", ierr
    end subroutine

    ! 4x3 problem, unconstrained solution positive
    subroutine test2()
        real(8) :: A(4,3), b(4), x(3), rnorm
        integer :: ierr

        A = reshape([1, 2, 3, 4, 1, 4, 9, 16, 1, 8, 27, 64],[4,3])
        b = [real(8) :: 0.73, 3.24, 8.31, 16.72]

        call nnls(A,b,x,rnorm=rnorm,ierr=ierr)
        print *, "x = ", x
        print *, "rnorm = ", rnorm
        print *, "ierr = ", ierr
    end subroutine

    ! Simple 4x4 problem, unconstrained solution non-negative
    subroutine test3()
        real(8) :: A(4,4), b(4), x(4), rnorm
        integer :: ierr

        A = reshape([1, 2, 3, 4, 1, 4, 9, 16, 1, 8, 27, 64,1,16,81,256],[4,4])
        b = [0.73d0,3.24d0,8.31d0,16.72d0]

        call nnls(A,b,x,rnorm=rnorm,ierr=ierr)
        print *, "x = ", x
        print *, "rnorm = ", rnorm
        print *, "ierr = ", ierr
    end subroutine

    ! Simple 4x3 problem, unconstrained solution non-negative
    subroutine test4()
        real(8) :: A(4,3), b(4), x(3), rnorm
        integer :: ierr

        A = reshape([1, 2, 3, 4, 1, 4, 9, 16, 1, 8, 27, 64],[4,3])
        b = [0.23d0,1.24d0,3.81d0,8.72d0]

        call nnls(A,b,x,rnorm=rnorm,ierr=ierr)
        print *, "x = ", x
        print *, "rnorm = ", rnorm
        print *, "ierr = ", ierr
    end subroutine

    ! Simple 4x3 problem, unconstrained solution indefinite
    subroutine test5()
        real(8) :: A(4,3), b(4), x(3), rnorm
        integer :: ierr

        A = reshape([1, 2, 3, 4, 1, 4, 9, 16, 1, 8, 27, 64],[4,3])
        b = [0.13d0,0.84d0,2.91d0,7.12d0]

        call nnls(A,b,x,rnorm=rnorm,ierr=ierr)
        print *, "x = ", x
        print *, "rnorm = ", rnorm
        print *, "ierr = ", ierr
    end subroutine

    ! 200x100 random problem
    subroutine test6()
        real(8) :: A(200,10), b(200), x(10), rnorm
        integer :: ierr

        call random_number(A)
        call random_number(b)

        call nnls(A,b,x,rnorm=rnorm,ierr=ierr)
        print *, "x = ", x
        print *, "rnorm = ", rnorm
        print *, "ierr = ", ierr
    end subroutine

    ! 3x2 problem taken from https://github.com/mlapshin/nnls
    subroutine test7()
        real(8) :: A(3,2), b(3), x(2), rnorm
        integer :: ierr

        A = reshape([0.5, 0.3, 0.2, 0.2, 0.7, 0.8],[3,2])
        b = [0.1,0.1,0.7]

        call nnls(A,b,x,rnorm=rnorm,ierr=ierr)
        print *, "x = ", x
        print *, "rnorm = ", rnorm
        print *, "ierr = ", ierr
    end subroutine

    ! Tests 8 and 9 are taken from https://github.com/stefanopalmieri/lsqnonneg/blob/master/lsqnonneg.py

    subroutine test8
        real(8) :: A(6,7), b(7), x(6), rnorm
        real(8) :: tA(7,6)
        integer :: ierr

        A(1,:) = [ real(8) :: 1, 0, 0.035, 0.09,  0, 0.125, 0.875]
        A(2,:) = [ real(8) :: 1, 0, 0.01,  0.95,  0, 0.960, 0.040]
        A(3,:) = [ real(8) :: 1, 0, 0.35,  0.058, 0, 0.408, 0.592]
        A(4,:) = [ real(8) :: 1, 1, 0,         0, 0, 1.000, 0.000]
        A(5,:) = [ real(8) :: 1, .96, 0,  0, 0, 0.96, 0.04]
        A(6,:) = [ real(8) :: 1, .40, 0,  0.40, 0.14, 0.94, 0.06]

        b = [ real(8) :: 1530, 278, 92.3, 150.5, 7, 527.8, 1002.2]

        tA = transpose(A)
        call nnls(tA,b,x,rnorm=rnorm,ierr=ierr)
        print *, "x = ", x
        print *, "rnorm = ", rnorm
        print *, "ierr = ", ierr
    end subroutine

    subroutine test9
        real(8) :: A(4,3), b(4), x(3), rnorm
        integer :: ierr

        A(1,:) = [0.0372, 0.2869, 0.4]
        A(2,:) = [0.6861, 0.7071, 0.3]
        A(3,:) = [0.6233, 0.6245, 0.1]
        A(4,:) = [0.6344, 0.6170, 0.5]

        b = [0.8587, 0.1781, 0.0747, 0.8405]

        call nnls(A,b,x,rnorm=rnorm,ierr=ierr)
        print *, "x = ", x
        print *, "rnorm = ", rnorm
        print *, "ierr = ", ierr
    end subroutine

end program