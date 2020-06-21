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