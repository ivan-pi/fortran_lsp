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
