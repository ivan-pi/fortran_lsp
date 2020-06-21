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