
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_mathlib

    implicit none

    real(kind=8), parameter :: TwoPi = 8.D0*atan(1.D0)

contains

    ! Gaussian Elimination method to solve the linear equation array AX=b
    subroutine AGAUS(A, B, N, X, L, JS)

        real(kind=8) :: A(N,N), B(N)
        integer, intent(in) :: N
        real(kind=8), intent(out) :: X(N)
        integer :: L, JS(N)

        real(kind=8) :: TT
        real(kind=8) :: D
        integer :: i, j, k, IS

        L = 1
        do k = 1, N-1
            D = 0.D0
            do i = k, N
                do j = k, N
                    if (abs(A(i,j)).GT.D) then
                        D = abs(A(i,j))
                        JS(k) = j
                        IS = i
                    end if
                end do
            end do
            if (D+1.D0.EQ.1.D0) then
                L = 0
            else
                if (JS(k).NE.k) then
                    do i = 1, N
                        TT = A(i,k)
                        A(i,k) = A(i,JS(k))
                        A(i,JS(k)) = TT
                    end do
                end if
                if (IS.NE.k) then
                    do j = k, N
                        TT = A(k,j)
                        A(k,j) = A(IS,j)
                        A(IS,j) = TT
                    end do
                    TT = B(k)
                    B(k) = B(IS)
                    B(IS) = TT
                end if
            end if
            if (L.EQ.0) then
                print *, 'Error: In AGAUS, L=0.'
                stop
            end if
            do j = k+1, N
                A(k,j) = A(k,j)/A(k,k)
            end do
            B(k) = B(k)/A(k,k)
            do i = k+1, N
                do j = k+1, N
                    A(i,j) = A(i,j)-A(i,k)*A(k,j)
                end do
                B(i) = B(i)-A(i,k)*B(k)
            end do
        end do

        if (abs(A(N,N))+1.D0.EQ.1.D0) then
            L = 0
            print *, 'Error: In AGAUS, L=0.'
            stop
        end if
    
        X(N) = B(N)/A(N,N)
        do i = N-1, 1, -1
            TT = 0.D0
            do j = i+1, N
                TT = TT+A(i,j)*X(j)
            end do
            X(i) = B(i)-TT
        end do

        JS(N) = N
        do k = N, 1, -1
            if (JS(k).NE.k) then
                TT = X(k)
                X(k) = X(JS(k))
                X(JS(k)) = TT
            end if
        end do

    end subroutine AGAUS

    ! get the root of the cubic equation
    subroutine getCubicRoot(PP, ZZ, size)

        real(kind=8), intent(in) :: PP(4)
        real(kind=8), intent(in out) :: ZZ(3)
        integer, intent(out) :: size

        real(kind=8) :: a, b, c, d, Alph, Beta, Delt, R1, R2, tht

        ZZ = 0.D0
        a = PP(1)
        b = PP(2)/(3.D0*a)
        c = PP(3)/(6.D0*a)
        d = PP(4)/(2.D0*a)

        Alph = -b*b*b + 3.D0*b*c - d
        Beta =  b*b - 2.D0*c
        Delt = Alph*Alph-Beta*Beta*Beta

        if(Delt > 0.D0) then
            tht = Alph+sqrt(Delt); R1 = sign(abs(tht)**(1.D0/3.D0), tht*1.D0)
            tht = Alph-sqrt(Delt); R2 = sign(abs(tht)**(1.D0/3.D0), tht*1.D0)
            ZZ(1) = -b+R1+R2
            size = 1
        else if(Delt == 0.D0) then
            R1 = sign(abs(Alph)**(1.D0/3.D0), Alph*1.D0)
            if(R1 == 0.D0) then
                ZZ(1) = -b
                size = 1
            else
                ZZ(1) = -b+2.D0*R1
                ZZ(2) = -b-R1
                size = 2
            end if
        else if(Delt < 0.D0) then
            tht = acos(Alph/(sqrt(Beta)*Beta))
            ZZ(1)  = -b+2.D0*sqrt(Beta)*cos(tht/3.D0)
            ZZ(2)  = -b+2.D0*sqrt(Beta)*cos((tht+TwoPi)/3.D0)
            ZZ(3)  = -b+2.D0*sqrt(Beta)*cos((tht-TwoPi)/3.D0)
            size = 3
        end if

    end subroutine getCubicRoot

    !Sort 'matrix' into ascending numerical order of the last column by straight insertion method.
    subroutine piksrt(row, column, matrix)

        integer :: row, column
        real(kind=8), dimension(:,:), pointer :: matrix

        integer :: i, j
        real(kind=8), dimension(:), pointer :: array
        logical :: jump

        allocate(array(column))

        do j = 2, row
            array = matrix(j,:)
            jump = .false.
            do i = j-1, 1, -1
                if(matrix(i,column) <= array(column)) then
                    matrix(i+1,:) = array
                    jump = .true.
                    exit
                end if
                matrix(i+1,:) = matrix(i,:)
            end do
            if(.not.jump) then
                i = 0
                matrix(i+1,:) = array
            end if
        end do

        deallocate(array)

    end subroutine piksrt

    function isPrime( num ) result(isP)

        integer, intent(in) :: num
        logical :: isP
        integer :: i

        do i = 2, floor(sqrt(num*1.D0))
            if(mod(num,i) == 0) then
                isP = .false.
                return
            end if
        end do

        isP = .true.

    end function isPrime

    subroutine genPrime(primesize)

        integer, intent(in) :: primesize !1.D4
        character(len=20) :: fprimetxt
        integer :: num, ierr, i

        fprimetxt = "../sg_hpc/prime.txt"
        open(unit=10, file=trim(adjustl(fprimetxt)), status='replace', iostat=ierr)
        if(ierr /= 0) then
            print *, 'open file error. ', ierr
            stop
        end if

        i = 1
        num = 2
        do while (i <= primesize)
            if(isPrime(num)) then
                write(10, fmt="(i6)") num
                i = i + 1
            end if
            num = num + 1
        end do

        close(10)

    end subroutine genPrime

end module RST_mathlib
