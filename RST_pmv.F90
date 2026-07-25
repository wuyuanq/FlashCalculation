
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_pmv

    use RST_globalFlashData
    use RST_fugacitycoef

    implicit none

contains

    function pmv(comp, P, local_z, Vold, K) result(local_v)

        integer, intent(in) :: comp
        real(kind=8), intent(in) :: P
        real(kind=8), dimension(:), pointer, intent(in) :: local_z
        real(kind=8), intent(in) :: Vold
        real(kind=8), dimension(:), pointer, intent(in out) :: K
        real(kind=8) :: local_v

        real(kind=8), dimension(:), pointer :: x, y
        real(kind=8) :: ZL, ZG
        real(kind=8), dimension(:), pointer :: am, bm
        real(kind=8) :: al, ag, bl, bg, bigAL, bigAG, bigBL, bigBG
        real(kind=8) :: xiL, xiG
        real(kind=8) :: rhoL, rhoG, CfL, CfG
        real(kind=8) :: Vnew
        real(kind=8) :: beta, lp, rp, hbeta, h0, h1, xt, yt
        real(kind=8) :: criteria, incmole, molenew
        real(kind=8), dimension(:), pointer :: znew
        real(kind=8) :: phasemole(2)
        integer :: i, local_t, nn1, nn2

        if(abs(local_z(comp)) < 1.D-12) then
            local_v = 0.D0
            return
        end if

        allocate(x(Nc))
        allocate(y(Nc))
        allocate(znew(Nc))
        allocate(am(Nc))
        allocate(bm(Nc))

        incmole = local_z(comp)*1.D-5 ! need adjust according to the case
        molenew = 1.D0+incmole

        do i = 1, Nc
            znew(i) = local_z(i)/molenew
        end do
        znew(comp) = (local_z(comp)+incmole)/molenew

        do local_t = 1, MAXTIME

            h0 = 0.D0
            do i = 1, Nc
                h0 = h0 + K(i)*znew(i)
            end do
            h0 = h0 - 1.D0
            if(h0 <= 0.D0) then
                do i = 1, Nc
                    x(i) = znew(i)
                end do
                yt = 0.D0
                do i = 1, Nc
                    yt = yt + K(i)*znew(i)
                end do
                do i = 1, Nc
                    y(i) = K(i)*znew(i)/yt
                end do
                go to 10
            end if

            h1 = 0.D0
            do i = 1, Nc
                h1 = h1 + znew(i)/K(i)
            end do
            h1 = 1.D0 - h1
            if(h1 >= 0.D0) then
                do i = 1, Nc
                    y(i) = znew(i)
                end do
                xt = 0.D0
                do i = 1, Nc
                    xt = xt + znew(i)/K(i)
                end do
                do i = 1, Nc
                    x(i) = znew(i)/K(i)/xt
                end do
                go to 10
            end if

            beta = 5.D-1
            lp = 0.D0
            rp = 1.D0
            do while(.true.)
                hbeta = 0.D0
                do i = 1, Nc
                    hbeta = hbeta + ((K(i)-1.D0)*znew(i))/(1.D0+beta*(K(i)-1.D0))
                end do
                if(abs(hbeta) < 1.D-12) then
                    exit
                end if

                if(hbeta > 0.D0) then
                    lp = beta
                else
                    rp = beta
                end if
                beta = (lp + rp)/2.D0

            end do

            do i = 1, Nc
                x(i) = znew(i)/(1.D0+(K(i)-1.D0)*beta)
                y(i) = K(i)*x(i)
            end do

10          call substeps( x, y, P, criteria, K, ZL, ZG, am, bm, al, ag, bl, bg, &!
                bigAL, bigAG, bigBL, bigBG, xiL, xiG, rhoL, rhoG, CfL, CfG )

            if(rhoL < rhoG) then
                call substeps( y, x, P, criteria, K, ZL, ZG, am, bm, al, ag, bl, bg, &!
                    bigAL, bigAG, bigBL, bigBG, xiL, xiG, rhoL, rhoG, CfL, CfG )
            end if

            if(criteria < 1.D-12) then
                do i = 1, Nc
                    if(znew(i) /= 0.D0) then
                        nn1 = i
                        exit
                    end if
                end do
                do i = nn1+1, Nc
                    if(znew(i) /= 0.D0) then
                        nn2 = i
                        exit
                    end if
                end do

                phasemole(1) = (y(nn2)*znew(nn1)-y(nn1)*znew(nn2))/(x(nn1)*y(nn2)-x(nn2)*y(nn1))*molenew
                phasemole(2) = (x(nn1)*znew(nn2)-x(nn2)*znew(nn1))/(x(nn1)*y(nn2)-x(nn2)*y(nn1))*molenew
                Vnew = phasemole(1)/xiL + phasemole(2)/xiG
                local_v = (Vnew-Vold)/incmole

                exit
            end if
        end do

        deallocate(x)
        deallocate(y)
        deallocate(znew)
        deallocate(am)
        deallocate(bm)

    end function pmv

    subroutine substeps(x, y, P, criteria, K, ZL, ZG, am, bm, al, ag, bl, bg, &!
        bigAL, bigAG, bigBL, bigBG, xiL, xiG, rhoL, rhoG, CfL, CfG)

        real(kind=8), dimension(:), pointer, intent(in) :: x
        real(kind=8), dimension(:), pointer, intent(in) :: y
        real(kind=8), intent(in) :: P
        real(kind=8), intent(out) :: criteria
        real(kind=8), dimension(:), pointer, intent(in out) :: K
        real(kind=8), intent(out) :: ZL, ZG
        real(kind=8), dimension(:), pointer, intent(in out) :: am, bm
        real(kind=8), intent(out) :: al, ag, bl, bg, bigAL, bigAG, bigBL, bigBG
        real(kind=8), intent(out) :: xiL
        real(kind=8), intent(out) :: xiG
        real(kind=8), intent(out) :: rhoL
        real(kind=8), intent(out) :: rhoG
        real(kind=8), intent(out) :: CfL, CfG

        real(kind=8), dimension(:), pointer :: phil, phig
        real(kind=8), dimension(:), pointer :: fl, fg
        integer :: i, nn

        allocate(phil(Nc))
        allocate(phig(Nc))

        call fugacitycoef(x, y, P, ZL, ZG, am, bm, al, ag, bl, bg, bigAL, &!
            bigAG, bigBL, bigBG, xiL, xiG, rhoL, rhoG, CfL, CfG, phil, phig)

        allocate(fl(Nc))
        allocate(fg(Nc))

        do i = 1, Nc
            fl(i) = x(i)*phil(i)*P
            fg(i) = y(i)*phig(i)*P
            K(i) = phil(i)/phig(i)
        end do

        criteria = 0.D0
        nn = 0 ! the number of existing components
        do i = 1, Nc
            if((x(i) /= 0.D0).or.(y(i) /= 0.D0)) then ! gurantee that component i exists
                criteria = criteria + (dlog(fg(i)/fl(i)*1.D0))**2.D0
                nn = nn + 1
            end if
        end do
        criteria = criteria / nn

        deallocate(fl)
        deallocate(fg)
        deallocate(phil)
        deallocate(phig)

    end subroutine substeps

end module RST_pmv
