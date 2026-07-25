
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_PREOS

    use RST_mathlib
    use RST_globalFlashData

    implicit none

contains

    subroutine PREOS( x, y, P, ZL, ZG, am, bm, al, ag, bl, bg, &!
        bigAL, bigAG, bigBL, bigBG, xiL, xiG, rhoL, rhoG, CfL, CfG )

        real(kind=8), dimension(:), pointer, intent(in) :: x
        real(kind=8), dimension(:), pointer, intent(in) :: y
        real(kind=8), intent(in) :: P
        real(kind=8), intent(out) :: ZL, ZG
        real(kind=8), dimension(:), pointer, intent(in out) :: am, bm
        real(kind=8), intent(out) :: al, ag, bl, bg, bigAL, bigAG, bigBL, bigBG
        real(kind=8), intent(out) :: xiL
        real(kind=8), intent(out) :: xiG
        real(kind=8), intent(out) :: rhoL
        real(kind=8), intent(out) :: rhoG
        real(kind=8), intent(out) :: CfL, CfG

        real(kind=8) :: mm
        real(kind=8), dimension(:), pointer :: lambda, alpha
        integer :: i, j

        allocate(lambda(Nc))
        do i = 1, Nc
            lambda(i) = 3.7464D-1 + 1.5423D0*af(i) - 2.6992D-1*af(i)**2.D0
        end do

        allocate(alpha(Nc))
        do i = 1, Nc
            alpha(i) = (1.D0+lambda(i)*(1.D0-sqrt(Temp/ct(i))))**2.D0
        end do

        do i = 1, Nc
            am(i) = 4.5724D-1*alpha(i)*R**2.D0*ct(i)**2.D0/cp(i)
            bm(i) = 7.7796D-2*R*ct(i)/cp(i)
        end do

        deallocate(lambda)
        deallocate(alpha)

        al = 0.D0
        do i = 1, Nc
            do j = 1, Nc
                al = al + x(i)*x(j)*(1.D0-delta(i,j))*sqrt(am(i)*am(j))
            end do
        end do

        ag = 0.D0
        do i = 1, Nc
            do j = 1, Nc
                ag = ag + y(i)*y(j)*(1.D0-delta(i,j))*sqrt(am(i)*am(j))
            end do
        end do

        bl = 0.D0
        do i = 1, Nc
            bl = bl + x(i)*bm(i)
        end do

        bg = 0.D0
        do i = 1, Nc
            bg = bg + y(i)*bm(i)
        end do

        call computeLiquidPhase(x, P, al, bl, bigAL, bigBL, ZL, xiL, rhoL, CfL)

        call computeGasPhase(y, P, ag, bg, bigAG, bigBG, ZG, xiG, rhoG, CfG)

    end subroutine PREOS

    subroutine computeLiquidPhase(x, P, al, bl, bigAL, bigBL, ZL, xiL, rhoL, CfL)

        real(kind=8), dimension(:), pointer, intent(in) :: x
        real(kind=8), intent(in) :: P        
        real(kind=8), intent(in) :: al, bl
        real(kind=8), intent(out) :: bigAL, bigBL
        real(kind=8), intent(out) :: ZL
        real(kind=8), intent(out) :: xiL
        real(kind=8), intent(out) :: rhoL
        real(kind=8), intent(out) :: CfL

        real(kind=8), dimension(:), pointer :: c, C1toC2, ZRP, ZRM
        integer :: i, j, m, number, size
        real(kind=8) :: deri_AL_p, deri_BL_p, deri_ZL_p, deri_xiL_p_bs, deri_xiL_p
        real(kind=8) :: C2, C3, ctotal, xiL_bs, PP(4), local_Z(3)

        bigAL = al*P/(R*Temp)**2.D0
        bigBL = bl*P/(R*Temp)
        PP(1) = 1.D0
        PP(2) = -(1.D0-bigBL)
        PP(3) = bigAL-3.D0*bigBL**2.D0-2.D0*bigBL
        PP(4) = -(bigAL*bigBL-bigBL**2.D0-bigBL**3.D0)
        call getCubicRoot(PP, local_Z, size)

        allocate(ZRP(size))
        number = 0
        do i = 1, size
            if(local_Z(i)>0.D0) then
                ZRP(number+1) = local_Z(i)
                number = number + 1
            end if
        end do
        if(number == 0) then
            print *, 'There is no positive root in the cubic equation of liquid phase.', local_Z
            stop
        end if
        allocate(ZRM(number))
        j = 1
        do i = 1, number
            if(ZRP(i) > bigBL) then ! otherwise, the phil will become NaN
                ZRM(j) = ZRP(i)
                j = j + 1
            end if
        end do
        ZL = ZRM(1)
        do i = 2, j-1
            if(ZRM(i) < ZL) then
                ZL = ZRM(i)
            end if
        end do

        deallocate(ZRM)
        deallocate(ZRP)

        xiL_bs = P/(R*Temp*ZL)

        allocate(c(Nc))
        allocate(C1toC2(Nc))
        do i = 1, Nc
            C1toC2(i) = 1.1007D2*af(i)**4.D0 - 8.3807D1*af(i)**3.D0 + 1.8926D1* &!
                af(i)**2.D0 - 1.6348D0*af(i) - 6.6D-3
        end do
        C2 = 2.013645D-3
        C3 = 8.9D-1
        do i = 1, Nc
            c(i) = C1toC2(i)*C2 + C2*(Temp/ct(i)-C3)**2.D0
        end do
        ctotal = 0.D0
        do i = 1, Nc
            ctotal = ctotal + x(i)*c(i)*mw(i)
        end do
        xiL = 1.D0/(1.D0/xiL_bs + ctotal)
        deallocate(c)
        deallocate(C1toC2)

        rhoL = 0.D0
        do m = 1, Nc
            rhoL = rhoL + x(m)*mw(m)
        end do
        rhoL = rhoL*xiL

        deri_AL_p = al/(R*Temp)**2.D0
        deri_BL_p = bl/(R*Temp)
        deri_ZL_p = -(deri_BL_p*ZL**2.D0+(deri_AL_p-2.D0*(1.D0+3.D0*bigBL)*deri_BL_p)*ZL - &!
            (deri_AL_p*bigBL+(bigAL-2.D0*bigBL-3.D0*bigBL**2.D0)*deri_BL_p))/ &!
            (3.D0*ZL**2.D0-2.D0*(1.D0-bigBL)*ZL+(bigAL-2.D0*bigBL-3.D0*bigBL**2.D0))
        deri_xiL_p_bs = 1.D0/(R*Temp*ZL) - P/(R*Temp*ZL**2.D0)*deri_ZL_p
        deri_xiL_p = (1.D0/(1.D0+ctotal*xiL_bs)**2.D0)*deri_xiL_p_bs
        CfL = deri_xiL_p/xiL

    end subroutine computeLiquidPhase

    subroutine computeGasPhase(y, P, ag, bg, bigAG, bigBG, ZG, xiG, rhoG, CfG)

        real(kind=8), dimension(:), pointer, intent(in) :: y
        real(kind=8), intent(in) :: P
        real(kind=8), intent(in) :: ag, bg        
        real(kind=8), intent(out) :: bigAG, bigBG
        real(kind=8), intent(out) :: ZG
        real(kind=8), intent(out) :: xiG
        real(kind=8), intent(out) :: rhoG
        real(kind=8), intent(out) :: CfG

        integer :: i, m, number, size
        real(kind=8) :: deri_AG_p, deri_BG_p, deri_ZG_p
        real(kind=8) :: PP(4), local_Z(3)

        bigAG = ag*P/(R*Temp)**2.D0
        bigBG = bg*P/(R*Temp)

        PP(1) = 1.D0
        PP(2) = -(1.D0-bigBG)
        PP(3) = bigAG-3.D0*bigBG**2.D0-2.D0*bigBG
        PP(4) = -(bigAG*bigBG-bigBG**2.D0-bigBG**3.D0)
        call getCubicRoot(PP, local_Z, size)

        ZG = local_Z(1)
        do i = 2, size
            if(local_Z(i) > ZG) then
                ZG = local_Z(i)
            end if
        end do
        if(ZG < 0.D0) then
            print *, 'There is no positive root in the cubic equation of gas phase.', local_Z
            stop
        end if

        xiG = P/(R*Temp*ZG)

        rhoG = 0.D0
        do m = 1, Nc
            rhoG = rhoG + y(m)*mw(m)
        end do
        rhoG = rhoG*xiG

        deri_AG_p = ag/(R*Temp)**2.D0
        deri_BG_p = bg/(R*Temp)
        deri_ZG_p = -(deri_BG_p*ZG**2.D0+(deri_AG_p-2.D0*(1.D0+3.D0*bigBG)*deri_BG_p)*ZG - &!
            (deri_AG_p*bigBG+(bigAG-2.D0*bigBG-3.D0*bigBG**2.D0)*deri_BG_p))/ &!
            (3.D0*ZG**2.D0-2.D0*(1.D0-bigBG)*ZG+(bigAG-2.D0*bigBG-3.D0*bigBG**2.D0))
        CfG = 1.D0/P - 1.D0/ZG*deri_ZG_p

    end subroutine computeGasPhase

end module RST_PREOS
