
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2022-9-17 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_flashcalculation_sparsegrid

    use RST_globalFlashData
    use RST_hashTable
    use RST_PREOS
    use RST_pmv
    use RST_comprefac

    implicit none

contains

    subroutine flashcalculation_sparsegrid(P, local_z, x, y, xiL, xiG, rhoL, rhoG, sL, local_v, local_Cf, isW, isN)

        real(kind=8), intent(in) :: P
        real(kind=8), dimension(:), pointer, intent(in) :: local_z
        real(kind=8), dimension(:), pointer, intent(in out) :: x
        real(kind=8), dimension(:), pointer, intent(in out) :: y
        real(kind=8), intent(out) :: xiL
        real(kind=8), intent(out) :: xiG
        real(kind=8), intent(out) :: rhoL
        real(kind=8), intent(out) :: rhoG
        real(kind=8), intent(out) :: sL
        real(kind=8), dimension(:), pointer, intent(in out) :: local_v
        real(kind=8), intent(out) :: local_Cf
        logical, intent(out) :: isW, isN

        integer :: purenum
        real(kind=8), dimension(:), pointer :: bigK
        real(kind=8) :: CfL, CfG
        logical :: isPure
        real(kind=8) :: ZL, ZG
        real(kind=8), dimension(:), pointer :: am, bm
        real(kind=8) :: al, ag, bl, bg, bigAL, bigAG, bigBL, bigBG
        real(kind=8) :: phasemole(2)
        real(kind=8), dimension(:), pointer :: psat
        real(kind=8) :: psatmix
        real(kind=8) :: molevg, molevl
        real(kind=8), dimension(:), pointer :: tempsum
        integer :: nn1, nn2
        real(kind=8) :: sum, sumx, sumy
        integer :: i, j, m

        ! decide whether it is pure substance
        isPure = .false.
        purenum = 0
        do m = 1, Nc
            if(abs(local_z(m)-1.D0) < EQUALPREC) then
                isPure = .true.
                purenum = m
                exit
            end if
        end do

        ! if it is pure substance
        if(isPure) then

            allocate(am(Nc))
            allocate(bm(Nc))
            call PREOS( local_z, local_z, P, ZL, ZG, am, bm, al, ag, bl, &!
                bg, bigAL, bigAG, bigBL, bigBG, xiL, xiG, rhoL, rhoG, CfL, CfG )
            deallocate(am)
            deallocate(bm)

            if(Temp > ct(purenum)) then
                x = 0.D0
                y = local_z
                xiL = 0.D0
                rhoL = 0.D0
                local_Cf = CfG
                local_v = 0.D0
                local_v(purenum) = 1.D0/xiG
                sL = 0.D0
                isW = .false.
                isN = .true.
            else
                ! Antoine Equation
                psatmix = dexp((psatA(purenum)-psatB(purenum)/((Temp-2.73D2)+psatC(purenum)))*1.D0)
                psatmix = psatmix*1.33322D2 !convert from mmHg to Pa
                if(P > psatmix) then
                    x = local_z
                    y = 0.D0
                    xiG = 0.D0
                    rhoG = 0.D0
                    local_Cf = CfL
                    local_v = 0.D0
                    local_v(purenum) = 1.D0/xiL
                    sL = 1.D0
                    isN = .false.
                    isW = .true.
                else
                    x = 0.D0
                    y = local_z
                    xiL = 0.D0
                    rhoL = 0.D0
                    local_Cf = CfG
                    local_v = 0.D0
                    local_v(purenum) = 1.D0/xiG
                    sL = 0.D0
                    isW = .false.
                    isN = .true.
                end if
            end if

            return
        end if

        if(isFirstSG) then
            isFirstSG = .false.
            call retrieveHashTable()
        end if

        ! It is noted that only x and y are used in the following code.
        call sgInterpo(P, local_z, x, y, xiL, xiG, rhoL, rhoG, sL, local_v, local_Cf)

        isW = .true.
        isN = .true.
        sumx = 0.D0
        sumy = 0.D0
        do m = 1, Nc
           sumx = sumx + x(m)
           sumy = sumy + y(m)
        end do
        if(sumx < EQUALPREC) then
            isW = .false.
        end if
        if(sumy < EQUALPREC) then
            isN = .false.
        end if
        if(.not.isW.and..not.isN) then
            if(sumx < sumy) then
                isN = .true.
            else
                isW = .true.
            end if
        end if

        allocate(am(Nc))
        allocate(bm(Nc))

        if((isW.and..not.isN).or.(isN.and..not.isW)) then

            call PREOS(local_z, local_z, P, ZL, ZG, am, bm, al, ag, bl, bg, bigAL, &!
                bigAG, bigBL, bigBG, xiL, xiG, rhoL, rhoG, CfL, CfG)

            allocate(tempsum(Nc))
            do i = 1, Nc
                tempsum(i) = 0.D0
                do j = 1, Nc
                    tempsum(i) = tempsum(i) + local_z(j)*(1.D0-delta(i,j))*dsqrt(am(i)*am(j)*1.D0)
                end do
                tempsum(i) = 2.D0*tempsum(i)
            end do

            if(isW.and..not.isN) then
                molevl = ZL*R*Temp/P
                do m = 1, Nc
                    local_v(m) = (R*Temp/(molevl-bl)*(1.D0+bm(m)/(molevl-bl))-(tempsum(m)-(2.D0*al*bm(m)* &!
                        (molevl-bl))/(molevl*(molevl+bl)+bl*(molevl-bl)))/(molevl*(molevl+bl)+bl* &!
                        (molevl-bl)))/(R*Temp/(molevl-bl)**2.D0-2.D0*al*(molevl+bl)/(molevl*(molevl+bl)+ &!
                        bl*(molevl-bl))**2.D0)
                    if(local_z(m) < EQUALPREC) then
                        local_v(m) = 0.D0
                    end if
                end do
                local_Cf = CfL
                x = local_z
                y = 0.D0
                xiG = 0.D0
                rhoG = 0.D0
                sL = 1.D0
            elseif(isN.and..not.isW) then
                molevg = ZG*R*Temp/P
                do m = 1, Nc
                    local_v(m) = (R*Temp/(molevg-bg)*(1.D0+bm(m)/(molevg-bg))-(tempsum(m)-(2.D0*ag*bm(m)*(molevg-bg)) &!
                        /(molevg*(molevg+bg)+bg*(molevg-bg)))/(molevg*(molevg+bg)+bg*(molevg-bg)))/(R*Temp &!
                        /(molevg-bg)**2.D0-2.D0*ag*(molevg+bg)/(molevg*(molevg+bg)+bg*(molevg-bg))**2.D0)
                    if(local_z(m) < EQUALPREC) then
                        local_v(m) = 0.D0
                    end if
                end do
                local_Cf = CfG
                x = 0.D0
                y = local_z
                xiL = 0.D0
                rhoL = 0.D0
                sL = 0.D0
            end if

            deallocate(tempsum)

        else

            do m = 1, Nc
                if(local_z(m) < EQUALPREC) then
                    x(m) = 0.D0
                    y(m) = 0.D0
                end if
            end do

            sumx = 0.D0
            sumy = 0.D0
            do m = 1, Nc
                sumx = sumx + x(m)
                sumy = sumy + y(m)
            end do
            do m = 1, Nc
                x(m) = x(m) / sumx
                y(m) = y(m) / sumy
            end do

            call PREOS(x, y, P, ZL, ZG, am, bm, al, ag, bl, bg, bigAL, &!
                bigAG, bigBL, bigBG, xiL, xiG, rhoL, rhoG, CfL, CfG)

            do m = 1, Nc
                if(local_z(m) > EQUALPREC) then
                    nn1 = m
                    exit
                end if
            end do
            do m = nn1+1, Nc
                if(local_z(m) > EQUALPREC) then
                    nn2 = m
                    exit
                end if
            end do

            sL = (y(nn2)*local_z(nn1)-y(nn1)*local_z(nn2))/(y(nn2)*local_z(nn1)-y(nn1)*local_z(nn2)+ &!
                (x(nn1)*local_z(nn2)-x(nn2)*local_z(nn1))*xiL/xiG)

            phasemole(1) = (y(nn2)*local_z(nn1)-y(nn1)*local_z(nn2))/(x(nn1)*y(nn2)-x(nn2)*y(nn1))* &!
                (sL*xiL+(1.D0-sL)*xiG)
            phasemole(2) = (x(nn1)*local_z(nn2)-x(nn2)*local_z(nn1))/(x(nn1)*y(nn2)-x(nn2)*y(nn1))* &!
                (sL*xiL+(1.D0-sL)*xiG)

            allocate(bigK(Nc))
            do m = 1, Nc
                bigK(m) = dexp(5.37D0*(1+af(m))*(1.D0-ct(m)/Temp)+dlog(cp(m)/P*1.D0))
            end do
            do m = 1, Nc
                if(local_z(m) < EQUALPREC) then
                    local_v(m) = 0.D0
                else
                    local_v(m) = pmv(m, P, local_z, 1.D0/(phasemole(1)+phasemole(2)), bigK)
                end if
            end do
            deallocate(bigK)

            local_Cf = comprefac(x, y, P, ZL, ZG, am, bm, al, ag, bl, bg, &!
                bigAL, bigAG, bigBL, bigBG, xiL, xiG, sL)

        end if

        deallocate(am)
        deallocate(bm)
        
    end subroutine flashcalculation_sparsegrid

end module RST_flashcalculation_sparsegrid
