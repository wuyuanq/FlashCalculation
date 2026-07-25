
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_flashcalculation

    use RST_globalFlashData
    use RST_PREOS
    use RST_pmv
    use RST_stability
    use RST_comprefac

    implicit none

    real(kind=8), parameter :: FEQUALPREC = 1.D-12

contains

    subroutine flashcalculation(P, local_z, x, y, xiL, xiG, rhoL, &!
        rhoG, sL, local_v, local_Cf, isW, isN, isRea)

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
        logical, intent(out) :: isW, isN, isRea

        integer :: purenum
        real(kind=8), dimension(:), pointer :: K
        real(kind=8) :: h0, h1, yt, xt, beta, lp, rp, hbeta
        real(kind=8) :: CfL, CfG, criteria, VL, VG, mt, mv, Vold
        logical :: isPure
        real(kind=8) :: ZL, ZG
        real(kind=8), dimension(:), pointer :: am, bm
        real(kind=8) :: al, ag, bl, bg, bigAL, bigAG, bigBL, bigBG
        real(kind=8) :: phasemole(2)
        integer :: size
        real(kind=8), dimension(:), pointer :: psat
        real(kind=8) :: psatmix
        real(kind=8) :: molevg, molevl
        real(kind=8), dimension(:), pointer :: tempsum
        integer :: i, j, local_t, m, nn1, nn2
        real(kind=8), dimension(:), pointer  :: ztest1f
        real(kind=8), dimension(:,:), pointer :: Kstab1ftry, K_stab_1F
        integer :: size1f_try
        logical :: SINGLE
        real(kind=8) :: sum

        ! decide whether it is pure substance
        isPure = .false.
        purenum = 0
        do m = 1, Nc
            if(abs(local_z(m)-1.D0) < FEQUALPREC) then
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

        allocate(K(Nc))
        do i = 1, Nc
            K(i) = dexp(5.37D0*(1+af(i))*(1.D0-ct(i)/Temp)+dlog(cp(i)/P*1.D0))
        end do

        allocate(ztest1f(Nc))
        ztest1f = local_z

        allocate(Kstab1ftry(Nc+8, Nc))
        allocate(K_stab_1F(Nc+8, Nc))

        Kstab1ftry(1,:) = K
        Kstab1ftry(2,:) = 1.D0/K
        Kstab1ftry(3,:) = K**(1.D0/3.D0)
        Kstab1ftry(4,:) = 1.D0/K**(1.D0/3.D0)
        do i = 1, Nc
            do j = 1, Nc
                if (i == j) then
                    Kstab1ftry(i+4,j) = 9.D-1/ztest1f(j)
                else
                    Kstab1ftry(i+4,j) = 1.D-1/dfloat(Nc-1)/ztest1f(j)
                endif
            end do
        end do
        Kstab1ftry(Nc+5:Nc+8,:) = 0.D0
        size1f_try = Nc+4

        call stability( P, local_z, ztest1f, Kstab1ftry, size1f_try, K_stab_1F, SINGLE )

        deallocate(Kstab1ftry)
        deallocate(ztest1f)

        allocate(am(Nc))
        allocate(bm(Nc))

        if(SINGLE) then

            call PREOS( local_z, local_z, P, ZL, ZG, am, bm, al, ag, bl, bg, bigAL, &!
                bigAG, bigBL, bigBG, xiL, xiG, rhoL, rhoG, CfL, CfG )

            allocate(tempsum(Nc))

            do i = 1, Nc
                tempsum(i) = 0.D0
                do j = 1, Nc
                    tempsum(i) = tempsum(i) + local_z(j)*(1.D0-delta(i,j))*dsqrt(am(i)*am(j)*1.D0)
                end do
                tempsum(i) = 2.D0*tempsum(i)
            end do

            mt = 0.D0
            do m = 1, Nc
                mt = mt + ct(m)*local_z(m)
            end do

            if(Temp > mt) then
                x = 0.D0
                y = local_z
                xiL = 0.D0
                rhoL = 0.D0
                sL = 0.D0
                local_Cf = CfG

                molevg = ZG*R*Temp/P
                do m = 1, Nc
                    local_v(m) = (R*Temp/(molevg-bg)*(1.D0+bm(m)/(molevg-bg))-(tempsum(m)-(2.D0*ag*bm(m)*(molevg-bg)) &!
                        /(molevg*(molevg+bg)+bg*(molevg-bg)))/(molevg*(molevg+bg)+bg*(molevg-bg)))/(R*Temp &!
                        /(molevg-bg)**2.D0-2.D0*ag*(molevg+bg)/(molevg*(molevg+bg)+bg*(molevg-bg))**2.D0)
                    if(abs(local_z(m)) < FEQUALPREC) then
                        local_v(m) = 0.D0
                    end if
                end do

                isW = .false.
                isN = .true.
            else
                !compute the saturation vapor pressure of each component at temperature Temp, Antoine equation
                allocate(psat(Nc))
                do m = 1, Nc
                    psat(m) = dexp((psatA(m)-psatB(m)/((Temp-2.73D2)+psatC(m)))*1.D0)
                    psat(m) = psat(m)*1.33322D2
                end do
                psatmix = 0.D0
                do m = 1, Nc
                    psatmix = psatmix + psat(m)*local_z(m)
                end do
                if(P > psatmix) then
                    x = local_z
                    y = 0.D0
                    xiG = 0.D0
                    rhoG = 0.D0
                    sL = 1.D0
                    local_Cf = CfL

                    molevl = ZL*R*Temp/P
                    do m = 1, Nc
                        local_v(m) = (R*Temp/(molevl-bl)*(1.D0+bm(m)/(molevl-bl))-(tempsum(m)-(2.D0*al*bm(m)* &!
                            (molevl-bl))/(molevl*(molevl+bl)+bl*(molevl-bl)))/(molevl*(molevl+bl)+bl* &!
                            (molevl-bl)))/(R*Temp/(molevl-bl)**2.D0-2.D0*al*(molevl+bl)/(molevl*(molevl+bl)+ &!
                            bl*(molevl-bl))**2.D0)
                        if(abs(local_z(m)) < FEQUALPREC) then
                            local_v(m) = 0.D0
                        end if
                    end do

                    isN = .false.
                    isW = .true.
                else
                    x = 0.D0
                    y = local_z
                    xiL = 0.D0
                    rhoL = 0.D0
                    sL = 0.D0
                    local_Cf = CfG

                    molevg = ZG*R*Temp/P
                    do m = 1, Nc
                        local_v(m) = (R*Temp/(molevg-bg)*(1.D0+bm(m)/(molevg-bg))-(tempsum(m)-(2.D0*ag*bm(m)* &!
                            (molevg-bg))/(molevg*(molevg+bg)+bg*(molevg-bg)))/(molevg*(molevg+bg)+bg* &!
                            (molevg-bg)))/(R*Temp/(molevg-bg)**2.D0-2.D0*ag*(molevg+bg)/(molevg*(molevg+bg)+ &!
                            bg*(molevg-bg))**2.D0)
                        if(abs(local_z(m)) < FEQUALPREC) then
                            local_v(m) = 0.D0
                        end if
                    end do

                    isW = .false.
                    isN = .true.
                end if
                deallocate(psat)
            end if

            deallocate(tempsum)

        else

            K = K_stab_1F(1,:)

20          do local_t = 1, MAXTIME

                h0 = 0.D0
                do i = 1, Nc
                    h0 = h0 + K(i)*local_z(i)
                end do
                h0 = h0 - 1.D0
                h1 = 0.D0
                do i = 1, Nc
                    h1 = h1 + local_z(i)/K(i)
                end do
                h1 = 1.D0 - h1

                if(h0 <= 0.D0) then
                    x = local_z
                    yt = 0.D0
                    do i = 1, Nc
                        yt = yt + K(i)*local_z(i)
                    end do
                    y = K*local_z/yt
                    go to 10
                end if
                
                if(h1 >= 0.D0) then
                    y = local_z
                    xt = 0.D0
                    do i = 1, Nc
                        xt = xt + local_z(i)/K(i)
                    end do
                    x = local_z/K/xt
                    go to 10
                end if

                beta = 5.D-1
                lp = 0.D0
                rp = 1.D0

                do while(.true.)
                    hbeta = 0.D0
                    do i = 1, Nc
                        hbeta = hbeta + ((K(i)-1.D0)*local_z(i))/(1.D0+beta*(K(i)-1.D0))
                    end do

                    if(abs(hbeta) < FEQUALPREC) then
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
                    x(i) = local_z(i)/(1.D0+(K(i)-1.D0)*beta)
                    y(i) = K(i)*x(i)
                end do

10              call substeps( x, y, P, criteria, K, ZL, ZG, am, bm, al, ag, bl, bg, &!
                    bigAL, bigAG, bigBL, bigBG, xiL, xiG, rhoL, rhoG, CfL, CfG )

                if(rhoL <= rhoG) then
                    call substeps( y, x, P, criteria, K, ZL, ZG, am, bm, al, ag, bl, bg, &!
                        bigAL, bigAG, bigBL, bigBG, xiL, xiG, rhoL, rhoG, CfL, CfG )
                end if

                if(criteria < FEQUALPREC) then

                    do i = 1, Nc
                        if(local_z(i) > FEQUALPREC) then
                            nn1 = i
                            exit
                        end if
                    end do
                    do i = nn1+1, Nc
                        if(local_z(i) > FEQUALPREC) then
                            nn2 = i
                            exit
                        end if
                    end do

                    sL = (y(nn2)*local_z(nn1)-y(nn1)*local_z(nn2))/(y(nn2)*local_z(nn1)-y(nn1)*local_z(nn2)+ &!
                        (x(nn1)*local_z(nn2)-x(nn2)*local_z(nn1))*xiL/xiG)

                    phasemole(1) = (y(nn2)*local_z(nn1)-y(nn1)*local_z(nn2))/(x(nn1)*y(nn2)-x(nn2)*y(nn1))* &!
                        (sL*xiL+(1.D0-sL)*xiG)
                    phasemole(2) = (x(nn1)*local_z(nn2)-x(nn2)*local_z(nn1))/(x(nn1)*y(nn2)-x(nn2)*y(nn1))* &!
                        (sL*xiL+(1.D0-sL)*xiG)

                    do i = 1, Nc
                        local_v(i) = pmv( i, P, local_z, 1.D0/(phasemole(1)+phasemole(2)), K )
                    end do

                    local_Cf = comprefac(x, y, P, ZL, ZG, am, bm, al, ag, bl, bg, &!
                        bigAL, bigAG, bigBL, bigBG, xiL, xiG, sL)

                    ! When sL is very small, the algorithm to compute local_Cf in 2 phases will return NaN,
                    ! so we have to process such condition.
                    if(local_Cf.ne.local_Cf) then
                        if(sL < FEQUALPREC) then
                            local_Cf = CfG
                        else
                            local_Cf = CfL
                        end if
                    end if

                    isW = .true.
                    isN = .true.

                    exit
                end if

            end do
        end if

        deallocate(K)
        deallocate(K_stab_1F)
        deallocate(am)
        deallocate(bm)

        ! Decide whether the results are reasonable.
        isRea = .true.
        do i = 1, Nc
            if((x(i)<0.D0).or.(x(i)>1.D0)) then
                isRea = .false.
                exit
            end if
            if((y(i)<0.D0).or.(y(i)>1.D0)) then
                isRea = .false.
                exit
            end if
        end do
        if(isRea) then
            if(xiL < 0.D0) then
                isRea = .false.
            end if
        end if
        if(isRea) then
            if(xiG < 0.D0) then
                isRea = .false.
            end if
        end if
        if(isRea) then
            if(isW.and.isN.and.(xiL<xiG)) then
                isRea = .false.
            end if
        end if
        if(isRea) then
            if(rhoL < 0.D0) then
                isRea = .false.
            end if
        end if
        if(isRea) then
            if(rhoG < 0.D0) then
                isRea = .false.
            end if
        end if
        if(isRea) then
            if(isW.and.isN.and.(rhoL<rhoG)) then
                isRea = .false.
            end if
        end if
        if(isRea) then
            if((sL<0.D0).or.(sL>1.D0)) then
                isRea = .false.
            end if
        end if
        if(isRea) then
            do i = 1, Nc
                if(abs(local_v(i))>1.D1) then
                    isRea = .false.
                    exit
                end if
            end do
        end if
        if(isRea) then
            if(abs(local_Cf)>1.D1) then
                isRea = .false.
            end if
        end if

    end subroutine flashcalculation

end module RST_flashcalculation
