
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_stability

    use RST_mathlib
    use RST_globalFlashData
    use RST_fugacitycoef

    implicit none

contains

    ! BFGS-Quasi-Newton-Method
    subroutine stability( Pk, ZI, ztest, Kstabtry, size_try, K_stab, judge )

        real(kind=8), intent(in) :: Pk
        real(kind=8), dimension(:), pointer, intent(in) :: ZI
        real(kind=8), dimension(:), pointer, intent(in) :: ztest
        real(kind=8), dimension(:,:), pointer, intent(in) :: Kstabtry
        integer, intent(in) :: size_try
        real(kind=8), dimension(:,:), pointer, intent(in out) :: K_stab
        logical, intent(in out) :: judge

        real(kind=8), dimension(:), pointer :: TPD_stab
        integer :: size_stab
        integer, dimension(:), pointer :: iter_stab
        real(kind=8), dimension(:), pointer :: XE, XOLD, XNEW, GOLD, FL, FG, FUGZ, S, YD, YE, GNEW, X_CAL
        real(kind=8) :: SUM0, TM
        integer :: i_stab, i_comp, j_comp, iteration_stab, m_stab
        integer :: num_judge_iter, num_judge_comp
        real(kind=8) :: err_stab
        real(kind=8), dimension(:,:), pointer :: matrix_stab
        real(kind=8), dimension(:), pointer :: phil, phig, x, y
        real(kind=8) :: xiL, xiG, rhoL, rhoG, CfL, CfG
        real(kind=8) :: ZL, ZG
        real(kind=8), dimension(:), pointer :: am, bm
        real(kind=8) :: al, ag, bl, bg, bigAL, bigAG, bigBL, bigBG

        allocate(XE(Nc))
        allocate(YE(Nc))
        allocate(XOLD(Nc))
        allocate(XNEW(Nc))
        allocate(GOLD(Nc))
        allocate(GNEW(Nc))
        allocate(FG(Nc))
        allocate(FL(Nc))
        allocate(FUGZ(Nc))
        allocate(S(Nc))
        allocate(YD(Nc))
        allocate(X_CAL(Nc))
        allocate(am(Nc))
        allocate(bm(Nc))
        allocate(phil(Nc))
        allocate(phig(Nc))
        allocate(x(Nc))
        allocate(y(Nc))
        allocate(iter_stab(Nc+8))
        allocate(TPD_stab(Nc+8))

        judge = .false.

        size_stab = 0
        num_judge_iter = 0
        TPD_stab = 0.D0
        K_stab = 0.D0

        y = 0.D0

        call fugacitycoef( y, ZI, Pk, ZL, ZG, am, bm, al, ag, bl, bg, bigAL, bigAG, &!
            bigBL, bigBG, xiL, xiG, rhoL, rhoG, CfL, CfG, phil, phig )

        FG = DLOG(phig*1.D0)

        do i_comp = 1, Nc
            if(ztest(i_comp) <= 0.D0) then
                cycle
            end if
            FUGZ(i_comp) = FG(i_comp)+DLOG(ztest(i_comp)*1.D0)
        end do

        do i_stab = 1, size_try

            do i_comp = 1, Nc
                XE(i_comp) = ztest(i_comp) / Kstabtry(i_stab,i_comp)
            end do

            X_CAL = XE/SUM(XE)

            x = 0.D0

            call fugacitycoef( X_CAL, x, Pk, ZL, ZG, am, bm, al, ag, bl, bg, bigAL, &!
                bigAG, bigBL, bigBG, xiL, xiG, rhoL, rhoG, CfL, CfG, phil, phig )

            FL = DLOG(phil*1.D0)

            do i_comp = 1, Nc
                if((ztest(i_comp) <= 0.D0).or.(xe(i_comp) <= 0.D0)) then
                    XOLD(i_comp) = 0.D0
                    GOLD(i_comp) = 0.D0
                    XNEW(i_comp) = 0.D0
                    cycle
                end if
                XOLD(i_comp) = 2.D0*DSQRT(XE(i_comp)*1.D0)
                TM = FUGZ(i_comp)-FL(i_comp)
                GOLD(i_comp) = XOLD(i_comp)/2.D0*(DLOG(XE(i_comp)*1.D0)-TM)
                XNEW(i_comp) = 2.D0*DEXP(TM/2.D0)
            end do

            iteration_stab = 0

            main_loopg: &
            do while(iteration_stab <= 100)

                XE(1:Nc) = 2.5D-1*XNEW(1:Nc)**2.D0
                X_CAL = XE/SUM(XE)

                x = 0.D0
                call fugacitycoef( X_CAL, x, Pk, ZL, ZG, am, bm, al, ag, bl, bg, bigAL, &!
                    bigAG, bigBL, bigBG, xiL, xiG, rhoL, rhoG, CfL, CfG, phil, phig )

                FL = DLOG(phil*1.D0)
                SUM0=1.D0
                do i_comp = 1, Nc
                    if((ztest(i_comp) <= 0.D0).or.(xe(i_comp) <= 0.D0)) then
                        GNEW(i_comp) = 0.D0
                        S(i_comp) = 0.D0
                        YD(i_comp) = 0.D0
                        cycle
                    end if

                    TM = DLOG(XE(i_comp)*1.D0)+FL(i_comp)-FUGZ(i_comp)
                    SUM0 = SUM0+XE(i_comp)*(TM-1.D0)
                    GNEW(i_comp) = XNEW(i_comp)/2.D0*TM
                    S(i_comp) = XNEW(i_comp)-XOLD(i_comp)
                    YD(i_comp) = GNEW(i_comp)-GOLD(i_comp)
                end do

                call UPDATE(Nc,YD,S,GNEW,YE)

                do i_comp = 1, Nc
                    GOLD(i_comp) = GNEW(i_comp)
                    XOLD(i_comp) = XNEW(i_comp)
                    XNEW(i_comp) = XNEW(i_comp)-YE(i_comp)
                end do

                err_stab = 0.D0
                do i_comp = 1, Nc
                    if(abs(xnew(i_comp)-xold(i_comp))>=err_stab) then
                        err_stab=abs(xnew(i_comp)-xold(i_comp))
                    end if
                end do

                if(err_stab < 1.D-10) then
                    exit main_loopg
                end if

                if(err_stab>1.D30) then
                    exit main_loopg
                end if

                iteration_stab = iteration_stab + 1

            end do main_loopg

            if(err_stab >= 1.D-10) then
                num_judge_iter = num_judge_iter + 1
            end if

            if(err_stab < 1.D-10) then                   !Convergence criteria
                num_judge_comp = 0
                do i_comp = 1, Nc
                    if (abs(FL(i_comp)-FG(i_comp)) <= 1.D-5) then
                        num_judge_comp = num_judge_comp + 1
                    end if
                end do
                if(num_judge_comp < Nc) then
                    size_stab = size_stab + 1

                    iter_stab(size_stab) = iteration_stab
                    TPD_stab(size_stab) = sum0
                    do i_comp = 1, Nc
                        K_stab(size_stab,i_comp) = dexp((FL(i_comp)-FG(i_comp))*1.D0)
                    enddo
                endif
            endif

        enddo

        if(size_stab == 0) then

            judge = .true.
            if(num_judge_iter == size_try) then
                print *, 'NO CONVERGENCE IN STABILITY: ASSUMED SINGLE PHASE'
            endif

        else if(size_stab >= 2) then

            allocate(matrix_stab(size_stab,Nc+2))
            matrix_stab(:,1:Nc) = K_stab(1:size_stab,1:Nc)
            matrix_stab(:,Nc+1) = dfloat(iter_stab(1:size_stab))
            matrix_stab(:,Nc+2) = TPD_stab(1:size_stab)
            call piksrt(size_stab,Nc+2,matrix_stab)
            K_stab(1:size_stab,1:Nc) = matrix_stab(:,1:Nc)
            iter_stab(1:size_stab) = dnint(matrix_stab(:,Nc+1)*1.D0)
            TPD_stab(1:size_stab) = matrix_stab(:,Nc+2)

            m_stab = 1
            do while(m_stab < size_stab)
                num_judge_comp = 0
                do i_comp = 1, Nc
                    if (abs(K_stab(m_stab,i_comp)/K_stab(m_stab+1,i_comp)-1.D0)<=1.D-5) then
                        num_judge_comp = num_judge_comp + 1
                    end if
                end do
                if (num_judge_comp == Nc) then
                    size_stab = size_stab - 1
                    TPD_stab(m_stab+1:size_stab) = TPD_stab(m_stab+2:size_stab+1)
                    iter_stab(m_stab+1:size_stab) = iter_stab(m_stab+2:size_stab+1)
                    K_stab(m_stab+1:size_stab,:) = K_stab(m_stab+2:size_stab+1,:)
                else
                    m_stab = m_stab + 1
                end if
            end do

            TPD_stab(size_stab+1:Nc+8) = 0.D0
            K_stab(size_stab+1:Nc+8,:) = 0.D0
            deallocate(matrix_stab)

        endif

        if (TPD_stab(1) >= -1.D-10) then
            judge = .true.
        endif

        deallocate(XE)
        deallocate(YE)
        deallocate(XOLD)
        deallocate(XNEW)
        deallocate(GOLD)
        deallocate(GNEW)
        deallocate(FG)
        deallocate(FL)
        deallocate(FUGZ)
        deallocate(S)
        deallocate(YD)
        deallocate(X_CAL)
        deallocate(am)
        deallocate(bm)
        deallocate(phil)
        deallocate(phig)
        deallocate(x)
        deallocate(y)
        deallocate(iter_stab)
        deallocate(TPD_stab)

    end subroutine stability

    subroutine update(N, Y1, SS, GD, X1)

        integer :: i, N
        real(kind=8), dimension(:), pointer :: Y1, SS, GD, X1
        real(kind=8) :: SY, HYY, XI, FAC, XG, SG0

        SY = 0.D0
        HYY = 0.D0
        do i = 1, N
            X1(I) = Y1(I)
            XI = X1(I)
            SY = SY+SS(I)*Y1(I)
            HYY = HYY+XI*Y1(I)
        end do

        FAC = (1.D0+HYY/SY)*5.D-1
        X1 = (X1-FAC*SS)/SY

        XG = 0.D0
        SG0 = 0.D0
        do i = 1, N
            XG = XG+X1(I)*GD(I)
            SG0 = SG0+SS(I)*GD(I)
        end do

        X1 = GD-XG*SS-SG0*X1

    end subroutine update

end module RST_stability
