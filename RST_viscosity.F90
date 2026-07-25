
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_viscosity

    use RST_globalFlashData

    implicit none

contains

    function viscosity( x, xi, P, phase ) result(mu)

        real(kind=8), dimension(:), pointer, intent(in) :: x
        real(kind=8), intent(in) :: xi
        real(kind=8), intent(in) :: P
        character, intent(in) :: phase
        real(kind=8) :: mu

        real(kind=8), dimension(:), pointer :: local_cp, local_mw
        real(kind=8), dimension(:), pointer :: CI, XMIU
        real(kind=8) :: RO, TR, S, S1, VS, SV, Z7PLUS, MW7PLUS, GM7PLUS, V7PLUS, ROC, ROR, ST, SM, SPI, F, COR
        logical :: isHeavy
        real(kind=8) :: local_A(16), local_B(4)
        real(kind=8) :: PPC, TPC, GML, VS1, PPR, TPR, CE, VISR
        integer :: i, j, m

        allocate(local_cp(Nc))
        allocate(local_mw(Nc))

        local_cp(1:Nc) = 9.8692D-6*cp(1:Nc)
        local_mw(1:Nc) = mw(1:Nc)*1.D3

        if(phase == 'l') then

            allocate(CI(Nc))
            do i = 1, Nc
                CI(i) = 1.D0/(ct(i)**(1.D0/6.D0)/sqrt(local_mw(i))/local_cp(i)**6.66666666666D-1)
            end do

            RO = xi*1.D-3

            allocate(XMIU(Nc))
            do i = 1, Nc
                TR = Temp/ct(i)
                if (TR <= 1.5D0) then
                    XMIU(i)=3.4D-4*CI(i)*TR**9.4D-1
                else
                    XMIU(i)=1.778D-4*CI(i)*(4.58D0*TR-1.67D0)**6.25D-1
                end if
            end do

            S = 0.D0
            S1 = 0.D0
            GML = 0.D0

            do i = 1, Nc
                S = S + x(i)*XMIU(i)*sqrt(local_mw(i))
                S1 = S1 + x(i)*sqrt(local_mw(i))
                GML = GML + x(i)*local_mw(i)
            end do
            RO = RO*GML
            VS = S/S1
            SV = 0.D0

            isHeavy = .false.
            if(.not.isHeavy) then
                do i = 1, Nc
                    SV = SV + x(i)*cv(i)*local_mw(i)
                end do
            else
                do i = 1, 6
                    SV = SV + x(i)*cv(i)
                end do

                Z7PLUS = 0.D0
                do i = 7, Nc
                    Z7PLUS = Z7PLUS + x(i)
                end do
                MW7PLUS = 4.307D2
                GM7PLUS = 9.88D-1 ! AS GIVEN IN THE HEAVY OIL CASE
                V7PLUS = 2.1573D1 + 1.5122D-2*MW7PLUS - 2.7656D1*GM7PLUS + 7.0615D-2*MW7PLUS*GM7PLUS
                V7PLUS = V7PLUS*3.048D-1**3/4.53592D-1 ! cub-ft/lb-mol -> cub-m/kg-mol

                SV = SV + Z7PLUS*V7PLUS
            end if

            ROC = 1.D0/SV
            if(.not.isHeavy) THEN
                ROR = RO/ROC/GML
            else
                ROR = RO/ROC
            end if

            ST = 0.D0
            SM = 0.D0
            SPI = 0.D0

            do i = 1, Nc
                SM = SM + x(i)*local_mw(i)
                ST = ST + x(i)*ct(i)
                SPI = SPI + x(i)*local_cp(i)
            end do

            CE = ST**(1.D0/6.D0)/sqrt(SM)/SPI**(2.D0/3.D0)

            ! AUGUST 16, 2007: If ROR is greater than 10, F is suspiciously large ...
            F = 1.023D-1 + 2.3364D-2*ROR + 5.8533D-2*ROR**2.D0 - 4.0758D-2*ROR**3.D0  + 9.3324D-3*ROR**4.D0

            COR = (F**4.D0-1.D-4)/CE
            mu = VS + COR

            deallocate(CI)
            deallocate(XMIU)

        else if(phase == 'g') then

            local_A(1) = -2.4621182D0
            local_A(2) = 2.97054714D0
            local_A(3) = -2.86264054D-1
            local_A(4) = 8.05420522*1.D-3
            local_A(5) = 2.80860949D0
            local_A(6) = -3.49803305D0
            local_A(7) = 3.6037302D-1
            local_A(8) = -1.04432413D-2
            local_A(9) = -7.93385684D-1
            local_A(10) = 1.39643306D0
            local_A(11) = -1.49144925D-1
            local_A(12) = 4.41015512*1.D-3
            local_A(13) = 8.39387178D-2
            local_A(14) = -1.86408848D-1
            local_A(15) = 2.03367881D-2
            local_A(16) = -6.09579263*1.D-4

            local_B(1) = 4.D0
            local_B(2) = 4.D0
            local_B(3) = 4.D0
            local_B(4) = 4.D0

            PPC = 0.D0
            TPC = 0.D0
            GML = 0.D0

            do m = 1, Nc
                PPC = PPC + x(m)*local_cp(m)
                TPC = TPC + x(m)*ct(m)
                GML = GML + x(m)*local_mw(m)
            end do

            VS1 = (7.43D0+1.33D-2*GML)*(1.8D0*Temp)**1.5D0 &!
                /(1.8D0*Temp+7.54D1+1.39D1*GML)*1.D0*1.D-4
            PPR = P*1.D-5/PPC
            TPR = Temp/TPC
            do i = 1, 4
                j = 4*i
                local_B(i) = local_A(j-3)+(local_A(j-2)+(local_A(j-1)+local_A(j)*PPR)*PPR)*PPR
            end do

            CE = local_B(1)+(local_B(2)+(local_B(3)+local_B(4)*TPR)*TPR)*TPR
            VISR = dexp(CE*1.D0)/TPR
            mu = VISR*VS1

        end if

        mu = mu * 1.D-3

        deallocate(local_cp)
        deallocate(local_mw)

    end function viscosity

end module RST_viscosity
