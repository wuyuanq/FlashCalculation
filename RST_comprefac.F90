
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_comprefac

    use RST_mathlib
    use RST_globalFlashData

    implicit none

contains

    function comprefac(x, y, P, ZL, ZG, am, bm, al, ag, bl, bg, &!
        bigAL, bigAG, bigBL, bigBG, xiL, xiG, sL) result(local_Cf)

        real(kind=8), dimension(:), pointer, intent(in) :: x, y
        real(kind=8), intent(in) :: P
        real(kind=8), intent(in) :: ZL, ZG
        real(kind=8), dimension(:), pointer, intent(in) :: am, bm
        real(kind=8), intent(in) :: al, ag, bl, bg, bigAL, bigAG, bigBL, bigBG
        real(kind=8), intent(in) :: xiL, xiG
        real(kind=8), intent(in) :: sL
        real(kind=8) :: local_Cf

        real(kind=8), dimension(:), pointer :: tempsumg, tempsuml
        real(kind=8) :: parZGbigAG, parZLbigAL, parZGbigBG, parZLbigBL
        real(kind=8) :: ng, nl
        real(kind=8), dimension(:), pointer :: paragy, paralx
        real(kind=8), dimension(:), pointer :: parbigAGy, parbigALx, parbigBGy, &!
            parbigBLx, parZGy, parZLx, par2gy, par2lx, par31gy, par31lx, par33gy, &!
            par33lx, rightgy, rightlx
        real(kind=8), dimension(:,:), pointer :: par1gy, par1lx, par32gy, par32lx, &!
            par3gy, par3lx, parfgy, parflx, parfgn, parfln
        real(kind=8) :: parZGp, parZLp, parbigAGp, parbigALp, parbigBGp, parbigBLp, &!
            par31gp, par31lp, par33gp, par33lp
        real(kind=8), dimension(:), pointer :: par1gp, par1lp, par2gp, par2lp, &!
            par3gp, par3lp, parfgp, parflp
        real(kind=8), dimension(:), pointer :: parnlp, parngp, parbigAGn, parbigALn, &!
            parbigBGn, parbigBLn, parZGn, parZLn
        real(kind=8) :: parZGPtotaln, parZLPtotaln, parngptotaln, parnlptotaln, parvfp
        real(kind=8), dimension(:,:), pointer :: local_A
        real(kind=8), dimension(:), pointer :: local_b
        real(kind=8), dimension(:), pointer :: solx
        integer, dimension(:), pointer :: JS
        integer :: L
        integer :: i, j

        allocate(solx(Nc))
        allocate(JS(Nc))
        allocate(local_A(Nc,Nc))
        allocate(local_b(Nc))
        local_A(:,:) = 0.D0
        local_b(:) = 0.D0

        allocate(tempsumg(Nc))
        allocate(tempsuml(Nc))
        allocate(paragy(Nc))
        allocate(paralx(Nc))
        allocate(parbigAGy(Nc))
        allocate(parbigALx(Nc))
        allocate(parbigBGy(Nc))
        allocate(parbigBLx(Nc))
        allocate(parZGy(Nc))
        allocate(parZLx(Nc))
        allocate(par2gy(Nc))
        allocate(par2lx(Nc))
        allocate(par31gy(Nc))
        allocate(par31lx(Nc))
        allocate(par33gy(Nc))
        allocate(par33lx(Nc))
        allocate(rightgy(Nc))
        allocate(rightlx(Nc))
        allocate(par1gy(Nc,Nc))
        allocate(par1lx(Nc,Nc))
        allocate(par32gy(Nc,Nc))
        allocate(par32lx(Nc,Nc))
        allocate(par3gy(Nc,Nc))
        allocate(par3lx(Nc,Nc))
        allocate(parfgy(Nc,Nc))
        allocate(parflx(Nc,Nc))
        allocate(parfgn(Nc,Nc))
        allocate(parfln(Nc,Nc))
        allocate(par1gp(Nc))
        allocate(par1lp(Nc))
        allocate(par2gp(Nc))
        allocate(par2lp(Nc))
        allocate(par3gp(Nc))
        allocate(par3lp(Nc))
        allocate(parfgp(Nc))
        allocate(parflp(Nc))
        allocate(parnlp(Nc))
        allocate(parngp(Nc))
        allocate(parbigAGn(Nc))
        allocate(parbigALn(Nc))
        allocate(parbigBGn(Nc))
        allocate(parbigBLn(Nc))
        allocate(parZGn(Nc))
        allocate(parZLn(Nc))

        do i = 1, Nc
            tempsumg(i) = 0.D0
            do j = 1, Nc
                tempsumg(i) = tempsumg(i) + y(j)*(1-delta(i,j))*dsqrt(am(i)*am(j)*1.D0)
            end do
            tempsumg(i) = 2.D0*tempsumg(i)
        end do

        do i = 1, Nc
            tempsuml(i) = 0.D0
            do j = 1, Nc
                tempsuml(i) = tempsuml(i) + x(j)*(1-delta(i,j))*dsqrt(am(i)*am(j)*1.D0)
            end do
            tempsuml(i) = 2.D0*tempsuml(i)
        end do

        parZGbigAG = (bigBG-ZG)/(3.D0*ZG**2.D0-2.D0*ZG*(1.D0-bigBG)+(bigAG-3.D0*bigBG**2.D0-2.D0*bigBG))
        parZLbigAL = (bigBL-ZL)/(3.D0*ZL**2.D0-2.D0*ZL*(1.D0-bigBL)+(bigAL-3.D0*bigBL**2.D0-2.D0*bigBL))

        parZGbigBG = (-ZG**2.D0+2.D0*(3.D0*bigBG+1.D0)*ZG+(bigAG-2.D0*bigBG-3.D0*bigBG**2.D0))/ &!
            (3.D0*ZG**2.D0-2.D0*(1.D0-bigBG)*ZG+(bigAG-3.D0*bigBG**2.D0-2.D0*bigBG))
        parZLbigBL = (-ZL**2.D0+2.D0*(3.D0*bigBL+1.D0)*ZL+(bigAL-2.D0*bigBL-3.D0*bigBL**2.D0))/ &!
            (3.D0*ZL**2.D0-2.D0*(1.D0-bigBL)*ZL+(bigAL-3.D0*bigBL**2.D0-2.D0*bigBL))

        ! suppose Vf = 1, then ng = xiG*(1-sL)*Vf = xiG*(1-sL)
        ng = xiG*(1.D0-sL)
        nl = xiL*sL

        do i = 1, Nc

            paragy(i) = 0.D0
            do j = 1, Nc
                if(j /= i) then
                    paragy(i) = paragy(i) + y(j)*(1.D0-delta(i,j))*dsqrt(am(i)*am(j)*1.D0)
                end if
            end do
            paragy(i) = 2.D0*paragy(i)

            paralx(i) = 0.D0
            do j = 1, Nc
                if(j /= i) then
                    paralx(i) = paralx(i) + x(j)*(1.D0-delta(i,j))*dsqrt(am(i)*am(j)*1.D0)
                end if
            end do
            paralx(i) = 2.D0*paralx(i)

            parbigAGy(i) = P/(R*Temp)**2.D0*paragy(i)
            parbigALx(i) = P/(R*Temp)**2.D0*paralx(i)

            parbigBGy(i) = P*bm(i)/R/Temp
            parbigBLx(i) = P*bm(i)/R/Temp

            parZGy(i) = parZGbigAG*P/(R*Temp)**2.D0*paragy(i) + parZGbigBG*P*bm(i)/R/Temp

            parZLx(i) = parZLbigAL*P/(R*Temp)**2.D0*paralx(i) + parZLbigBL*P*bm(i)/R/Temp

            par2gy(i) = 1/(ZG-bigBG)*(parZGy(i)-parbigBGy(i))
            par2lx(i) = 1/(ZL-bigBL)*(parZLx(i)-parbigBLx(i))

            par31gy(i) = parbigAGy(i)/(2.D0*sqrt(2.D0)*bigBG)-bigAG/(2.D0*sqrt(2.D0)*bigBG**2.D0) &!
                *parbigBGy(i)
            par31lx(i) = parbigALx(i)/(2.D0*sqrt(2.D0)*bigBL)-bigAL/(2.D0*sqrt(2.D0)*bigBL**2.D0) &!
                *parbigBLx(i)

            par33gy(i) = ((parZGy(i)+2.414D0*parbigBGy(i))*(ZG-4.14D-1*bigBG)-(ZG+2.414D0*bigBG)* &!
                (parZGy(i)-4.14D-1*parbigBGy(i)))/(ZG+2.414D0*bigBG)/(ZG-4.14D-1*bigBG)
            par33lx(i) = ((parZLx(i)+2.414D0*parbigBLx(i))*(ZL-4.14D-1*bigBL)-(ZL+2.414D0*bigBL)* &!
                (parZLx(i)-4.14D-1*parbigBLx(i)))/(ZL+2.414D0*bigBL)/(ZL-4.14D-1*bigBL)

            rightgy(i) = bm(i)/bg*(ZG-1.D0)-dlog((ZG-bigBG)*1.D0) - bigAG/(2.D0*sqrt(2.D0)*&!
                bigBG)*(tempsumg(i)-bm(i)/bg)*dlog((ZG+2.414D0*bigBG)/(ZG-4.14D-1*bigBG)*1.D0)
            rightlx(i) = bm(i)/bl*(ZL-1.D0)-dlog((ZL-bigBL)*1.D0) - bigAL/(2.D0*sqrt(2.D0)*&!
                bigBL)*(tempsuml(i)-bm(i)/bl)*dlog((ZL+2.414D0*bigBL)/(ZL-4.14D-1*bigBL)*1.D0)

        end do

        do j = 1, Nc
            do i = 1, Nc

                par1gy(i,j) = bm(i)/bg*parZGy(j) - bm(i)*bm(j)/bg**2.D0*(ZG-1.D0)
                par1lx(i,j) = bm(i)/bl*parZLx(j) - bm(i)*bm(j)/bl**2.D0*(ZL-1.D0)

                par32gy(i,j) = 2.D0/ag*(1.D0-delta(i,j))*dsqrt(am(i)*am(j)*1.D0) -paragy(j)/ &!
                    ag**2.D0*tempsumg(i) + bm(i)*bm(j)/bg**2.D0

                par32lx(i,j) = 2.D0/al*(1.D0-delta(i,j))*dsqrt(am(i)*am(j)*1.D0) -paralx(j)/ &!
                    al**2.D0*tempsuml(i) + bm(i)*bm(j)/bl**2.D0

                par3gy(i,j) = bigAG/(2.D0*sqrt(2.D0)*bigBG)*(tempsumg(i)/ag-bm(i)/bg)*par33gy(j) &!
                    + bigAG/(2.D0*sqrt(2.D0)*bigBG)*dlog((ZG+2.414D0*bigBG)/(ZG-4.14D-1*bigBG)*1.D0)* &!
                    par32gy(i,j)+(tempsumg(i)/ag-bm(i)/bg)*dlog((ZG+2.414D0*bigBG)/ &!
                    (ZG-4.14D-1*bigBG)*1.D0)*par31gy(j)

                par3lx(i,j) = bigAL/(2.D0*sqrt(2.D0)*bigBL)*(tempsuml(i)/al-bm(i)/bl)*par33lx(j) &!
                    + bigAL/(2.D0*sqrt(2.D0)*bigBL)*dlog((ZL+2.414D0*bigBL)/(ZL-4.14D-1*bigBL)*1.D0)* &!
                    par32lx(i,j)+(tempsuml(i)/al-bm(i)/bl)*dlog((ZL+2.414D0*bigBL)/ &!
                    (ZL-4.14D-1*bigBL)*1.D0)*par31lx(j)

                parfgy(i,j) = y(i)*P*dexp(rightgy(i)*1.D0)*(par1gy(i,j)-par2gy(j)- par3gy(i,j))
                parflx(i,j) = x(i)*P*dexp(rightlx(i)*1.D0)*(par1lx(i,j)-par2lx(j)- par3lx(i,j))

                if(i == j) then
                    parfgy(i,j) = parfgy(i,j) + P*dexp(rightgy(i)*1.D0)
                    parflx(i,j) = parflx(i,j) + P*dexp(rightlx(i)*1.D0)
                end if

                parfgn(i,j) = parfgy(i,j)/ng
                parfln(i,j) = parflx(i,j)/nl
            end do
        end do

        parZGp = parZGbigAG*ag/(R*Temp)**2.D0 + parZGbigBG*bg/R/Temp
        parZLp = parZLbigAL*al/(R*Temp)**2.D0 + parZLbigBL*bl/R/Temp

        parbigAGp = ag/(R*Temp)**2.D0
        parbigALp = al/(R*Temp)**2.D0

        parbigBGp = bg/R/Temp
        parbigBLp = bl/R/Temp

        par31gp = parbigAGp/(2.D0*sqrt(2.D0)*bigBG) - bigAG*parbigBGp/(2.D0*sqrt(2.D0)*bigBG**2.D0)
        par31lp = parbigALp/(2.D0*sqrt(2.D0)*bigBL) - bigAL*parbigBLp/(2.D0*sqrt(2.D0)*bigBL**2.D0)

        par33gp = ((parZGp-4.14D-1*parbigBGp)*(ZG+2.414D0*bigBG)-(ZG-4.14D-1*bigBG)* &!
            (parZGp+2.414D0*parbigBGp))/((ZG+2.414D0*bigBG)*(ZG-4.14D-1*bigBG))
        par33lp = ((parZLp-4.14D-1*parbigBLp)*(ZL+2.414D0*bigBL)-(ZL-4.14D-1*bigBL)* &!
            (parZLp+2.414D0*parbigBLp))/((ZL+2.414D0*bigBL)*(ZL-4.14D-1*bigBL))

        do i = 1, Nc

            par1gp(i) = bm(i)/bg*parZGp
            par1lp(i) = bm(i)/bl*parZLp

            par2gp(i) = 1.D0/(ZG-bigBG)*(parZGp-bg/R/Temp)
            par2lp(i) = 1.D0/(ZL-bigBL)*(parZLp-bl/R/Temp)

            par3gp(i) = bigAG/(2.D0*sqrt(2.D0)*bigBG)*(tempsumg(i)/ag-bm(i)/bg)*par33gp + &!
                (tempsumg(i)/ag-bm(i)/bg)*dlog((ZG+2.414D0*bigBG)/(ZG-4.14D-1*bigBG)*1.D0)*par31gp
            par3lp(i) = bigAL/(2.D0*sqrt(2.D0)*bigBL)*(tempsuml(i)/al-bm(i)/bl)*par33lp + &!
                (tempsuml(i)/al-bm(i)/bl)*dlog((ZL+2.414D0*bigBL)/(ZL-4.14D-1*bigBL)*1.D0)*par31lp

            parfgp(i) = y(i)*dexp(rightgy(i)*1.D0) + y(i)*P*dexp(rightgy(i)*1.D0)*(par1gp(i)-par2gp(i)-par3gp(i))
            parflp(i) = x(i)*dexp(rightlx(i)*1.D0) + x(i)*P*dexp(rightlx(i)*1.D0)*(par1lp(i)-par2lp(i)-par3lp(i))

        end do

        do i = 1, Nc
            do j = 1, Nc
                local_A(i,j) = parfgn(i,j) + parfln(i,j)
            end do
        end do

        do i = 1, Nc
            local_b(i) = parfgp(i) - parflp(i)
        end do

        !call dgesv(Nc, 1, local_A, Nc, IPIV, local_b, Nc, INFO) ! use lapack to solve AX=b
        call AGAUS(local_A, local_b, Nc, solx, L, JS)

        do i = 1, Nc
            parnlp(i) = solx(i)
            parngp(i) = -solx(i)
        end do

        do i = 1, Nc

            parbigAGn(i) = ag/(R*Temp)**2.D0/parngp(i)+P/(R*Temp)**2.D0/ng*paragy(i)
            parbigALn(i) = al/(R*Temp)**2.D0/parnlp(i)+P/(R*Temp)**2.D0/nl*paralx(i)

            parbigBGn(i) = bg/R/Temp/parngp(i)+P/R/Temp*bm(i)/ng
            parbigBLn(i) = bl/R/Temp/parnlp(i)+P/R/Temp*bm(i)/nl

            parZGn(i) = parZGbigAG*parbigAGn(i)+parZGbigBG*parbigBGn(i)
            parZLn(i) = parZLbigAL*parbigALn(i)+parZLbigBL*parbigBLn(i)

        end do

        parZGPtotaln = 0.D0
        do i = 1, Nc

            if(.not.(parZGn(i).ne.parZGn(i))) then
                parZGPtotaln = parZGPtotaln + parZGn(i)*parngp(i)
            end if
        end do
        parZGPtotaln = parZGPtotaln + parZGp

        parZLPtotaln = 0.D0
        do i = 1, Nc
            if(.not.(parZLn(i).ne.parZLn(i))) then
                parZLPtotaln = parZLPtotaln + parZLn(i)*parnlp(i)
            end if
        end do
        parZLPtotaln = parZLPtotaln + parZLp

        parngptotaln = 0.D0
        do i = 1, Nc
            parngptotaln = parngptotaln + parngp(i)
        end do

        parnlptotaln = 0.D0
        do i = 1, Nc
            parnlptotaln = parnlptotaln + parnlp(i)
        end do

        parvfp = -R*Temp/P**2.D0*(ZG*ng+ZL*nl)+R*Temp/P*(ng*parZGPtotaln+ZG*parngptotaln + &!
            nl*parZLPtotaln+ZL*parnlptotaln)

        ! since Vf=1, Cf = -parvfp/Vf = -parvfp
        local_Cf = -parvfp

        deallocate(local_A)
        deallocate(local_b)
        deallocate(solx)
        deallocate(JS)
        deallocate(tempsumg)
        deallocate(tempsuml)
        deallocate(paragy)
        deallocate(paralx)
        deallocate(parbigAGy)
        deallocate(parbigALx)
        deallocate(parbigBGy)
        deallocate(parbigBLx)
        deallocate(parZGy)
        deallocate(parZLx)
        deallocate(par2gy)
        deallocate(par2lx)
        deallocate(par31gy)
        deallocate(par31lx)
        deallocate(par33gy)
        deallocate(par33lx)
        deallocate(rightgy)
        deallocate(rightlx)
        deallocate(par1gy)
        deallocate(par1lx)
        deallocate(par32gy)
        deallocate(par32lx)
        deallocate(par3gy)
        deallocate(par3lx)
        deallocate(parfgy)
        deallocate(parflx)
        deallocate(parfgn)
        deallocate(parfln)
        deallocate(par1gp)
        deallocate(par1lp)
        deallocate(par2gp)
        deallocate(par2lp)
        deallocate(par3gp)
        deallocate(par3lp)
        deallocate(parfgp)
        deallocate(parflp)
        deallocate(parnlp)
        deallocate(parngp)
        deallocate(parbigAGn)
        deallocate(parbigALn)
        deallocate(parbigBGn)
        deallocate(parbigBLn)
        deallocate(parZGn)
        deallocate(parZLn)

    end function comprefac

end module RST_comprefac
