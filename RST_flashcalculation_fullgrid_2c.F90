
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_flashcalculation_fullgrid_2c

    use RST_globalFlashData

    implicit none

contains

    subroutine flashcalculation_fullgrid_2c(P, local_z, x, y, xiL, xiG, rhoL, rhoG, sL, local_v, local_Cf, isW, isN)

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
        
        real(kind=8) :: pmin
        real(kind=8) :: pmax
        real(kind=8) :: pinterval
        integer :: np
        real(kind=8) :: pleft, pright
        integer :: pleftindex
        real(kind=8) :: z1min
        real(kind=8) :: z1max
        real(kind=8) :: z1interval
        integer :: nz1
        real(kind=8) :: z1left, z1right
        integer :: z1leftindex

        real(kind=8), dimension(4) :: xinterpo
        real(kind=8), dimension(4) :: yinterpo
        real(kind=8), dimension(4) :: xiLinterpo
        real(kind=8), dimension(4) :: xiGinterpo
        real(kind=8), dimension(4) :: rhoLinterpo
        real(kind=8), dimension(4) :: rhoGinterpo
        real(kind=8), dimension(4) :: sLinterpo
        real(kind=8), dimension(4) :: vinterpo
        real(kind=8), dimension(4) :: Cfinterpo
        integer, dimension(4) :: interpoindex

        real(kind=8) :: buf1, buf2, dp, dz1, p_grid_coord, z1_grid_coord
        real(kind=8) :: sum
        integer :: num, i

        np = 129
        pmin = 1.9D6
        pmax = 2.1D6
        pinterval = (pmax-pmin)/(np-1)
        
        nz1 = 129
        z1min = 0.0
        z1max = 1.0
        z1interval = (z1max-z1min)/(nz1-1)

        p_grid_coord = (P-pmin)/pinterval+1
        z1_grid_coord = (local_z(1)-z1min)/z1interval+1

        pleftindex = floor(p_grid_coord)
        pleft = pmin + (pleftindex-1)*pinterval
        pright = pleft + pinterval

        z1leftindex = floor(z1_grid_coord)
        z1left = z1min + (z1leftindex-1)*z1interval
        z1right = z1left + z1interval

        dp = p_grid_coord - pleftindex
        dz1 = z1_grid_coord - z1leftindex

        interpoindex(1) = (pleftindex-1)*nz1 + z1leftindex
        if (z1leftindex /= nz1) then
            interpoindex(2) = interpoindex(1) + 1
        end if
        interpoindex(3) = interpoindex(1) + nz1
        if (z1leftindex /= nz1) then
            interpoindex(4) = interpoindex(3) + 1
        end if

        if (z1leftindex == nz1) then
            if (isWtable(interpoindex(1))*isWtable(interpoindex(3)) == 0) then
                isW = .false.
            else
                isW = .true.
            end if
            if (isNtable(interpoindex(1))*isNtable(interpoindex(3)) == 0) then
                isN = .false.
            else
                isN = .true.
            end if
        else
            if (isWtable(interpoindex(1))*isWtable(interpoindex(2))*isWtable(interpoindex(3))*isWtable(interpoindex(4)) == 0) then
                isW = .false.
            else
                isW = .true.
            end if
            if (isNtable(interpoindex(1))*isNtable(interpoindex(2))*isNtable(interpoindex(3))*isNtable(interpoindex(4)) == 0) then
                isN = .false.
            else
                isN = .true.
            end if
        end if

        !-------------------------x------------------
        if (.not.isW) then
            x(1) = 0
            x(2) = 0
        elseif (z1leftindex /= nz1) then
            xinterpo(1) = xtable(1,interpoindex(1))
            xinterpo(2) = xtable(1,interpoindex(2))
            xinterpo(3) = xtable(1,interpoindex(3))
            xinterpo(4) = xtable(1,interpoindex(4))

            buf1 = (1-dp)*xinterpo(1) + dp*xinterpo(3)
            buf2 = (1-dp)*xinterpo(2) + dp*xinterpo(4)
            x(1) = (1-dz1)*buf1 + dz1*buf2
            x(2) = 1 - x(1)
        else
            xinterpo(1) = xtable(1,interpoindex(1))
            xinterpo(3) = xtable(1,interpoindex(3))

            x(1) = (xinterpo(1)+xinterpo(3))/2
            x(2) = 1 - x(1)
        end if

        !-------------------y------------------------
        if (.not.isN) then
            y(1) = 0
            y(2) = 0
        elseif (z1leftindex /= nz1) then
            yinterpo(1) = ytable(1,interpoindex(1))
            yinterpo(2) = ytable(1,interpoindex(2))
            yinterpo(3) = ytable(1,interpoindex(3))
            yinterpo(4) = ytable(1,interpoindex(4))

            buf1 = (1-dp)*yinterpo(1) + dp*yinterpo(3)
            buf2 = (1-dp)*yinterpo(2) + dp*yinterpo(4)
            y(1) = (1-dz1)*buf1 + dz1*buf2
            y(2) = 1 - y(1)
        else
            yinterpo(1) = ytable(1,interpoindex(1))
            yinterpo(3) = ytable(1,interpoindex(3))
            y(1) = (yinterpo(1)+yinterpo(3))/2
            y(2) = 1 - y(1)
        end if

        !--------------------xiL---------------------
        if (.not.isW) then
            xiL = 0
        elseif (z1leftindex /= nz1) then
            xiLinterpo(1) = xiLtable(interpoindex(1))
            xiLinterpo(2) = xiLtable(interpoindex(2))
            xiLinterpo(3) = xiLtable(interpoindex(3))
            xiLinterpo(4) = xiLtable(interpoindex(4))

            buf1 = (1-dp)*xiLinterpo(1) + dp*xiLinterpo(3)
            buf2 = (1-dp)*xiLinterpo(2) + dp*xiLinterpo(4)
            xiL = (1-dz1)*buf1 + dz1*buf2
        else
            xiLinterpo(1) = xiLtable(interpoindex(1))
            xiLinterpo(3) = xiLtable(interpoindex(3))

            xiL = (xiLinterpo(1)+xiLinterpo(3))/2
        end if

        !-------------------xiG----------------------
        if (.not.isN) then
            xiG = 0
        elseif (z1leftindex /= nz1) then
            xiGinterpo(1) = xiGtable(interpoindex(1))
            xiGinterpo(2) = xiGtable(interpoindex(2))
            xiGinterpo(3) = xiGtable(interpoindex(3))
            xiGinterpo(4) = xiGtable(interpoindex(4))

            buf1 = (1-dp)*xiGinterpo(1) + dp*xiGinterpo(3)
            buf2 = (1-dp)*xiGinterpo(2) + dp*xiGinterpo(4)
            xiG = (1-dz1)*buf1 + dz1*buf2
        else
            xiGinterpo(1) = xiGtable(interpoindex(1))
            xiGinterpo(3) = xiGtable(interpoindex(3))

            xiG = (xiGinterpo(1)+xiGinterpo(3))/2
        end if

        !----------------rhoL--------------------------
        if (.not.isW) then
            rhoL = 0
        elseif (z1leftindex /= nz1) then
            rhoLinterpo(1) = rhoLtable(interpoindex(1))
            rhoLinterpo(2) = rhoLtable(interpoindex(2))
            rhoLinterpo(3) = rhoLtable(interpoindex(3))
            rhoLinterpo(4) = rhoLtable(interpoindex(4))

            buf1 = (1-dp)*rhoLinterpo(1) + dp*rhoLinterpo(3)
            buf2 = (1-dp)*rhoLinterpo(2) + dp*rhoLinterpo(4)
            rhoL = (1-dz1)*buf1 + dz1*buf2
        else
            rhoLinterpo(1) = rhoLtable(interpoindex(1))
            rhoLinterpo(3) = rhoLtable(interpoindex(3))

            rhoL = (rhoLinterpo(1)+rhoLinterpo(3))/2
        end if

        !-------------------rhoG-----------------------
        if (.not.isN) then
            rhoG = 0
        elseif (z1leftindex /= nz1) then
            rhoGinterpo(1) = rhoGtable(interpoindex(1))
            rhoGinterpo(2) = rhoGtable(interpoindex(2))
            rhoGinterpo(3) = rhoGtable(interpoindex(3))
            rhoGinterpo(4) = rhoGtable(interpoindex(4))

            buf1 = (1-dp)*rhoGinterpo(1) + dp*rhoGinterpo(3)
            buf2 = (1-dp)*rhoGinterpo(2) + dp*rhoGinterpo(4)
            rhoG = (1-dz1)*buf1 + dz1*buf2
        else
            rhoGinterpo(1) = rhoGtable(interpoindex(1))
            rhoGinterpo(3) = rhoGtable(interpoindex(3))

            rhoG = (rhoGinterpo(1)+rhoGinterpo(3))/2
        end if

        !------------------sL------------------
        if (.not.isW) then
            sL = 0
        elseif (.not.isN) then
            sL = 1
        elseif (z1leftindex /= nz1) then
            sLinterpo(1) = sLtable(interpoindex(1))
            sLinterpo(2) = sLtable(interpoindex(2))
            sLinterpo(3) = sLtable(interpoindex(3))
            sLinterpo(4) = sLtable(interpoindex(4))

            buf1 = (1-dp)*sLinterpo(1) + dp*sLinterpo(3)
            buf2 = (1-dp)*sLinterpo(2) + dp*sLinterpo(4)
            sL = (1-dz1)*buf1 + dz1*buf2
        else
            sLinterpo(1) = sLtable(interpoindex(1))
            sLinterpo(3) = sLtable(interpoindex(3))

            sL = (sLinterpo(1)+sLinterpo(3))/2
        end if

        !-------------------local_v-------------------
        do i = 1, 2
            if (.not.isW) then
                if (z1leftindex /= nz1) then
                    vinterpo(1) = vtable(i,interpoindex(1))
                    vinterpo(2) = vtable(i,interpoindex(2))
                    vinterpo(3) = vtable(i,interpoindex(3))
                    vinterpo(4) = vtable(i,interpoindex(4))

                    sum = 0
                    num = 0
                    if(isWtable(interpoindex(1)) == 0) then
                        sum = sum + vinterpo(1)
                        num = num + 1
                    end if
                    if(isWtable(interpoindex(2)) == 0) then
                        sum = sum + vinterpo(2)
                        num = num + 1
                    end if
                    if(isWtable(interpoindex(3)) == 0) then
                        sum = sum + vinterpo(3)
                        num = num + 1
                    end if
                    if(isWtable(interpoindex(4)) == 0) then
                        sum = sum + vinterpo(4)
                        num = num + 1
                    end if
                    local_v(i) = sum/num
                else
                    vinterpo(1) = vtable(i,interpoindex(1))
                    vinterpo(3) = vtable(i,interpoindex(3))

                    sum = 0
                    num = 0
                    if(isWtable(interpoindex(1)) == 0) then
                        sum = sum + vinterpo(1)
                        num = num + 1
                    end if
                    if(isWtable(interpoindex(3)) == 0) then
                        sum = sum + vinterpo(3)
                        num = num + 1
                    end if
                    local_v(i) = sum/num
                end if
            elseif (.not.isN) then
                if (z1leftindex /= nz1) then
                    vinterpo(1) = vtable(i,interpoindex(1))
                    vinterpo(2) = vtable(i,interpoindex(2))
                    vinterpo(3) = vtable(i,interpoindex(3))
                    vinterpo(4) = vtable(i,interpoindex(4))

                    sum = 0
                    num = 0
                    if(isNtable(interpoindex(1)) == 0) then
                        sum = sum + vinterpo(1)
                        num = num + 1
                    end if
                    if(isNtable(interpoindex(2)) == 0) then
                        sum = sum + vinterpo(2)
                        num = num + 1
                    end if
                    if(isNtable(interpoindex(3)) == 0) then
                        sum = sum + vinterpo(3)
                        num = num + 1
                    end if
                    if(isNtable(interpoindex(4)) == 0) then
                        sum = sum + vinterpo(4)
                        num = num + 1
                    end if
                    local_v(i) = sum/num
                else
                    vinterpo(1) = vtable(i,interpoindex(1))
                    vinterpo(3) = vtable(i,interpoindex(3))

                    sum = 0
                    num = 0
                    if(isNtable(interpoindex(1)) == 0) then
                        sum = sum + vinterpo(1)
                        num = num + 1
                    end if
                    if(isNtable(interpoindex(3)) == 0) then
                        sum = sum + vinterpo(3)
                        num = num + 1
                    end if
                    local_v(i) = sum/num
                end if
            else
                if (z1leftindex /= nz1) then
                    vinterpo(1) = vtable(i,interpoindex(1))
                    vinterpo(2) = vtable(i,interpoindex(2))
                    vinterpo(3) = vtable(i,interpoindex(3))
                    vinterpo(4) = vtable(i,interpoindex(4))

                    buf1 = (1-dp)*vinterpo(1) + dp*vinterpo(3)
                    buf2 = (1-dp)*vinterpo(2) + dp*vinterpo(4)
                    local_v(i) = (1-dz1)*buf1 + dz1*buf2
                else
                    vinterpo(1) = vtable(i,interpoindex(1))
                    vinterpo(3) = vtable(i,interpoindex(3))

                    local_v(i) = (vinterpo(1)+vinterpo(3))/2
                end if
            end if
            if(abs(local_z(i)) < 1.D-12) then
                local_v(i) = 0.0
            end if
        end do

        !------------local_Cf------------------------
        if (.not.isW) then
            if (z1leftindex /= nz1) then
                Cfinterpo(1) = Cftable(interpoindex(1))
                Cfinterpo(2) = Cftable(interpoindex(2))
                Cfinterpo(3) = Cftable(interpoindex(3))
                Cfinterpo(4) = Cftable(interpoindex(4))

                sum = 0
                num = 0
                if(isWtable(interpoindex(1)) == 0) then
                    sum = sum + Cfinterpo(1)
                    num = num + 1
                end if
                if(isWtable(interpoindex(2)) == 0) then
                    sum = sum + Cfinterpo(2)
                    num = num + 1
                end if
                if(isWtable(interpoindex(3)) == 0) then
                    sum = sum + Cfinterpo(3)
                    num = num + 1
                end if
                if(isWtable(interpoindex(4)) == 0) then
                    sum = sum + Cfinterpo(4)
                    num = num + 1
                end if
                local_Cf = sum/num
            else
                Cfinterpo(1) = Cftable(interpoindex(1))
                Cfinterpo(3) = Cftable(interpoindex(3))

                sum = 0
                num = 0
                if(isWtable(interpoindex(1)) == 0) then
                    sum = sum + Cfinterpo(1)
                    num = num + 1
                end if
                if(isWtable(interpoindex(3)) == 0) then
                    sum = sum + Cfinterpo(3)
                    num = num + 1
                end if
                local_Cf = sum/num
            end if
        elseif (.not.isN) then
            if (z1leftindex /= nz1) then
                Cfinterpo(1) = Cftable(interpoindex(1))
                Cfinterpo(2) = Cftable(interpoindex(2))
                Cfinterpo(3) = Cftable(interpoindex(3))
                Cfinterpo(4) = Cftable(interpoindex(4))

                sum = 0
                num = 0
                if(isNtable(interpoindex(1)) == 0) then
                    sum = sum + Cfinterpo(1)
                    num = num + 1
                end if
                if(isNtable(interpoindex(2)) == 0) then
                    sum = sum + Cfinterpo(2)
                    num = num + 1
                end if
                if(isNtable(interpoindex(3)) == 0) then
                    sum = sum + Cfinterpo(3)
                    num = num + 1
                end if
                if(isNtable(interpoindex(4)) == 0) then
                    sum = sum + Cfinterpo(4)
                    num = num + 1
                end if
                local_Cf = sum/num
            else
                Cfinterpo(1) = Cftable(interpoindex(1))
                Cfinterpo(3) = Cftable(interpoindex(3))

                sum = 0
                num = 0
                if(isNtable(interpoindex(1)) == 0) then
                    sum = sum + Cfinterpo(1)
                    num = num + 1
                end if
                if(isNtable(interpoindex(3)) == 0) then
                    sum = sum + Cfinterpo(3)
                    num = num + 1
                end if
                local_Cf = sum/num
            end if
        else
            if (z1leftindex /= nz1) then
                Cfinterpo(1) = Cftable(interpoindex(1))
                Cfinterpo(2) = Cftable(interpoindex(2))
                Cfinterpo(3) = Cftable(interpoindex(3))
                Cfinterpo(4) = Cftable(interpoindex(4))

                buf1 = (1-dp)*Cfinterpo(1) + dp*Cfinterpo(3)
                buf2 = (1-dp)*Cfinterpo(2) + dp*Cfinterpo(4)
                local_Cf = (1-dz1)*buf1 + dz1*buf2
            else
                Cfinterpo(1) = Cftable(interpoindex(1))
                Cfinterpo(3) = Cftable(interpoindex(3))

                local_Cf = (Cfinterpo(1)+Cfinterpo(3))/2
            end if
        end if

    end subroutine flashcalculation_fullgrid_2c

end module RST_flashcalculation_fullgrid_2c
