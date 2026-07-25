
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_flashcalculation_fullgrid_3c

    use RST_globalFlashData

    implicit none

contains

    subroutine flashcalculation_fullgrid_3c(P, local_z, x, y, xiL, xiG, rhoL, rhoG, sL, local_v, local_Cf, isW, isN)

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
        real(kind=8) :: z2min
        real(kind=8) :: z2max
        real(kind=8) :: z2interval
        integer :: nz2
        real(kind=8) :: z2left, z2right
        integer :: z2leftindex

        real(kind=8), dimension(8) :: xinterpo
        real(kind=8), dimension(8) :: yinterpo
        real(kind=8), dimension(8) :: xiLinterpo
        real(kind=8), dimension(8) :: xiGinterpo
        real(kind=8), dimension(8) :: rhoLinterpo
        real(kind=8), dimension(8) :: rhoGinterpo
        real(kind=8), dimension(8) :: sLinterpo
        real(kind=8), dimension(8) :: vinterpo
        real(kind=8), dimension(8) :: Cfinterpo
        integer, dimension(8) :: interpoindex

        real(kind=8) :: buf1, buf2, buf3, buf4, buf5, buf6, dp, dz1, &!
            dz2, p_grid_coord, z1_grid_coord, z2_grid_coord
        real(kind=8) :: sum
        integer :: num, i

        np = 17!129
        pmin = 1.9D6!1.984D6
        pmax = 1.6D8!2.112D6
        pinterval = (pmax-pmin)/(np-1)
        
        nz1 = 257!129
        z1min = 0.0
        z1max = 1.0
        z1interval = (z1max-z1min)/(nz1-1)

        nz2 = 257!129
        z2min = 0.0
        z2max = 1.0
        z2interval = (z2max-z2min)/(nz2-1)

        p_grid_coord = (P-pmin)/pinterval+1
        z1_grid_coord = (local_z(1)-z1min)/z1interval+1
        z2_grid_coord = (local_z(2)-z2min)/z2interval+1

        pleftindex = floor(p_grid_coord)
        pleft = pmin + (pleftindex-1)*pinterval
        pright = pleft + pinterval

        z1leftindex = floor(z1_grid_coord)
        z1left = z1min + (z1leftindex-1)*z1interval
        z1right = z1left + z1interval

        z2leftindex = floor(z2_grid_coord)
        z2left = z2min + (z2leftindex-1)*z2interval
        z2right = z2left + z2interval

        dp = p_grid_coord - pleftindex
        dz1 = z1_grid_coord - z1leftindex
        dz2 = z2_grid_coord - z2leftindex

        ! half of the column
        interpoindex(1) = (pleftindex-1)*(1+nz1)*nz2/2 + (nz1+nz2-z1leftindex+2)*(z1leftindex-1)/2 + z2leftindex
        interpoindex(2) = interpoindex(1) + 1
        interpoindex(3) = interpoindex(1) + (nz1-(z1leftindex-1))
        if(z1right+z2right <= 1) then
            interpoindex(4) = interpoindex(3) + 1
        else
            interpoindex(4) = -1
        end if
        
        interpoindex(5) = interpoindex(1) + (1+nz1)*nz2/2
        interpoindex(6) = interpoindex(5) + 1
        interpoindex(7) = interpoindex(5) + (nz1-(z1leftindex-1))
        if(z1right+z2right <= 1) then
            interpoindex(8) = interpoindex(7) + 1
        else
            interpoindex(8) = -1
        end if

        if (interpoindex(4) == -1) then
            if (isWtable(interpoindex(1))*isWtable(interpoindex(2))*isWtable(interpoindex(3))*isWtable(interpoindex(5))* &!
                isWtable(interpoindex(6))*isWtable(interpoindex(7)) == 0) then
                isW = .false.
            else
                isW = .true.
            end if
            if (isNtable(interpoindex(1))*isNtable(interpoindex(2))*isNtable(interpoindex(3))*isNtable(interpoindex(5))* &!
                isNtable(interpoindex(6))*isNtable(interpoindex(7)) == 0) then
                isN = .false.
            else
                isN = .true.
            end if
        else
            if (isWtable(interpoindex(1))*isWtable(interpoindex(2))*isWtable(interpoindex(3))*isWtable(interpoindex(4))* &!
                isWtable(interpoindex(5))*isWtable(interpoindex(6))*isWtable(interpoindex(7))*isWtable(interpoindex(8)) == 0) then
                isW = .false.
            else
                isW = .true.
            end if
            if (isNtable(interpoindex(1))*isNtable(interpoindex(2))*isNtable(interpoindex(3))*isNtable(interpoindex(4))* &!
                isNtable(interpoindex(5))*isNtable(interpoindex(6))*isNtable(interpoindex(7))*isNtable(interpoindex(8)) == 0) then
                isN = .false.
            else
                isN = .true.
            end if
        end if

        !-------------------------x------------------
        if (.not.isW) then
            x(1) = 0
            x(2) = 0
            x(3) = 0
        elseif(interpoindex(4) /= -1) then
            xinterpo(1) = xtable(1,interpoindex(1))
            xinterpo(2) = xtable(1,interpoindex(2))
            xinterpo(3) = xtable(1,interpoindex(3))
            xinterpo(4) = xtable(1,interpoindex(4))
            xinterpo(5) = xtable(1,interpoindex(5))
            xinterpo(6) = xtable(1,interpoindex(6))
            xinterpo(7) = xtable(1,interpoindex(7))
            xinterpo(8) = xtable(1,interpoindex(8))

            buf1 = (1-dp)*xinterpo(1) + dp*xinterpo(5)
            buf2 = (1-dp)*xinterpo(3) + dp*xinterpo(7)
            buf3 = (1-dp)*xinterpo(2) + dp*xinterpo(6)
            buf4 = (1-dp)*xinterpo(4) + dp*xinterpo(8)

            buf5 = (1-dz1)*buf1 + dz1*buf2
            buf6 = (1-dz1)*buf3 + dz1*buf4

            x(1) = (1-dz2)*buf5 + dz2*buf6

            xinterpo(1) = xtable(2,interpoindex(1))
            xinterpo(2) = xtable(2,interpoindex(2))
            xinterpo(3) = xtable(2,interpoindex(3))
            xinterpo(4) = xtable(2,interpoindex(4))
            xinterpo(5) = xtable(2,interpoindex(5))
            xinterpo(6) = xtable(2,interpoindex(6))
            xinterpo(7) = xtable(2,interpoindex(7))
            xinterpo(8) = xtable(2,interpoindex(8))

            buf1 = (1-dp)*xinterpo(1) + dp*xinterpo(5)
            buf2 = (1-dp)*xinterpo(3) + dp*xinterpo(7)
            buf3 = (1-dp)*xinterpo(2) + dp*xinterpo(6)
            buf4 = (1-dp)*xinterpo(4) + dp*xinterpo(8)

            buf5 = (1-dz1)*buf1 + dz1*buf2
            buf6 = (1-dz1)*buf3 + dz1*buf4

            x(2) = (1-dz2)*buf5 + dz2*buf6

            x(3) = 1 - x(1) - x(2)
        else
            xinterpo(1) = xtable(1,interpoindex(1))
            xinterpo(2) = xtable(1,interpoindex(2))
            xinterpo(3) = xtable(1,interpoindex(3))
            xinterpo(5) = xtable(1,interpoindex(5))
            xinterpo(6) = xtable(1,interpoindex(6))
            xinterpo(7) = xtable(1,interpoindex(7))

            sum = xinterpo(1)+xinterpo(2)+xinterpo(3)+xinterpo(5)+xinterpo(6)+xinterpo(7)
            x(1) = sum/6

            xinterpo(1) = xtable(2,interpoindex(1))
            xinterpo(2) = xtable(2,interpoindex(2))
            xinterpo(3) = xtable(2,interpoindex(3))
            xinterpo(5) = xtable(2,interpoindex(5))
            xinterpo(6) = xtable(2,interpoindex(6))
            xinterpo(7) = xtable(2,interpoindex(7))

            sum = xinterpo(1)+xinterpo(2)+xinterpo(3)+xinterpo(5)+xinterpo(6)+xinterpo(7)
            x(2) = sum/6

            x(3) = 1 - x(1) - x(2)

        end if

        !-------------------y------------------------
        if (.not.isN) then
            y(1) = 0
            y(2) = 0
            y(3) = 0
        elseif(interpoindex(4) /= -1) then
            yinterpo(1) = ytable(1,interpoindex(1))
            yinterpo(2) = ytable(1,interpoindex(2))
            yinterpo(3) = ytable(1,interpoindex(3))
            yinterpo(4) = ytable(1,interpoindex(4))
            yinterpo(5) = ytable(1,interpoindex(5))
            yinterpo(6) = ytable(1,interpoindex(6))
            yinterpo(7) = ytable(1,interpoindex(7))
            yinterpo(8) = ytable(1,interpoindex(8))

            buf1 = (1-dp)*yinterpo(1) + dp*yinterpo(5)
            buf2 = (1-dp)*yinterpo(3) + dp*yinterpo(7)
            buf3 = (1-dp)*yinterpo(2) + dp*yinterpo(6)
            buf4 = (1-dp)*yinterpo(4) + dp*yinterpo(8)

            buf5 = (1-dz1)*buf1 + dz1*buf2
            buf6 = (1-dz1)*buf3 + dz1*buf4

            y(1) = (1-dz2)*buf5 + dz2*buf6

            yinterpo(1) = ytable(2,interpoindex(1))
            yinterpo(2) = ytable(2,interpoindex(2))
            yinterpo(3) = ytable(2,interpoindex(3))
            yinterpo(4) = ytable(2,interpoindex(4))
            yinterpo(5) = ytable(2,interpoindex(5))
            yinterpo(6) = ytable(2,interpoindex(6))
            yinterpo(7) = ytable(2,interpoindex(7))
            yinterpo(8) = ytable(2,interpoindex(8))

            buf1 = (1-dp)*yinterpo(1) + dp*yinterpo(5)
            buf2 = (1-dp)*yinterpo(3) + dp*yinterpo(7)
            buf3 = (1-dp)*yinterpo(2) + dp*yinterpo(6)
            buf4 = (1-dp)*yinterpo(4) + dp*yinterpo(8)

            buf5 = (1-dz1)*buf1 + dz1*buf2
            buf6 = (1-dz1)*buf3 + dz1*buf4

            y(2) = (1-dz2)*buf5 + dz2*buf6

            y(3) = 1 - y(1) - y(2)
        else
            yinterpo(1) = ytable(1,interpoindex(1))
            yinterpo(2) = ytable(1,interpoindex(2))
            yinterpo(3) = ytable(1,interpoindex(3))
            yinterpo(5) = ytable(1,interpoindex(5))
            yinterpo(6) = ytable(1,interpoindex(6))
            yinterpo(7) = ytable(1,interpoindex(7))

            sum = yinterpo(1)+yinterpo(2)+yinterpo(3)+yinterpo(5)+yinterpo(6)+yinterpo(7)
            y(1) = sum/6

            yinterpo(1) = ytable(2,interpoindex(1))
            yinterpo(2) = ytable(2,interpoindex(2))
            yinterpo(3) = ytable(2,interpoindex(3))
            yinterpo(5) = ytable(2,interpoindex(5))
            yinterpo(6) = ytable(2,interpoindex(6))
            yinterpo(7) = ytable(2,interpoindex(7))

            sum = yinterpo(1)+yinterpo(2)+yinterpo(3)+yinterpo(5)+yinterpo(6)+yinterpo(7)
            y(2) = sum/6

            y(3) = 1 - y(1) - y(2)

        end if

        !--------------------xiL---------------------
        if (.not.isW) then
            xiL = 0
        elseif(interpoindex(4) /= -1) then
            xiLinterpo(1) = xiLtable(interpoindex(1))
            xiLinterpo(2) = xiLtable(interpoindex(2))
            xiLinterpo(3) = xiLtable(interpoindex(3))
            xiLinterpo(4) = xiLtable(interpoindex(4))
            xiLinterpo(5) = xiLtable(interpoindex(5))
            xiLinterpo(6) = xiLtable(interpoindex(6))
            xiLinterpo(7) = xiLtable(interpoindex(7))
            xiLinterpo(8) = xiLtable(interpoindex(8))

            buf1 = (1-dp)*xiLinterpo(1) + dp*xiLinterpo(5)
            buf2 = (1-dp)*xiLinterpo(3) + dp*xiLinterpo(7)
            buf3 = (1-dp)*xiLinterpo(2) + dp*xiLinterpo(6)
            buf4 = (1-dp)*xiLinterpo(4) + dp*xiLinterpo(8)

            buf5 = (1-dz1)*buf1 + dz1*buf2
            buf6 = (1-dz1)*buf3 + dz1*buf4

            xiL = (1-dz2)*buf5 + dz2*buf6
        else
            xiLinterpo(1) = xiLtable(interpoindex(1))
            xiLinterpo(2) = xiLtable(interpoindex(2))
            xiLinterpo(3) = xiLtable(interpoindex(3))
            xiLinterpo(5) = xiLtable(interpoindex(5))
            xiLinterpo(6) = xiLtable(interpoindex(6))
            xiLinterpo(7) = xiLtable(interpoindex(7))

            sum = xiLinterpo(1)+xiLinterpo(2)+xiLinterpo(3)+xiLinterpo(5)+xiLinterpo(6)+xiLinterpo(7)
            xiL = sum/6
        end if

        !-------------------xiG----------------------
        if (.not.isN) then
            xiG = 0
        elseif(interpoindex(4) /= -1) then
            xiGinterpo(1) = xiGtable(interpoindex(1))
            xiGinterpo(2) = xiGtable(interpoindex(2))
            xiGinterpo(3) = xiGtable(interpoindex(3))
            xiGinterpo(4) = xiGtable(interpoindex(4))
            xiGinterpo(5) = xiGtable(interpoindex(5))
            xiGinterpo(6) = xiGtable(interpoindex(6))
            xiGinterpo(7) = xiGtable(interpoindex(7))
            xiGinterpo(8) = xiGtable(interpoindex(8))

            buf1 = (1-dp)*xiGinterpo(1) + dp*xiGinterpo(5)
            buf2 = (1-dp)*xiGinterpo(3) + dp*xiGinterpo(7)
            buf3 = (1-dp)*xiGinterpo(2) + dp*xiGinterpo(6)
            buf4 = (1-dp)*xiGinterpo(4) + dp*xiGinterpo(8)

            buf5 = (1-dz1)*buf1 + dz1*buf2
            buf6 = (1-dz1)*buf3 + dz1*buf4

            xiG = (1-dz2)*buf5 + dz2*buf6
        else
            xiGinterpo(1) = xiGtable(interpoindex(1))
            xiGinterpo(2) = xiGtable(interpoindex(2))
            xiGinterpo(3) = xiGtable(interpoindex(3))
            xiGinterpo(5) = xiGtable(interpoindex(5))
            xiGinterpo(6) = xiGtable(interpoindex(6))
            xiGinterpo(7) = xiGtable(interpoindex(7))

            sum = xiGinterpo(1)+xiGinterpo(2)+xiGinterpo(3)+xiGinterpo(5)+xiGinterpo(6)+xiGinterpo(7)
            xiG = sum/6
        end if

        !----------------rhoL--------------------------
        if (.not.isW) then
            rhoL = 0
        elseif(interpoindex(4) /= -1) then
            rhoLinterpo(1) = rhoLtable(interpoindex(1))
            rhoLinterpo(2) = rhoLtable(interpoindex(2))
            rhoLinterpo(3) = rhoLtable(interpoindex(3))
            rhoLinterpo(4) = rhoLtable(interpoindex(4))
            rhoLinterpo(5) = rhoLtable(interpoindex(5))
            rhoLinterpo(6) = rhoLtable(interpoindex(6))
            rhoLinterpo(7) = rhoLtable(interpoindex(7))
            rhoLinterpo(8) = rhoLtable(interpoindex(8))

            buf1 = (1-dp)*rhoLinterpo(1) + dp*rhoLinterpo(5)
            buf2 = (1-dp)*rhoLinterpo(3) + dp*rhoLinterpo(7)
            buf3 = (1-dp)*rhoLinterpo(2) + dp*rhoLinterpo(6)
            buf4 = (1-dp)*rhoLinterpo(4) + dp*rhoLinterpo(8)

            buf5 = (1-dz1)*buf1 + dz1*buf2
            buf6 = (1-dz1)*buf3 + dz1*buf4

            rhoL = (1-dz2)*buf5 + dz2*buf6
        else
            rhoLinterpo(1) = rhoLtable(interpoindex(1))
            rhoLinterpo(2) = rhoLtable(interpoindex(2))
            rhoLinterpo(3) = rhoLtable(interpoindex(3))
            rhoLinterpo(5) = rhoLtable(interpoindex(5))
            rhoLinterpo(6) = rhoLtable(interpoindex(6))
            rhoLinterpo(7) = rhoLtable(interpoindex(7))

            sum = rhoLinterpo(1)+rhoLinterpo(2)+rhoLinterpo(3)+rhoLinterpo(5)+rhoLinterpo(6)+rhoLinterpo(7)
            rhoL = sum/6
        end if

        !-------------------rhoG-----------------------
        if (.not.isN) then
            rhoG = 0
        elseif(interpoindex(4) /= -1) then
            rhoGinterpo(1) = rhoGtable(interpoindex(1))
            rhoGinterpo(2) = rhoGtable(interpoindex(2))
            rhoGinterpo(3) = rhoGtable(interpoindex(3))
            rhoGinterpo(4) = rhoGtable(interpoindex(4))
            rhoGinterpo(5) = rhoGtable(interpoindex(5))
            rhoGinterpo(6) = rhoGtable(interpoindex(6))
            rhoGinterpo(7) = rhoGtable(interpoindex(7))
            rhoGinterpo(8) = rhoGtable(interpoindex(8))

            buf1 = (1-dp)*rhoGinterpo(1) + dp*rhoGinterpo(5)
            buf2 = (1-dp)*rhoGinterpo(3) + dp*rhoGinterpo(7)
            buf3 = (1-dp)*rhoGinterpo(2) + dp*rhoGinterpo(6)
            buf4 = (1-dp)*rhoGinterpo(4) + dp*rhoGinterpo(8)

            buf5 = (1-dz1)*buf1 + dz1*buf2
            buf6 = (1-dz1)*buf3 + dz1*buf4

            rhoG = (1-dz2)*buf5 + dz2*buf6
        else
            rhoGinterpo(1) = rhoGtable(interpoindex(1))
            rhoGinterpo(2) = rhoGtable(interpoindex(2))
            rhoGinterpo(3) = rhoGtable(interpoindex(3))
            rhoGinterpo(5) = rhoGtable(interpoindex(5))
            rhoGinterpo(6) = rhoGtable(interpoindex(6))
            rhoGinterpo(7) = rhoGtable(interpoindex(7))

            sum = rhoGinterpo(1)+rhoGinterpo(2)+rhoGinterpo(3)+rhoGinterpo(5)+rhoGinterpo(6)+rhoGinterpo(7)
            rhoG = sum/6
        end if

        !------------------sL------------------
        if (.not.isW) then
            sL = 0
        elseif (.not.isN) then
            sL = 1
        elseif(interpoindex(4) /= -1) then
            sLinterpo(1) = sLtable(interpoindex(1))
            sLinterpo(2) = sLtable(interpoindex(2))
            sLinterpo(3) = sLtable(interpoindex(3))
            sLinterpo(4) = sLtable(interpoindex(4))
            sLinterpo(5) = sLtable(interpoindex(5))
            sLinterpo(6) = sLtable(interpoindex(6))
            sLinterpo(7) = sLtable(interpoindex(7))
            sLinterpo(8) = sLtable(interpoindex(8))

            buf1 = (1-dp)*sLinterpo(1) + dp*sLinterpo(5)
            buf2 = (1-dp)*sLinterpo(3) + dp*sLinterpo(7)
            buf3 = (1-dp)*sLinterpo(2) + dp*sLinterpo(6)
            buf4 = (1-dp)*sLinterpo(4) + dp*sLinterpo(8)

            buf5 = (1-dz1)*buf1 + dz1*buf2
            buf6 = (1-dz1)*buf3 + dz1*buf4

            sL = (1-dz2)*buf5 + dz2*buf6
        else
            sLinterpo(1) = sLtable(interpoindex(1))
            sLinterpo(2) = sLtable(interpoindex(2))
            sLinterpo(3) = sLtable(interpoindex(3))
            sLinterpo(5) = sLtable(interpoindex(5))
            sLinterpo(6) = sLtable(interpoindex(6))
            sLinterpo(7) = sLtable(interpoindex(7))

            sum = sLinterpo(1)+sLinterpo(2)+sLinterpo(3)+sLinterpo(5)+sLinterpo(6)+sLinterpo(7)
            sL = sum/6
        end if

        !-------------------local_v-------------------
        do i = 1, 3
            if (.not.isW) then
                if (interpoindex(4) /= -1) then
                    vinterpo(1) = vtable(i,interpoindex(1))
                    vinterpo(2) = vtable(i,interpoindex(2))
                    vinterpo(3) = vtable(i,interpoindex(3))
                    vinterpo(4) = vtable(i,interpoindex(4))
                    vinterpo(5) = vtable(i,interpoindex(5))
                    vinterpo(6) = vtable(i,interpoindex(6))
                    vinterpo(7) = vtable(i,interpoindex(7))
                    vinterpo(8) = vtable(i,interpoindex(8))

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
                    if(isWtable(interpoindex(5)) == 0) then
                        sum = sum + vinterpo(5)
                        num = num + 1
                    end if
                    if(isWtable(interpoindex(6)) == 0) then
                        sum = sum + vinterpo(6)
                        num = num + 1
                    end if
                    if(isWtable(interpoindex(7)) == 0) then
                        sum = sum + vinterpo(7)
                        num = num + 1
                    end if
                    if(isWtable(interpoindex(8)) == 0) then
                        sum = sum + vinterpo(8)
                        num = num + 1
                    end if
                    local_v(i) = sum/num
                else
                    vinterpo(1) = vtable(i,interpoindex(1))
                    vinterpo(2) = vtable(i,interpoindex(2))
                    vinterpo(3) = vtable(i,interpoindex(3))
                    vinterpo(5) = vtable(i,interpoindex(5))
                    vinterpo(6) = vtable(i,interpoindex(6))
                    vinterpo(7) = vtable(i,interpoindex(7))

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
                    if(isWtable(interpoindex(5)) == 0) then
                        sum = sum + vinterpo(5)
                        num = num + 1
                    end if
                    if(isWtable(interpoindex(6)) == 0) then
                        sum = sum + vinterpo(6)
                        num = num + 1
                    end if
                    if(isWtable(interpoindex(7)) == 0) then
                        sum = sum + vinterpo(7)
                        num = num + 1
                    end if
                    local_v(i) = sum/num
                end if
            elseif (.not.isN) then
                if (interpoindex(4) /= -1) then
                    vinterpo(1) = vtable(i,interpoindex(1))
                    vinterpo(2) = vtable(i,interpoindex(2))
                    vinterpo(3) = vtable(i,interpoindex(3))
                    vinterpo(4) = vtable(i,interpoindex(4))
                    vinterpo(5) = vtable(i,interpoindex(5))
                    vinterpo(6) = vtable(i,interpoindex(6))
                    vinterpo(7) = vtable(i,interpoindex(7))
                    vinterpo(8) = vtable(i,interpoindex(8))

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
                    if(isNtable(interpoindex(5)) == 0) then
                        sum = sum + vinterpo(5)
                        num = num + 1
                    end if
                    if(isNtable(interpoindex(6)) == 0) then
                        sum = sum + vinterpo(6)
                        num = num + 1
                    end if
                    if(isNtable(interpoindex(7)) == 0) then
                        sum = sum + vinterpo(7)
                        num = num + 1
                    end if
                    if(isNtable(interpoindex(8)) == 0) then
                        sum = sum + vinterpo(8)
                        num = num + 1
                    end if
                    local_v(i) = sum/num
                else
                    vinterpo(1) = vtable(i,interpoindex(1))
                    vinterpo(2) = vtable(i,interpoindex(2))
                    vinterpo(3) = vtable(i,interpoindex(3))
                    vinterpo(5) = vtable(i,interpoindex(5))
                    vinterpo(6) = vtable(i,interpoindex(6))
                    vinterpo(7) = vtable(i,interpoindex(7))

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
                    if(isNtable(interpoindex(5)) == 0) then
                        sum = sum + vinterpo(5)
                        num = num + 1
                    end if
                    if(isNtable(interpoindex(6)) == 0) then
                        sum = sum + vinterpo(6)
                        num = num + 1
                    end if
                    if(isNtable(interpoindex(7)) == 0) then
                        sum = sum + vinterpo(7)
                        num = num + 1
                    end if
                    local_v(i) = sum/num
                end if
            else
                if (interpoindex(4) /= -1) then
                    vinterpo(1) = vtable(i,interpoindex(1))
                    vinterpo(2) = vtable(i,interpoindex(2))
                    vinterpo(3) = vtable(i,interpoindex(3))
                    vinterpo(4) = vtable(i,interpoindex(4))
                    vinterpo(5) = vtable(i,interpoindex(5))
                    vinterpo(6) = vtable(i,interpoindex(6))
                    vinterpo(7) = vtable(i,interpoindex(7))
                    vinterpo(8) = vtable(i,interpoindex(8))

                    buf1 = (1-dp)*vinterpo(1) + dp*vinterpo(5)
                    buf2 = (1-dp)*vinterpo(3) + dp*vinterpo(7)
                    buf3 = (1-dp)*vinterpo(2) + dp*vinterpo(6)
                    buf4 = (1-dp)*vinterpo(4) + dp*vinterpo(8)

                    buf5 = (1-dz1)*buf1 + dz1*buf2
                    buf6 = (1-dz1)*buf3 + dz1*buf4

                    local_v(i) = (1-dz2)*buf5 + dz2*buf6
                else
                    vinterpo(1) = vtable(i,interpoindex(1))
                    vinterpo(2) = vtable(i,interpoindex(2))
                    vinterpo(3) = vtable(i,interpoindex(3))
                    vinterpo(5) = vtable(i,interpoindex(5))
                    vinterpo(6) = vtable(i,interpoindex(6))
                    vinterpo(7) = vtable(i,interpoindex(7))

                    sum = vinterpo(1)+vinterpo(2)+vinterpo(3)+vinterpo(5)+vinterpo(6)+vinterpo(7)
                    local_v(i) = sum/6
                end if
            end if
            if(abs(local_z(i)) < 1.D-12) then
                local_v(i) = 0.0
            end if
        end do

        !------------local_Cf------------------------
        if (.not.isW) then
            if (interpoindex(4) /= -1) then
                Cfinterpo(1) = Cftable(interpoindex(1))
                Cfinterpo(2) = Cftable(interpoindex(2))
                Cfinterpo(3) = Cftable(interpoindex(3))
                Cfinterpo(4) = Cftable(interpoindex(4))
                Cfinterpo(5) = Cftable(interpoindex(5))
                Cfinterpo(6) = Cftable(interpoindex(6))
                Cfinterpo(7) = Cftable(interpoindex(7))
                Cfinterpo(8) = Cftable(interpoindex(8))

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
                if(isWtable(interpoindex(5)) == 0) then
                    sum = sum + Cfinterpo(5)
                    num = num + 1
                end if
                if(isWtable(interpoindex(6)) == 0) then
                    sum = sum + Cfinterpo(6)
                    num = num + 1
                end if
                if(isWtable(interpoindex(7)) == 0) then
                    sum = sum + Cfinterpo(7)
                    num = num + 1
                end if
                if(isWtable(interpoindex(8)) == 0) then
                    sum = sum + Cfinterpo(8)
                    num = num + 1
                end if
                local_Cf = sum/num
            else
                Cfinterpo(1) = Cftable(interpoindex(1))
                Cfinterpo(2) = Cftable(interpoindex(2))
                Cfinterpo(3) = Cftable(interpoindex(3))
                Cfinterpo(5) = Cftable(interpoindex(5))
                Cfinterpo(6) = Cftable(interpoindex(6))
                Cfinterpo(7) = Cftable(interpoindex(7))

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
                if(isWtable(interpoindex(5)) == 0) then
                    sum = sum + Cfinterpo(5)
                    num = num + 1
                end if
                if(isWtable(interpoindex(6)) == 0) then
                    sum = sum + Cfinterpo(6)
                    num = num + 1
                end if
                if(isWtable(interpoindex(7)) == 0) then
                    sum = sum + Cfinterpo(7)
                    num = num + 1
                end if
                local_Cf = sum/num
            end if
        elseif (.not.isN) then
            if (interpoindex(4) /= -1) then
                Cfinterpo(1) = Cftable(interpoindex(1))
                Cfinterpo(2) = Cftable(interpoindex(2))
                Cfinterpo(3) = Cftable(interpoindex(3))
                Cfinterpo(4) = Cftable(interpoindex(4))
                Cfinterpo(5) = Cftable(interpoindex(5))
                Cfinterpo(6) = Cftable(interpoindex(6))
                Cfinterpo(7) = Cftable(interpoindex(7))
                Cfinterpo(8) = Cftable(interpoindex(8))

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
                if(isNtable(interpoindex(5)) == 0) then
                    sum = sum + Cfinterpo(5)
                    num = num + 1
                end if
                if(isNtable(interpoindex(6)) == 0) then
                    sum = sum + Cfinterpo(6)
                    num = num + 1
                end if
                if(isNtable(interpoindex(7)) == 0) then
                    sum = sum + Cfinterpo(7)
                    num = num + 1
                end if
                if(isNtable(interpoindex(8)) == 0) then
                    sum = sum + Cfinterpo(8)
                    num = num + 1
                end if
                local_Cf = sum/num
            else
                Cfinterpo(1) = Cftable(interpoindex(1))
                Cfinterpo(2) = Cftable(interpoindex(2))
                Cfinterpo(3) = Cftable(interpoindex(3))
                Cfinterpo(5) = Cftable(interpoindex(5))
                Cfinterpo(6) = Cftable(interpoindex(6))
                Cfinterpo(7) = Cftable(interpoindex(7))

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
                if(isNtable(interpoindex(5)) == 0) then
                    sum = sum + Cfinterpo(5)
                    num = num + 1
                end if
                if(isNtable(interpoindex(6)) == 0) then
                    sum = sum + Cfinterpo(6)
                    num = num + 1
                end if
                if(isNtable(interpoindex(7)) == 0) then
                    sum = sum + Cfinterpo(7)
                    num = num + 1
                end if
                local_Cf = sum/num
            end if
        else
            if (interpoindex(4) /= -1) then
                Cfinterpo(1) = Cftable(interpoindex(1))
                Cfinterpo(2) = Cftable(interpoindex(2))
                Cfinterpo(3) = Cftable(interpoindex(3))
                Cfinterpo(4) = Cftable(interpoindex(4))
                Cfinterpo(5) = Cftable(interpoindex(5))
                Cfinterpo(6) = Cftable(interpoindex(6))
                Cfinterpo(7) = Cftable(interpoindex(7))
                Cfinterpo(8) = Cftable(interpoindex(8))

                buf1 = (1-dp)*Cfinterpo(1) + dp*Cfinterpo(5)
                buf2 = (1-dp)*Cfinterpo(3) + dp*Cfinterpo(7)
                buf3 = (1-dp)*Cfinterpo(2) + dp*Cfinterpo(6)
                buf4 = (1-dp)*Cfinterpo(4) + dp*Cfinterpo(8)

                buf5 = (1-dz1)*buf1 + dz1*buf2
                buf6 = (1-dz1)*buf3 + dz1*buf4

                local_Cf = (1-dz2)*buf5 + dz2*buf6
            else
                Cfinterpo(1) = Cftable(interpoindex(1))
                Cfinterpo(2) = Cftable(interpoindex(2))
                Cfinterpo(3) = Cftable(interpoindex(3))
                Cfinterpo(5) = Cftable(interpoindex(5))
                Cfinterpo(6) = Cftable(interpoindex(6))
                Cfinterpo(7) = Cftable(interpoindex(7))

                sum = Cfinterpo(1)+Cfinterpo(2)+Cfinterpo(3)+Cfinterpo(5)+Cfinterpo(6)+Cfinterpo(7)
                local_Cf = sum/6
            end if
        end if

    end subroutine flashcalculation_fullgrid_3c

end module RST_flashcalculation_fullgrid_3c
