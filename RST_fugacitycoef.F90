
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_fugacitycoef

    use RST_globalFlashData
    use RST_PREOS

    implicit none

contains

    subroutine fugacitycoef(x, y, P, ZL, ZG, am, bm, al, ag, bl, bg, bigAL, bigAG, bigBL, &!
        bigBG, xiL, xiG, rhoL, rhoG, CfL, CfG, phil, phig)

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
        real(kind=8), dimension(:), pointer, intent(in out) :: phil, phig

        real(kind=8) :: CL, CG
        integer :: i, j

        call PREOS( x, y, P, ZL, ZG, am, bm, al, ag, bl, bg, bigAL, bigAG, bigBL, &!
            bigBG, xiL, xiG, rhoL, rhoG, CfL, CfG )

        do i = 1, Nc
            CL = 0.D0
            do j = 1, Nc
                CL = CL + x(j)*(1.D0-delta(i,j))*dsqrt(am(i)*am(j)*1.D0)
            end do
            CL = CL * 2.D0/al

            CG = 0.D0
            do j = 1, Nc
                CG = CG + y(j)*(1.D0-delta(i,j))*dsqrt(am(i)*am(j)*1.D0)
            end do
            CG = CG * 2.D0/ag

            phil(i) = dexp(bm(i)/bl*(ZL-1.D0)-dlog((ZL-bigBL)*1.D0) - bigAL/(2.D0*sqrt(2.D0)*&!
                bigBL)*(CL-bm(i)/bl)*dlog((ZL+2.414D0*bigBL)/(ZL-4.14D-1*bigBL)*1.D0))
            phig(i) = dexp(bm(i)/bg*(ZG-1.D0)-dlog((ZG-bigBG)*1.D0) - bigAG/(2.D0*sqrt(2.D0)*&!
                bigBG)*(CG-bm(i)/bg)*dlog((ZG+2.414D0*bigBG)/(ZG-4.14D-1*bigBG)*1.D0))

        end do

    end subroutine fugacitycoef

end module RST_fugacitycoef
