
!!$ Author:
!!$   Yuanqing Wu, KAUST, Saudi Arabia
!!$
!!$ History:
!!$   2015-2-9 by Yuanqing Wu
!!$
!!$ Support:
!!$   wuyuanq@gmail.com

module RST_globalFlashData

    implicit none

    ! the program parameters
    integer, parameter :: MAXTIME = 1000
#ifdef FULL_2c
    integer, parameter :: TABLESIZE = 129*129
    character(len=*), parameter :: FULLGRIDPREFIX = '/Users/yuanqingwu/research/CompositionalFlow/Fortran/Fullgrid_2c/'
#elif FULL_3c
    integer, parameter :: TABLESIZE = 17*(1+257)/2*257
    character(len=*), parameter :: FULLGRIDPREFIX = '/Users/yuanqingwu/research/CompositionalFlow/Fortran/Fullgrid_3c/'
#endif

    ! the physical parameters
    real(kind=8), parameter :: R = 8.314

    ! the physical variables
    real(kind=8) :: Temp
    real(kind=8), dimension(:), pointer :: ct 
    real(kind=8), dimension(:), pointer :: cp 
    real(kind=8), dimension(:), pointer :: af 
    real(kind=8), dimension(:), pointer :: mw 
    real(kind=8), dimension(:), pointer :: cv 
    real(kind=8), dimension(:), pointer :: psatA 
    real(kind=8), dimension(:), pointer :: psatB
    real(kind=8), dimension(:), pointer :: psatC
    real(kind=8), dimension(:,:), pointer :: delta

    ! the global variables
    integer :: Nc
#if defined(FULL_2c) || defined(FULL_3c)
    real(kind=8), dimension(:,:), pointer :: xtable
    real(kind=8), dimension(:,:), pointer :: ytable
    real(kind=8), dimension(:), pointer :: xiLtable
    real(kind=8), dimension(:), pointer :: xiGtable
    real(kind=8), dimension(:), pointer :: rhoLtable
    real(kind=8), dimension(:), pointer :: rhoGtable
    real(kind=8), dimension(:), pointer :: sLtable
    real(kind=8), dimension(:,:), pointer :: vtable
    real(kind=8), dimension(:), pointer :: Cftable
    integer, dimension(:), pointer :: isWtable
    integer, dimension(:), pointer :: isNtable
#elif SPARSE
    logical :: isFirstSG
#elif NN
    logical :: isFirstNN
#endif

end module RST_globalFlashData

