module model_mod
   implicit none
    integer :: imodel   ! 2:ttm2f, 21: ttm21f,  3:ttm3f
    double precision :: vdwA, vdwB, vdwC, vdwD, vdwE
    double precision :: aDD, aCCaCD
    double precision :: polarO, polarH, polarM
    double precision :: polfacO, polfacH, polfacM
    double precision :: gammaM
    double precision :: dms_param1, dms_param2, dms_param3

integer, parameter          :: MAXITER  = 400
double precision, parameter :: diptol   = 1.d-20
double precision, parameter :: dmix     = 0.7d0
logical :: debug = .false.
!----------------------------------------------------------------------------!
! CONSTANTS                                                                  !
!----------------------------------------------------------------------------!
!double precision, parameter :: CHARGECON = 18.22261720426243437986d0
double precision, parameter :: CHARGECON = 18.22234397655801030455d0
double precision, parameter :: DEBYE  = 4.8033324d0

end module model_mod
