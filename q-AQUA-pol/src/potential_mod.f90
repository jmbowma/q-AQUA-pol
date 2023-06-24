module potential_mod
use model_mod

Contains
   subroutine potential(nw, rr, drr, en)
   use ttm2f_mod
   use ttm3f_mod
   implicit none
   integer :: nw
   double precision, dimension(3, 3*nw) :: rr
   double precision, dimension(3, 3*nw) :: drr
   double precision :: en

   if (imodel==2 .or. imodel==21) then
      call ttm2f(Nw,RR,dRR,En)
   else if (imodel==3) then
      call ttm3f(Nw,RR,dRR,En)
   endif

   end subroutine potential
end module potential_mod
