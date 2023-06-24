!*****************************************************************************
!*** The following subroutines used for the calculation of the Gamma function
!*** have been taken from the "Numerical Recipes"   
!*****************************************************************************
!--------------------------------------------------
function gammln(XX)

!***   taken from the "Numerical Recipes"

!--------------------------------------------------
implicit none
double precision :: XX
double precision :: GAMMLN
integer :: J
double precision :: SER,STP,TMP,X,Y
double precision, DIMENSION(6) :: COF
SAVE cof,stp
data COF,STP/76.18009172947146D0,-86.50532032941677d0,24.01409824083091D0,  &
            -1.231739572450155D0,.1208650973866179D-2,-.5395239384953D-5, &
            2.5066282746310005D0/
!
!-------- Executable code
!
X=XX
Y=X
TMP=X+5.5D0
TMP=(X+0.5D0)*DLOG(TMP)-TMP
SER=1.000000000190015D0
do j=1,6
   Y=Y+1.D0
   SER=SER+COF(J)/Y
enddo
GAMMLN=TMP+DLOG(STP*SER/X)
RETURN
END FUNCTION gammln

!--------------------------------------------------
function gammq(a,x)

!***   taken from the "Numerical Recipes"

!--------------------------------------------------
implicit none
double precision :: gammq
double precision :: a, x
double precision :: gamser, gammcf, gln

if (x<0.d0 .or. a <= 0.d0) stop 'bad arguments in gammq'
if (x<a+1.d0) then
   call gser(gamser, a, x, gln)
   gammq = 1.d0 - gamser
else
   call gcf(gammcf, a, x, gln)
   gammq = gammcf
endif
return
end function gammq

!--------------------------------------------------
subroutine gser(gamser, a, x, gln)

!***   taken from the "Numerical Recipes"

!--------------------------------------------------
implicit none
double precision :: gamser
double precision :: a, x, gln
integer          :: n
integer, parameter :: ITMAX = 800
double precision, parameter :: EPS = 3.d-7
double precision :: ap, ssum, del, gammln

gln = gammln(a)
if (x <= 0.d0) then
   if (x < 0.d0) stop 'x<0 in gser'
   gamser = 0.d0
   return
endif

ap = a
ssum = 1.d0 / a
del =ssum

do n=1, ITMAX
   ap = ap+1.d0
   del = del*x/ap
   ssum = ssum+del
   if (dabs(del) .lt. dabs(ssum)*EPS) goto 1
enddo
stop 'a too large, ITMAX too small in gser'
1 gamser = ssum*dexp(-x+a*dlog(x)-gln)
return
end subroutine gser

!--------------------------------------------------
subroutine gcf(gammcf, a, x, gln)

!***   taken from the "Numerical Recipes"

!--------------------------------------------------
implicit none
double precision :: gammcf
double precision :: a, x, gln
integer          :: i
double precision :: an, b, c, d, del, h, gammln
integer,parameter :: ITMAX = 100
double precision, parameter :: EPS = 3.d-7, FPMIN=1.d-30

gln = gammln(a)
b = x+1.d0-a
c = 1.d0/FPMIN
d = 1.d0/b
h=d
do i=1, ITMAX
   an = -i*(i-a)
   b = b+2.d0
   d = an*d+b
   if (dabs(d) < FPMIN) d=FPMIN
   c = b+an/c
   if (dabs(c) < FPMIN) c=FPMIN
   d = 1.d0/d
   del = d*c
   h = h*del
   if (dabs(del-1.d0) < EPS) goto 1
enddo
stop ' a too large ITMAX too small in gcf'
1 gammcf=dexp(-x+a*dlog(x)-gln)*h
return
end subroutine gcf
