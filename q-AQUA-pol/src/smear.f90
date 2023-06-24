subroutine smear01(drsq, pol12, a, ts0, ts1)
implicit none
double precision, parameter :: g23=1.3541179394264d0
double precision, intent(in) :: drsq, pol12, a
double precision, intent(out) :: ts0, ts1
double precision :: dd,dri,drsqi,AA, rA,rA3,exp1, a_sqrt3 
double precision, external :: gammq

dd = dsqrt(drsq)
dri = 1.d0/dd
drsqi = dri*dri

AA = (pol12)**(1.d0/6.d0)
rA = dd/AA
rA3 = rA**3
exp1 = dexp(-a*rA3)
a_sqrt3 = a**(1.d0/3.d0)

ts0 = (1.d0-exp1 + a_sqrt3*rA*g23*gammq(2.d0/3.d0,a*ra3))*dri
ts1 = (1.d0 -exp1)*dri*drsqi
end subroutine smear01
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine smear1(drsq, pol12, a, ts1)
implicit none
double precision, intent(in) :: drsq, pol12, a
double precision, intent(out) :: ts1
double precision :: dd,dri,drsqi,AA, rA,rA3,exp1
double precision, external :: gammq

dd = dsqrt(drsq)
dri = 1.d0/dd
drsqi = dri*dri

AA = (pol12)**(1.d0/6.d0)
rA = dd/AA
rA3 = rA**3
exp1 = dexp(-a*rA3)
ts1 = (1.d0 -exp1)*dri*drsqi
end subroutine smear1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine smear2(drsq, pol12, a, ts1, ts2)
implicit none
double precision, intent(in) :: drsq, pol12, a
double precision, intent(out) :: ts1, ts2
double precision :: dd,dri,drsqi,AA, rA,rA3,exp1

dd = dsqrt(drsq)
dri = 1.d0/dd
drsqi = dri*dri

AA = (pol12)**(1.d0/6.d0)
rA = dd/AA
rA3 = rA**3
exp1 = dexp(-a*rA3)
ts1 = (1.d0 -exp1)*dri*drsqi
ts2 = (ts1 - exp1*a/AA**3) *drsqi
end subroutine smear2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine smear3(drsq, pol12, a, ts1, ts2, ts3)
implicit none
double precision, intent(in) :: drsq, pol12, a
double precision, intent(out) :: ts1, ts2, ts3
double precision :: dd,dri,drsqi,AA, rA,rA3,exp1

dd = dsqrt(drsq)
dri = 1.d0/dd
drsqi = dri*dri

AA = (pol12)**(1.d0/6.d0)
rA = dd/AA
rA3 = rA**3
exp1 = dexp(-a*rA3)
ts1 = (1.d0 -exp1)*dri*drsqi
ts2 = (ts1 - exp1*a/AA**3) *drsqi
ts3 = (ts2 - 0.6d0*exp1*dd*a*a/AA**6) * drsqi

end subroutine smear3

