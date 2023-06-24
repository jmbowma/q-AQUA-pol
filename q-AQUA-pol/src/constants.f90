module constants
  implicit none

  ! Define the mass of different atoms
  real,parameter::c_mass= 12.0000000  !21874.66
  real,parameter::d_mass=  2.0135532127
  real,parameter::h_mass=  1.0078250  !1837.15
  real,parameter::o_mass= 15.9949146  !29156.95
  real,parameter::ar_mass=39.9623831  !72846.97
  real,parameter::cl_mass=34.968853388
  real,parameter::pi=acos(-1.0)

  ! Define constants
  real,parameter::raddeg=57.2957795
  real,parameter::auang=0.5291772083
  real,parameter::aucm=219474.6313710
  real,parameter::cmau=0.00000456
  real,parameter::aukcal=627.51
  real,parameter::cmhz=29979250000.0
  real,parameter::auhz=41341370000000000.0

  real,parameter::emass=1822.88848

end module constants
