subroutine vibpot(rij,v,n)
  use pot_monomer_mod
!!$c
!!$c     pes for h2o,
!!$c     Harry Partridge and David W. Schwenke, J. Chem. Phys.,
!!$c     submitted Nov. 8, 1996.
!!$c     rij(i,1)& rij(i,2) are oh distances in au
!!$c     rij(i,3) is hoh angle in rad
!!$c     v(i) is pes in au
!!$c     n is number of geometries
!!$c     mass dependent factors are included. the nuclear masses
!!$c     should be passed to this program using the array xm in
!!$c     common potmcm. xm(1) is the
!!$c     mass of the hydrogen associated with rij(i,1), and xm(2)
!!$c     is the mass of the hydrogen associated with rij(i,2).
!!$c     all masses are in au.
  integer,intent(in)::n
  real,dimension(n,3),intent(in)::rij
  real,dimension(n),intent(inout)::v
  !::::::::::::::::::::
  integer::i,j
  real::x1,x2,x3,rhh,voh1,voh2,ex
  real::fmat(15,3),v1,v2,term,vhh
  
  v=0.d0
  do i=1,n
     x1=(rij(i,1)-reoh)/reoh
     x2=(rij(i,2)-reoh)/reoh
     x3=cos(rij(i,3))-ce
     rhh=sqrt(rij(i,1)**2+rij(i,2)**2-2.d0*rij(i,1)*rij(i,2)*cos(rij(i,3)))
     vhh=phh1*exp(-phh2*rhh)
     ex=exp(-alphaoh*(rij(i,1)-roh))
     voh1=deoh*ex*(ex-2.d0)
     ex=exp(-alphaoh*(rij(i,2)-roh))
     voh2=deoh*ex*(ex-2.d0)

     fmat(1,1)=1d0
     fmat(1,2)=1d0
     fmat(1,3)=1d0

     do j=2,15
        fmat(j,1)=fmat(j-1,1)*x1
        fmat(j,2)=fmat(j-1,2)*x2
        fmat(j,3)=fmat(j-1,3)*x3
     end do

     do j=2,245
        term=c5z(j)*(fmat(idx(j,1),1)*fmat(idx(j,2),2) &
             +fmat(idx(j,2),1)*fmat(idx(j,1),2))*fmat(idx(j,3),3)
        v(i)=v(i)+term
     end do

     v1=0d0
     v2=0d0
     do j=1,9
        v1=v1+cmass(j)*fmat(idxm(j,1),1)*fmat(idxm(j,2),2)*fmat(idxm(j,3),3)
        v2=v2+cmass(j)*fmat(idxm(j,2),1)*fmat(idxm(j,1),2)*fmat(idxm(j,3),3)
     end do
     v(i)=v(i)+xm1*v1+xm2*v2
     v(i)=v(i)*exp(-b1*((rij(i,1)-reoh)**2+(rij(i,2)-reoh)**2))+c5z(1)+voh1+voh2+vhh
  end do
  return
end subroutine vibpot
