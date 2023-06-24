module ttm3f_mod
use model_mod
!----------------------------------------------------------------------------!
! Parameters of the TTM3F potential                                         !
!----------------------------------------------------------------------------!
! (12-10-6) inverse polynomial + exponential parameters for the vdw ineractions
!double precision, parameter ::   &
!                                 vdwC=-0.72298855D+03,vdwD=0.10211829D+06, vdwE=0.37170376D+01
!! smearing factors for dipole-dipole(aDD), charge-charge/charge-dipole(aCCaCD)
!double precision, parameter :: aDD=0.175d0, aCCaCD=0.175d0
!!   ........... polarizabilities
!double precision, parameter :: polarM=1.444d0
!double precision, parameter :: polfacO=0.837d0, polfacH = 0.496d0, polfacM=0.837d0
!double precision, parameter :: dms_param1 =0.5d0,dms_param2 = 0.9578d0, dms_param3=0.012d0
!!   ...........   gammaM
!double precision, parameter :: gammaM=0.46d0
!!----------------------------------------------------------------------------!
!! Parameters related to the accuracy of the energy/derivatives calculation   !
!!----------------------------------------------------------------------------!
!integer, parameter          :: MAXITER  = 400
!double precision, parameter :: diptol   = 1.d-15 
!double precision, parameter :: dmix     = 0.7d0
!!----------------------------------------------------------------------------!
!! CONSTANTS                                                                  !
!!----------------------------------------------------------------------------!
!double precision, parameter :: CHARGECON = 18.22234397655801030455d0
!double precision, parameter :: DEBYE  = 4.8033324d0
!----------------------------------------------------------------------------!
! Variables and allocatable arrays needed by the "ttm3f" subroutine         !
!----------------------------------------------------------------------------!
integer :: fO, lO, fH, lH, fM, lM, fO3, lO3, fH3, lH3, fM3, lM3
integer :: Nw_old
integer :: Nats, Natsd, Natsq
double precision, dimension(:,:), allocatable :: RM
double precision, dimension(:,:), allocatable :: dRM
double precision, dimension(:,:), allocatable :: DDT
double precision, dimension( : ), allocatable :: dip
double precision, dimension( : ), allocatable :: pr_dip
double precision, dimension( : ), allocatable :: Efq
double precision, dimension( : ), allocatable :: Efd
double precision, dimension( : ), allocatable :: charge
double precision, dimension( : ), allocatable :: phi
double precision, dimension(:,:,:,:), allocatable :: grdq

Contains

   !***
   !***  allocate necessery arrays needed for the calculation of the potenial
   !***
   subroutine init_ttm3f(Nw)
!   use ttm3f_mod
   implicit none
   integer, intent(in) :: Nw
   logical :: alloc

   vdwC=-0.72298855D+03
   vdwD=0.10211829D+06
   vdwE=0.37170376D+01
   aDD=0.175d0 
   aCCaCD=0.175d0
   polarM=1.444d0
   polfacO=0.837d0
   polfacH = 0.496d0
   polfacM=0.837d0
   dms_param1 =0.5d0
   dms_param2 = 0.9578d0
   dms_param3=0.012d0
   gammaM=0.46d0

   
   
   if (.not. allocated(RM)) then
      alloc = .true.
   else if (Nw/=Nw_old) then
      alloc = .true.
      deallocate( RM )    ! temporary array keeping the coordinates of the Msite
      deallocate( dRM )   ! temporary array keeping the derivatives of the Msite
      deallocate( DDT )   ! dipole tensor
      deallocate( dip )   ! array containg the induced dipoles
      deallocate( pr_dip )  ! 
      deallocate( phi )    ! electrostatic potential
      deallocate( Efq )     ! electric field from charges
      deallocate( Efd )     ! electric field from dipoles
      deallocate( charge )  ! charges 
      deallocate( grdq )    ! dericatives of charge wrt monomer geometry
   else
      alloc = .false.
   endif
   
   if (alloc) then
      Natsq = 3*Nw       ! # of atoms with charge   (including M-sites)
      Natsd =   Nw       ! # of atoms with dipole   (including M-sites)
      Nats  = 4*Nw       ! # of total atoms         (including M-sites)
      Nw_old = Nw
      fO = 1             ! index on the first oxygen
      lO = Nw            ! index on the last oxygen        
      fH = Nw+1          ! index on the first hydrogen
      lH = 3*Nw          ! index on the last hydrogen
      fM = 3*Nw+1        ! index on the first M-site
      lM = 4*Nw          ! index on the lst M-site
      fO3 = 3*fO-2  ! index on the x-comp. (out of the x,y,z) of the first oxygen
      lO3 = 3*lO    ! index on the x-comp. (out of the x,y,z) of the last oxygen
      fH3 = 3*fH-2  ! index on the x-comp. (out of the x,y,z) of the first hydrogen
      lH3 = 3*lH    ! index on the x-comp. (out of the x,y,z) of the last hydrogen
      fM3 = 3*fM-2  ! index on the x-comp. (out of the x,y,z) of the first M-site
      lM3 = 3*lM    ! index on the x-comp. (out of the x,y,z) of the last M-site
   
      allocate(charge(fO:lM))          
      allocate(RM (3, fM:lM))          
      allocate(dRM(3, fM:lM))         
      allocate(dip(fM3:lM3))
      allocate(pr_dip(fM3:lM3))
      allocate(phi(fH:lM))
      allocate(Efq(fO3:lM3))
      allocate(Efd(fM3:lM3))
      allocate( DDT(fM3:lM3, fM3:lM3) )
      allocate( grdq(Nw, 3, 3, 3) )
   endif
   
   end subroutine init_ttm3f



   subroutine ttm3f(Nw,RR,dRR,En)
   !** for a system of "Nw" water molecules returns the potential energy (En) and 
   !** the derivatives of the energy wrt the atomic displacements (dRR).
   !** The coordinates of the atoms are found in the array "RR".
   !** NOTE. It is very important the input coordinates to be stored in the
   !** array RR in the proper way. The first column of RR has dimension 3 that
   !** corresponds to the x,y,z- cartesian coordinates, while in the second one 
   !** all the oxygens of the system should be saved before the hydrogens.
   !** It follows an example for a system with 3-water molecules (Nw=3)
   !**   RR(1:3, 1)   = x,y,z-coordinates of the Oxygen of the first molecule
   !**   RR(1:3, 2)   = x,y,z-coordinates of the Oxygen of the second molecule
   !**   RR(1:3, 3)   = x,y,z-coordinates of the Oxygen of the third molecule
   !**   RR(1:3, 4)   = x,y,z-coordinates of the hydrogen-1 of the first molecule
   !**   RR(1:3, 5)   = x,y,z-coordinates of the hydrogen-2 of the first molecule
   !**   RR(1:3, 6)   = x,y,z-coordinates of the hydrogen-1 of the second molecule
   !**   RR(1:3, 7)   = x,y,z-coordinates of the hydrogen-2 of the second molecule
   !**   RR(1:3, 8)   = x,y,z-coordinates of the hydrogen-1 of the third molecule
   !**   RR(1:3, 9)   = x,y,z-coordinates of the hydrogen-2 of the third molecule
!   use ttm3f_mod
   implicit none
   integer, intent(in) :: Nw
   double precision, dimension(3, 3*Nw), intent(in)  ::  RR
   double precision, dimension(3, 3*Nw), intent(out) :: dRR
   double precision, intent(out) :: En
   !
   integer :: iter, i3, j3, isp, jsp, ix, iy, iw, jw, iat, jat, iO, iH1, iH2, iM
   double precision, dimension(3) :: Ri, Rij
   double precision, dimension(3,3) :: r1, dr1, dd3
   double precision, dimension(3) :: q3, di, dj, derij
   double precision, dimension(3,3,3) :: dq3
   double precision :: dRijsq, dRij, dR6, dR10, dR12
   double precision :: Eint, Evdw, Eelec, Eind, e1
   double precision :: deltadip,tmp,expon, stath, didj, dir, djr,qi,qj
   double precision :: ts0, ts1, ts2, ts3
   
   !-------------------------------------------------------------------------!
   ! initialize (if neccesary) some arrays, needed for following calculations!
   !-------------------------------------------------------------------------!
   call init_ttm3f(Nw)
   dRR = 0.d0
   dRM = 0.d0
   Eint = 0.d0
   Evdw = 0.d0
   Eelec = 0.d0
   Eind = 0.d0
   Efq = 0.d0
   phi = 0.d0
   ddt = 0.d0
   !-------------------------------------------------------------------------!
   ! calculate the coordinates of the Msites and store them in RM            !
   !-------------------------------------------------------------------------!
   do iw=1, Nw
      iO  = fO + iw-1
      iH1 = fH + 2*iw - 2
      iH2 = fH + 2*iw - 1
      iM  = fM + iw-1
      RM(1:3,iM) = RR(1:3,iO)*(1.d0-gammaM) + &
                              0.5d0*gammaM*( RR(1:3,iH1) + RR(1:3,iH2) )
   enddo
   !-------------------------------------------------------------------------!
   ! calculates the INTRA-molecular energy and Dipole Moment Surface         !
   !                     according the Partridge Sw  (JCP     )              !
   !-------------------------------------------------------------------------!
   tmp = 0.5d0*gammaM/(1.d0-gammaM)
   do iw=1, Nw
      iO  = fO + iw-1
      iH1 = fH+2*iw-2
      iH2 = fH+2*iw-1
      iM  = fM + iw-1
      r1(1:3, 1:3) = RR(1:3, (/iO, iH1, iH2/) )
      call pot_nasa(r1, dr1, e1)
      call dms_nasa(r1, q3, dq3)
      Eint = Eint + e1
      dRR(1:3, (/iO,iH1,iH2/) ) = dr1
      charge(iH1)  = q3(2)+tmp*(q3(2)+q3(3) )
      charge(iH2)  = q3(3)+tmp*(q3(2)+q3(3) )
      charge(iM )  = q3(1) / (1.d0-gammaM)
      grdq(iw,1,1,:)= dq3(1,1,:) + tmp*(dq3(1,1,:)+dq3(1,2,:))
      grdq(iw,2,1,:)= dq3(2,1,:) + tmp*(dq3(2,1,:)+dq3(2,2,:))
      grdq(iw,3,1,:)= dq3(3,1,:) + tmp*(dq3(3,1,:)+dq3(3,2,:))
   
      grdq(iw,1,2,:)= dq3(1,2,:) + tmp*(dq3(1,1,:)+dq3(1,2,:))
      grdq(iw,2,2,:)= dq3(2,2,:) + tmp*(dq3(2,1,:)+dq3(2,2,:))
      grdq(iw,3,2,:)= dq3(3,2,:) + tmp*(dq3(3,1,:)+dq3(3,2,:))
   
      grdq(iw,1,3,:)= dq3(1,3,:)-2.d0*tmp*(dq3(1,1,:)+dq3(1,2,:))
      grdq(iw,2,3,:)= dq3(2,3,:)-2.d0*tmp*(dq3(2,1,:)+dq3(2,2,:))
      grdq(iw,3,3,:)= dq3(3,3,:)-2.d0*tmp*(dq3(3,1,:)+dq3(3,2,:))
   enddo
   charge = charge*CHARGECON
   grdq   = grdq*CHARGECON
   !-------------------------------------------------------------------------!
   ! Calculate the CHARGE-CHARGE interactions for all atoms                  !
   !-------------------------------------------------------------------------!
   do iw=1, Nw-1
      do jw=iw+1, Nw
         !Oxygen-Oxygen interactions
         iat=fO+iw-1          ! iat=Oxygen-1   
         jat=fO+jw-1          ! jat=Oxygen-2
         Rij=RR(:,iat) - RR(:,jat)
         dRijsq=Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
         dRij = dsqrt(dRijsq)
         !... vdw interactions
         dR6 = dRijsq**3
         dR10=dR6*dRijsq*dRijsq
         dR12=dR6*dR6
         expon = vdwD*dexp(-vdwE*dRij)
         Evdw = Evdw + vdwC/dR6 + expon
         tmp=-(6.d0*vdwC/dR6)/dRijsq -vdwE*expon/dRij
         dRR(:,iat) = dRR(:,iat) + tmp * Rij
         dRR(:,jat) = dRR(:,jat) - tmp * Rij
         !Hydrogen-Hydrogen interactions
         do isp=1, 2
            iat=fH+2*(iw-1)+ isp-1   ! iat=Hydrogen-1A/Hydrogen-1B
            Ri = RR(1:3, iat)
            do jsp=1, 2
               jat=fH+2*(jw-1)+ jsp-1   ! jat=Hydrogen-2A/Hydrogen-2B
               Rij = Ri - RR(1:3, jat)
               dRijsq=Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
               call smear01(dRijsq, polfacH**2, aCCaCD, ts0, ts1)
               phi(iat) = phi(iat) + ts0*charge(jat)
               phi(jat) = phi(jat) + ts0*charge(iat)
               Efq(3*iat-2:3*iat) = Efq(3*iat-2:3*iat) + ts1*charge(jat)*Rij
               Efq(3*jat-2:3*jat) = Efq(3*jat-2:3*jat) - ts1*charge(iat)*Rij
            enddo
         enddo
         !Msite-Msite interactions
         iat=fM+iw-1          ! iat=Msite-1   
         jat=fM+jw-1          ! jat=Msite-2
         Rij=RM(:,iat) - RM(:,jat)
         dRijsq=Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
         call smear01(dRijsq, polfacM**2, aCCaCD, ts0, ts1)
         phi(iat) = phi(iat) + ts0*charge(jat)
         phi(jat) = phi(jat) + ts0*charge(iat)
         Efq(3*iat-2:3*iat) = Efq(3*iat-2:3*iat) + ts1*charge(jat)*Rij
         Efq(3*jat-2:3*jat) = Efq(3*jat-2:3*jat) - ts1*charge(iat)*Rij
      enddo
   enddo
   do iw=1, Nw
      do jw=1, Nw
         if (iw/=jw) then
            !Msite-Hydrogen interactions
            iat = fM + iw -1   ! iat=Msite
            Ri = RM(1:3, iat)
            do jsp=1,2
               jat = fH+2*(jw-1) + jsp-1  !jat = hydrogen-2
               Rij = Ri - RR(1:3, jat)
               dRijsq=Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
               call smear01(dRijsq, polfacH*polfacM, aCCaCD, ts0, ts1)
               phi(iat) = phi(iat) + ts0*charge(jat)
               phi(jat) = phi(jat) + ts0*charge(iat)
               Efq(3*iat-2:3*iat) = Efq(3*iat-2:3*iat) + ts1*charge(jat)*Rij
               Efq(3*jat-2:3*jat) = Efq(3*jat-2:3*jat) - ts1*charge(iat)*Rij
            enddo
         endif 
      enddo
   enddo
   !-------------------------------------------------------------------------!
   ! Calculate the DIPOLE-DIPOLE TENSOR Array.  (according to the Tholes   --!  
   !   model the intra-molecular interactions  should be also considered)  --!
   !-------------------------------------------------------------------------!
   do iat=fM, lM-1
      i3=3*iat-2
      Ri = RM(:, iat)
      do jat=iat+1, lM
         j3=3*jat-2
         Rij = Ri - RM(:, jat)
         dRijsq=Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
         call smear2(dRijsq, polfacM**2, aDD, ts1, ts2)
         dd3(1,1) = 3.d0*ts2*Rij(1)*Rij(1) - ts1
         dd3(2,2) = 3.d0*ts2*Rij(2)*Rij(2) - ts1
         dd3(3,3) = 3.d0*ts2*Rij(3)*Rij(3) - ts1
         dd3(1,2) = 3.d0*ts2*Rij(1)*Rij(2)
         dd3(1,3) = 3.d0*ts2*Rij(1)*Rij(3)
         dd3(2,3) = 3.d0*ts2*Rij(2)*Rij(3)
         dd3(2,1)=dd3(1,2); dd3(3,1)=dd3(1,3); dd3(3,2)=dd3(2,3)
         ddt(i3:i3+2, j3:j3+2) = dd3
         ddt(j3:j3+2, i3:i3+2) = dd3
      enddo  ! do jat=iat+1, Natsd
   enddo  ! do iat=1, Natsd
   !-------------------------------------------------------------------------!
   ! Calculate the induced Electric Field using an iterative proced.         !
   !-------------------------------------------------------------------------!
   dip(fM3:lM3) = polarM*Efq(fM3:lM3)
   pr_dip(fM3:lM3) = dip(fM3:lM3)   ! keep the previous dipole
   
   stath = DEBYE/CHARGECON/dsqrt(dble(Natsd))
   do iter=1, MAXITER
      Efd = matmul(ddt, dip)
      dip(fM3:lM3) = polarM*( Efq(fM3:lM3) + Efd(fM3:lM3) )
      dip(fM3:lM3) = dmix*dip(fM3:lM3) + (1.d0-dmix)*pr_dip(fM3:lM3)
      deltadip = sum(  (dip(fM3:lM3)-pr_dip(fM3:lM3))**2  )
      deltadip = dsqrt(deltadip)*stath
   !   print*,'iter=',iter, deltadip
      if (deltadip<diptol) then
         goto 100
      else
         pr_dip(fM3:lM3) = dip(fM3:lM3)
      endif
   enddo
   100 Continue
   Eelec = 0.5d0*sum(charge(fH:lM)*phi(fH:lM))

   Eind = -0.5d0*sum(dip(fM3:lM3)*Efq(fM3:lM3))
   En = Eint + Evdw + Eelec + Eind
   debug = .false.
   if (debug) then
      print*,'Eint=',Eint
      print*,'Evdw', Evdw
      print*,'Eelec=',Eelec
      print*,'Eind=',Eind
   endif

!   print*,'Eint=',Eint; print*,'Evdw', Evdw; print*,'Eelec=',Eelec; print*,'Eind=',Eind
   !-------------------------------------------------------------------------!
   !---------  Calculate the remaining part of the derivatives   ------------!
   !-------------------------------------------------------------------------!
   !....... derivatives due to charge-charge interaction
   do iat=fH, lH
      dRR(:,iat) = dRR(:,iat) - charge(iat)*Efq(3*iat-2:3*iat)
   enddo
   do iat=fM, lM
      dRM(:,iat) = dRM(:,iat) - charge(iat)*Efq(3*iat-2:3*iat)
   enddo
   !....... derivatives due to charge-dipole and dipole-dipole interaction
   do iw=1, Nw-1
      do jw=iw+1, Nw
         !Msite-Msite interactions
         iat=fM+iw-1          ! iat=Msite-1   
         jat=fM+jw-1          ! jat=Msite-2
         Rij=RM(:,iat) - RM(:,jat)
         dRijsq=Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
         call smear3(dRijsq, polfacM**2, aDD, ts1, ts2, ts3)
         qi = charge(iat)
         qj = charge(jat)
         di = dip(3*iat-2:3*iat)  ! dipole-I
         dj = dip(3*jat-2:3*jat)  ! dipole-J
         didj = di(1)*dj(1) + di(2)*dj(2) + di(3)*dj(3)
         dir = di(1)*Rij(1) + di(2)*Rij(2) + di(3)*Rij(3)
         djr = dj(1)*Rij(1) + dj(2)*Rij(2) + dj(3)*Rij(3)
         derij=-3.d0*ts2*(didj*Rij+djr*di+dir*dj) + 15.d0*ts3*dir*djr*Rij 
         derij=derij-3.d0*ts2*qi*djr*Rij + ts1*qi*dj
         derij=derij+3.d0*ts2*qj*dir*Rij - ts1*qj*di
         dRM(:,iat) = dRM(:,iat) + derij 
         dRM(:,jat) = dRM(:,jat) - derij 
         phi(iat) = phi(iat) + ts1*djr
         phi(jat) = phi(jat) - ts1*dir
      enddo
   enddo
   !Hydrogen-Msite interactions
   do iw=1, Nw
      iat = fM + iw -1   ! iat=Msite
      qi = charge(iat)         ! charge-I
      di = dip(3*iat-2:3*iat)
      Ri = RM(1:3, iat)
      do jw=1, Nw
         if (iw/=jw) then
            do jsp=1,2
               jat = fH+2*(jw-1) + jsp-1  !jat = hydrogen-2
               qj = charge(jat)         ! charge-J
               Rij = Ri - RR(1:3, jat)
               dRijsq=Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
               call smear2(dRijsq, polfacH*polfacM, aCCaCD, ts1, ts2)
               dir  = di(1)*Rij(1) + di(2)*Rij(2) + di(3)*Rij(3)
               derij=   3.d0*ts2*qj*dir*Rij - ts1*qj*di
               dRM(:,iat) = dRM(:,iat) + derij
               dRR(:,jat) = dRR(:,jat) - derij
               phi(jat) = phi(jat)-ts1*dir
            enddo
         endif
      enddo
   enddo
   !----derivatives from the adjustable charges of the NASA PES
   do iw=1,Nw
      iO  = fO + iw-1
      iH1 = fH + 2*iw-2
      iH2 = fH + 2*iw-1
      iM  = fM + iw-1
      dRR(:,iH1)=dRR(:,iH1)+(grdq(iw,1,1,:)*phi(iH1)+  &
                             grdq(iw,1,2,:)*phi(iH2)+grdq(iw,1,3,:)*phi(iM))
      dRR(:,iH2)=dRR(:,iH2)+(grdq(iw,2,1,:)*phi(iH1)+  & 
                             grdq(iw,2,2,:)*phi(iH2)+grdq(iw,2,3,:)*phi(iM))
      dRR(:,iO )=dRR(:,iO )+(grdq(iw,3,1,:)*phi(iH1)+  &  
                             grdq(iw,3,2,:)*phi(iH2)+grdq(iw,3,3,:)*phi(iM))
   enddo
   !-------------------------------------------------------------------------!
   !-- Redistribute the ders from the sites (H,H,M) to the atoms (O,H,H)     !
   !-------------------------------------------------------------------------!
   do iw=1, Nw
      iO  = fO + iw-1
      iH1 = fH + 2*iw-2
      iH2 = fH + 2*iw-1
      iM  = fM + iw-1
      dRR(:,iH1) = dRR(:,iH1) +  0.5d0*gammaM*dRM(:,iM)
      dRR(:,iH2) = dRR(:,iH2) +  0.5d0*gammaM*dRM(:,iM)
      dRR(:, iO) = dRR(:, iO) + (1.d0-gammaM)*dRM(:,iM)
   enddo
   
   end subroutine ttm3f
end module ttm3f_mod
