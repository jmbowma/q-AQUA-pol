module nma_proc
  use pes_shell
  implicit none

  ! Define globe variables
  real,dimension(:),allocatable::nma_x     ! nma geom in bohr
  real,dimension(:),allocatable::nma_mass  ! sqrt root of atom mass
  character(len=2),dimension(:),allocatable :: nma_symb
  
contains
  !==========================================
  ! calculate the mass weighted hessian at p        
  !==========================================
  subroutine mw_hessian(p,H)
    real,dimension(:),intent(in)::p
    real,dimension(:,:),intent(inout)::H
    !:::::::::::::::::::
    real::rmass
    integer::dim,i,j

    dim=size(p)
    call hessian(p,H)

    do i=1,dim
       do j=1,dim
          rmass=sqrt(nma_mass(ceiling(i/3.0)))*sqrt(nma_mass(ceiling(j/3.0)))
          H(i,j)=H(i,j)/rmass
       end do
    end do

    return 
  end subroutine mw_hessian

  !==================================================
  ! diagonalize the hessian matrix and return        
  ! the eigen value and eigenvectors.         
  ! The original Hessian matrix  will be destroied
  !==================================================
  subroutine diag_hessian(H,w)
    real,dimension(:,:),intent(inout)::H
    real,dimension(:),intent(out)::w
    ! ::::::::::::::::::::
    real,dimension(:),allocatable::work
    integer::dim,lwork,info,i,j
    
    dim=size(H,1)
    lwork=dim*dim*10;
    allocate(work(1:lwork))
    
    call dsyev('v','u',dim,H,dim,w,work,lwork,info) 
    
    do i=1,dim
       w(i)=sign(sqrt(abs(w(i)))*aucm,w(i))
    end do
    
    return
  end subroutine diag_hessian

  !==================================================
  ! print vectors that can be visualized by xmakemole
  !==================================================
  subroutine prtxvec(x,q,symbs,mode,f)
    real,dimension(:),intent(in)::x
    real,dimension(:),intent(in)::q
    character(len=2),dimension(:),intent(in)::symbs
    integer,intent(in)::mode
    integer,intent(in)::f
    
    integer::i,natm,dim

    natm=size(symbs,1)

    write(f,'(I2)') natm
    write(f,'(A,I8)') "Mode",mode 

    do i=1,natm
       write(f,'(1x,A,3F13.8,2x,A,3F10.5)') symbs(i),x(3*i-2),x(3*i-1),x(3*i),&
            'atom_vector',q(3*i-2),q(3*i-1),q(3*i)
    end do

    return
  end subroutine prtxvec
  
end module nma_proc
