program main
  use pes_shell

  implicit none
  real,dimension(:),allocatable::x,tmpx
  real,dimension(:,:),allocatable::gd1,gd2,H
  character(len=32)::filename
  integer::i,j,natm,ierr
  character::symb
  real::eng,time_start,time_end,eng_ccsd,engf,engb
  real::eps
 
  call getarg(1,filename)
  open(21,file=trim(filename),status="old")

  read(21,*) natm
  allocate(x(3*natm), gd1(3,natm), gd2(3,natm),H(3*natm,3*natm),tmpx(3*natm))
  rewind(21)

  call pes_init(natm/3)
  eps = 0.001

  do
     read(21,*,iostat=ierr) natm
     if (ierr < 0) exit
     read(21,*) 
     do i=1,natm
        read(21,*) symb,x(3*i-2:3*i)
     end do
     x = x / auang
    
     call pot_gd(x,eng,gd1)
     write(*,*) eng*627.51
     do i=1,natm
      write(*,'(I4,3F15.8)') i,gd1(:,i)
     end do
  end do

end program main
