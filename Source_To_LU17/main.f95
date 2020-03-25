! A fortran95 program for G95
! By WQY
program main

  implicit none
  integer:: JMAX,IMAX,KBMAX,KMAX
  integer:: IOSTAT
  logical:: OPENED,EXIST,NAMED
  character(len=10):: ACCESS,SEQUENTIAL,DIRECT
  real,allocatable :: C(:,:,:,:),CS(:,:,:,:)
  character(len=32) :: InputSource
  
  call getarg(1,InputSource)
  read(InputSource,'(i)') KMAX

  call getarg(2,InputSource)
  read(InputSource,'(i)') IMAX
  
  call getarg(3,InputSource)
  read(InputSource,'(i)') JMAX
  
  call getarg(4,InputSource) 
  read(InputSource,'(i)') KBMAX
  
  write(*,*) InputSource
  write(*,*) KMAX,KBMAX,IMAX,JMAX
  
  allocate(C(JMAX,IMAX,KBMAX,KMAX))
  allocate(CS(JMAX,IMAX,KBMAX,KMAX))
  inquire(17,IOSTAT=IOSTAT,OPENED=OPENED,EXIST=EXIST,NAMED=NAMED,ACCESS=ACCESS,SEQUENTIAL=SEQUENTIAL,DIRECT=DIRECT)
  
  call ReadSource(C,JMAX,IMAX,KBMAX,KMAX)
  call ToLU17(C,JMAX,IMAX,KBMAX,KMAX)
  call FromLU17(CS,JMAX,IMAX,KBMAX,KMAX)
  print*,KMAX
  print*, CS
  print*, IOSTAT,OPENED,EXIST,NAMED,ACCESS,SEQUENTIAL,DIRECT

  contains
  subroutine ReadSource(C,JMAX,IMAX,KBMAX,KMAX)
!   programa para lectura de los precursores de archivo externo
    implicit none
    integer,intent(in):: JMAX,IMAX,KBMAX,KMAX
    real,intent(inout):: C(JMAX,IMAX,KBMAX,KMAX)
    character(len=1024) :: format_len,reading_format = "E15.7)"
    integer :: J,I,KB,k
    write(format_len,'(A1,I5,A7)') '(',JMAX*IMAX*KBMAX,reading_format
    open(unit=10,file='source.dat',status='old',access='sequential',form='formatted',action='read')
    do k=1,KMAX
        read(10,fmt=format_len) (((C(J,I,KB,k),J=1,JMAX),I=1,IMAX),KB=1,KBMAX)
    end do
    close(10)
    end subroutine

  subroutine ToLU17(C,JMAX,IMAX,KBMAX,KMAX)
!   una vez leidos, pasan al archivo fort.17
    implicit none
    integer,intent(in):: JMAX,IMAX,KBMAX,KMAX
    real,intent(in):: C(JMAX,IMAX,KBMAX,KMAX)

    integer :: J,I,KB,k
    open(unit=17,form='unformatted',action='write')
    do k=1,KMAX
        write(17) (((C(J,I,KB,k),J=1,JMAX),I=1,IMAX),KB=1,KBMAX)
    end do
	endfile 17
    close(17)
    end subroutine

  subroutine FromLU17(C,JMAX,IMAX,KBMAX,KMAX)
!   solo para chequeo de lectura de archivo fort.17
    implicit none
    integer,intent(in):: JMAX,IMAX,KBMAX,KMAX
!    real,intent(inout):: C(JMAX,IMAX,KBMAX,KMAX)
    real,intent(inout):: C(JMAX,IMAX,KBMAX)
    integer :: J,I,KB,k
    do k=1,KMAX
!        read(17) (((C(J,I,KB,k),J=1,JMAX),I=1,IMAX),KB=1,KBMAX)
        read(17) C
    end do
    end subroutine
end

