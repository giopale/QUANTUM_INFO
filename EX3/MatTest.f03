! **********************************************************************
! Project           : Quantum information, Ex3
! 
! Program name      : MatTest.f03
! 
! Author            : Giorgio Palermo
! 
! Date created      : 20201007
! 
! Purpose           : To test some operations with matrices
! 
! Revision History  :
!
! Date        Author      Ref    Revision (Date in YYYYMMDD format) 
!
! 20201020    G. Palermo         Debug subroutine implementation
! **********************************************************************

module Functions

contains
    subroutine LoopMult(A,B,C)
        real,intent(in)                 ::  A(:,:),B(:,:)
        real,intent(out),allocatable    ::  C(:,:)
        integer                         ::  nn
        integer                         ::  ii,jj,kk
        if(size(A,dim=2) /= size(B,dim=1)) stop "Input array sizes do not match"
        mm=size(A,dim=1)
        nn=size(B,dim=2)
        allocate(C(mm,nn))
        C=0.

        do ii=1,mm
            do jj=1,nn
                do kk=1,size(A,dim=2)
                    C(ii,jj)=C(ii,jj)+A(ii,kk)*B(kk,jj)
                end do      
            end do
        end do
    end subroutine

    subroutine LoopMultColumns(A,B,C)
        real,intent(in)                 ::  A(:,:),B(:,:)
        real,intent(out),allocatable    ::  C(:,:)
        integer                         ::  nn
        integer                         ::  ii,jj,kk
        if(size(A,dim=2) /= size(B,dim=1)) stop "Input array sizes do not match"
        mm=size(A,dim=1)
        nn=size(B,dim=2)
        allocate(C(mm,nn))
        C=0.

        do jj=1,nn
            do kk=1,size(A,dim=2)
                do ii=1,mm
                    C(ii,jj)=C(ii,jj)+A(ii,kk)*B(kk,jj)
                end do      
            end do
        end do
    end subroutine

    subroutine IntrinsicMult(A,B,C)
        real,intent(in)                 ::  A(:,:),B(:,:)
        real,intent(out),allocatable    ::  C(:,:)
        integer                         ::  mm,nn
        if(size(A,dim=2) /= size(B,dim=1)) stop "Input array sizes do not match"
        mm=size(A,dim=1)
        nn=size(B,dim=2)
        allocate(C(mm,nn))
        C=matmul(A,B)

    end subroutine

    subroutine debug0(integer, status)
        integer :: activation, status, oo
        oo=6
        if(activation==1) then
            select case (status)
            case(1)
               write(oo,*) "Primo caso: suca"
            case(2)
                write(oo,*) "Secondo caso: sparisco"
            case default
                write(oo,*) "Ti voglio bene"
            end select
            return
        else
            return
        end if
    end subroutine
end module



program MatTest
use Functions
implicit none
integer :: nn=5, ii, jj, kk, ierror, size, num_args
integer :: act=1, status=0 ! Debug activation variables
integer :: oo=6
real, dimension(:,:), allocatable :: A,B,C,C1,C2
real :: start=0, finish=0, sum=0
real, dimension(:,:), allocatable :: time
character(len=30), dimension(:), allocatable :: args
character(len=50) :: filename
write(oo,*)
write(oo,*) "	*** Matrix multiplication test program ***	"


num_args = command_argument_count()
allocate(args(num_args))
if(num_args==1) then
    ! Case: input filename
    call get_command_argument(1,args(1))
    filename = args(1)
    allocate(time(3,nn))
else if(num_args==2) then
    ! Case: input filename and matrix size
    call get_command_argument(1,args(1))
    call get_command_argument(2,args(2))
    filename = args(1)
    read (args(2),'(I10)') nn
    allocate(time(3,nn))
else if(num_args>=2) then
    write(oo,*) "Error: too many input arguments"
    stop
else
    ! Case: no input arguments
    allocate(time(3,nn))
    write(oo, *) "Enter the filename:"
    read(*,*) filename
end if
filename= trim(filename) // ".txt"
open(unit=1,file=filename,action='write',status='unknown')
write(oo, '("Filename: ",a20," Testing up to size ", i6 )') filename, nn*100
write(oo,*)
write(oo,*)
do ii=1,nn,1
    size=ii*100
    write(oo,'("Running test on size ", i6, " matrix...")') size
    allocate(A(size,size),B(size,size),C(size,size),C1(size,size),C2(size,size))

    call random_number(A)
    call random_number(B)

    call cpu_time(start)
    call LoopMult(A,B,C)
    call cpu_time(finish)
    time(1,ii)=finish-start

    call cpu_time(start)
    call LoopMultColumns(A,B,C1)
    call cpu_time(finish)
    time(2,ii)=finish-start

    call cpu_time(start)
    call IntrinsicMult(A,B,C2)
    call cpu_time(finish)
    time(3,ii)=finish-start
    deallocate(A,B,C,C1,C2)

    write(oo,*) "Done"
    write(oo,*) "       Size   ", "    ByRows[s]  ","     ByCols[s]   ","    Intrinsic[s]   "
    write(oo,*) size,time(1,ii),time(2,ii),time(3,ii)
    write(1,*) size, time(:,ii)
end do
close(1)
call debug0(act,status)

write(oo,*) "End of the program reached, good job!"


end program