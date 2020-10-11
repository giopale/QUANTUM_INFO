! **********************************************************************
! Project           : Quantum information, Ex1
! 
! Program name      : MatTest.f90
! 
! Author            : Giorgio Palermo
! 
! Date created      : 20201007
! 
! Purpose           : Test some operations with matrices
! 
! Revision History  :
!
! Date        Author      Ref    Revision (Date in YYYYMMDD format) 
!
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

end module



program MatTest
use Functions
implicit none
integer :: nn=0, ii, jj, kk, ierror, size
real, dimension(:,:), allocatable :: A,B,C,C1,C2
real :: start=0, finish=0, sum=0
real, dimension(:,:), allocatable :: time

CALL chdir("/Users/gpalermo/Google\ Drive/QUANTUM\ INFO/QUANTUM_INFO/EX1 ")

write(*,*) "	*** Matrix multiplication test program ***	"

do
    write(*,*) "Enter an integer number or press enter to run an automatic test"
    read(*,'(i10)',iostat=ierror) nn

    if ( ierror == 0 ) then
      exit
    endif
    write(*,*) 'An error occured - please try again'

enddo


if(nn/=0) then
    allocate(A(nn,nn),B(nn,nn),C(nn,nn),C1(nn,nn),C2(nn,nn))

    call random_number(A)
    call random_number(B)


    call cpu_time(start)
    call LoopMult(A,B,C)
    call cpu_time(finish)
    print '("Time = ",f6.5," seconds.")',finish-start

    call cpu_time(start)
    call LoopMultColumns(A,B,C1)
    call cpu_time(finish)
    print '("Time = ",f6.5," seconds.")',finish-start

    call cpu_time(start)
    call IntrinsicMult(A,B,C2)
    call cpu_time(finish)
    print '("Time = ",f6.5," seconds.")',finish-start

    deallocate(A,B,C,C1,C2)
else
    allocate(time(3,10))
    open(unit=1,file="time_results.txt",action='write',status='unknown')
    do ii=1,3,1
        size=ii*100
        print '("Running test on size ", i6, " matrix...")', size
        allocate(A(size,size),B(size,size),C(size,size),C1(size,size),C2(size,size))

        call random_number(A)
        call random_number(B)


        call cpu_time(start)
        call LoopMult(A,B,C)
        call cpu_time(finish)
        time(1,ii)=finish-start
        ! print '("Time Loop= ",f6.5," seconds.")',finish-start

        call cpu_time(start)
        call LoopMultColumns(A,B,C1)
        call cpu_time(finish)
        time(2,ii)=finish-start
        ! print '("Time = ",f6.5," seconds.")',finish-start

        call cpu_time(start)
        call IntrinsicMult(A,B,C2)
        call cpu_time(finish)
        time(3,ii)=finish-start
        ! print '("Time = ",f6.5," seconds.")',finish-start
        ! print '("Times [Loop,ColumnsLoop,Intrinsic], seconds = ",f6.5,f6.5,f6.5)',time(1,ii),time(2,ii),time(3,ii)
        deallocate(A,B,C,C1,C2)

        print*, "Done"

        print*,size,time(1,ii),time(2,ii),time(3,ii)
        write(1,*) size, time(:,ii)
    end do
    close(1)
end if


end program