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

program MatTest

integer :: nn=0
real, dimension(:,:), allocatable :: A,B,C, At,Bt,Ct
real :: start=0, finish=0, sum=0

real, dimension(:,:), allocatable :: D, E, F

write(*,*) "	*** Matrix multiplication test program ***	"
do
    write(*,*) "Enter an integer number"
    read(*,'(i10)',iostat=ierror) nn

    if ( ierror == 0 ) then
      exit
    endif
    write(*,*) 'An error occured - please try again'

enddo

allocate(A(nn,nn))
allocate(B(nn,nn))
allocate(C(nn,nn))

call random_number(A)
call random_number(B)


!   TEST STUFF

! allocate(D(2,2))
! allocate(E(2,2))
! allocate(F(2,2))
! D(1,1)=1
! D(1,2)=2
! D(2,1)=3
! D(2,2)=4
! E(1,1)=3
! E(1,2)=6
! E(2,1)=7
! E(2,2)=1

! call cpu_time(start)
! do i=1,2       
!     do j=1,2
!         do k=1,2
!             sum=sum+D(i,k)*E(k,j)
!             if (k==2) then
!                 F(i,j) = sum
!                 sum=0
!                 print*, F(i,j)
!             end if
!         end do
!     end do
! end do
! call cpu_time(finish)

! Direct multiplication
call cpu_time(start)
do i=1,nn      
    do j=1,nn
        do k=1,nn
            sum=sum+A(i,k)*B(k,j)
            if (k==nn) then
                C(i,j) = sum
                sum=0
            end if
        end do
    end do
end do
call cpu_time(finish)
print '("Time = ",f6.5," seconds.")',finish-start
! Transposed multiplication
call cpu_time(start)
do j=1,nn      
    do i=1,nn
        do k=1,nn
            sum=sum+A(j,k)*B(i,j)
            if (k==nn) then
                C(i,j) = sum
                sum=0
            end if
        end do
    end do
end do
! C=transpose(C)
call cpu_time(finish)
print '("Time = ",f6.5," seconds.")',finish-start
deallocate(A)
deallocate(B)
deallocate(C)
end program