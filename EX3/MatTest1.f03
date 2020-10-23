! **********************************************************************
! Project           : Quantum information, Ex3
! 
! Program name      : MatTest1.f03
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
! 20201021    G. Palermo         Change all reals to doubles
! **********************************************************************

module Functions

contains
    subroutine LoopMult(A,B,C)
        ! input variables
        double precision,intent(in)                 ::  A(:,:),B(:,:)
        double precision,intent(out),allocatable    ::  C(:,:)
        integer                         ::  nn
        integer                         ::  ii,jj,kk
        ! input check
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

    subroutine LoopMultColumns(A,B,C)]
        ! input variables
        double precision,intent(in)                 ::  A(:,:),B(:,:)
        double precision,intent(out),allocatable    ::  C(:,:)
        integer                         ::  nn
        integer                         ::  ii,jj,kk
        ! input check
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
        ! input variables
        double precision,intent(in)                 ::  A(:,:),B(:,:)
        double precision,intent(out),allocatable    ::  C(:,:)
        integer                         ::  mm,nn
        ! input check
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

module debug
    contains
    function CheckDim(activation,A,B,C)
        ! Checks dimensions for matrix multiplication
        ! 0 = OK, 1 = ERROR
        integer :: activation, a2, b1, c1, c2
        double precision, dimension(:,:), allocatable :: A,B,C
        logical :: CheckDim
        a2=size(A,2)
        b1=size(B,1)
        c1=size(C,1)
        c2=size(C,2)
        if(allocated(C))then
            CheckDim = .not.((a2==b1).and.(c1==a2).and.(c2==b1))
        else
            CheckDim = .not.((a2==b1))
        end if
        if(CheckDim) then
            write(*,*) "DEBUG *** Matrix multiplication shape check: ", CheckDim
        end if
    end function CheckDim

    function CheckTrace(activation,A)
        integer :: activation,ii
        logical :: CheckTrace 
        double precision, dimension(:,:) :: A
        double precision :: trace =0
        if(size(A,1)/=size(A,2)) then
            CheckTrace = .true.
            write(*,*) "DEBUG *** Matrix trace > 0 check: ", CheckTrace
            write(*,*) "DEBUG *** the matrix isn't square."
        else if(size(A,1)==size(A,2)) then
            CheckTrace = .true.
            do ii=1,size(A,1)
                trace=trace + A(ii,ii)
            end do
            if(trace>=0) then
                CheckTrace = .false.
            end if
        end if
    end function CheckTrace




end module debug



program MatTest
use Functions
use debug
implicit none
! default max_size here! -> nn
integer :: nn=5, ii, jj, kk, aux, ierror, size, num_args 
integer :: act=1, status=0 ! Debug activation variables
logical :: deb1, deb2
integer :: oo=6 ! to suppress screen output set oo=something
                ! this won't work for subroutines
double precision, dimension(:,:), allocatable :: A,B,C,C1,C2
double precision :: start=0, finish=0, sum=0
double precision, dimension(:,:), allocatable :: time
character(len=30), dimension(:), allocatable :: args
character(len=50) :: filename
write(oo,*)
write(oo,*) "	*** Matrix multiplication test program ***	"

! Input interface: initializes nn and filename
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
    read (args(2),'(I10)') aux
    nn=aux/100
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

! Testing section 
open(unit=1,file=filename,action='write',status='unknown')
write(oo, '("Filename: ",a20," Testing up to size ", i6 )') filename, nn*100
write(oo,*)
do ii=1,nn,1
    size=ii*100
    write(oo,'("Running test on size ", i6, " matrix...")') size ! courtesy message
    allocate(A(size,size),B(size,size),C(size,size),C1(size,size),C2(size,size))
    ! deb1=CheckDim(act,A,B,C)
    call random_number(A)
    call random_number(B)
    ! deb2= CheckTrace(act,A)
    
    call cpu_time(start)
    call LoopMult(A,B,C)    ! by rows
    call cpu_time(finish)
    time(1,ii)=finish-start

    call cpu_time(start)
    call LoopMultColumns(A,B,C1) ! by columns
    call cpu_time(finish)
    time(2,ii)=finish-start

    call cpu_time(start)
    call IntrinsicMult(A,B,C2)  ! Intrinsic
    call cpu_time(finish)
    time(3,ii)=finish-start
    deallocate(A,B,C,C1,C2)

    ! print results on screen at each iteration
    write(oo,*) "Done"
    write(oo,*) "Size   ", "    ByRows[s]  ","     ByCols[s]   ","    Intrinsic[s]   "
    write(oo,'(i5,3G15.5)') size,time(1,ii),time(2,ii),time(3,ii)
    write(1,*) size, time(:,ii)
end do
close(1)

write(oo,*) "End of the program reached, good job!"


end program