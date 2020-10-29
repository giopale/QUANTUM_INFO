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
! 20201026    G. Palermo         New functions in Debug module,
!                                comments added
! **********************************************************************
! This file includes two modules and one program:
!     module FUNCTIONS
!     module DEBUG
!     program MATTEST1
! 
! *********   module FUNCTIONS    **************
! This module contains all the functions used in the program.
! 
!     subroutine LOOPMULT(A,B,C)
!         Arguments:
!             A(:,:),B(:,:)    double precision,intent(in)
!             C(:,:)           double precision,intent(out),allocatable
! 
!         Performs multiplication of double precision matrices
!         A,B and writes the result on C.
!         Multiplication is performed by rows.
! 
!         Preallocation is needed only for A,B (input).
! 
!         The subroutine checks the match of size(A,2) and size(B,1)
!         and stops the execution of the program if the check fails.
! 
!     subroutine LOOPMULTCOLUMNS(A,B,C)
!         Arguments:
!             A(:,:),B(:,:)    double precision,intent(in)
!             C(:,:)           double precision,intent(out),allocatable
! 
!         Performs multiplication of double precision matrices
!         A,B and writes the result on C.
!         Multiplication is performed by columns.
! 
!         Preallocation is needed only for A,B (input).
! 
!         The subroutine checks the match of size(A,2) and size(B,1)
!         and stops the execution of the program if the check fails.
! 
!     subroutine INTRINSICMULT(A,B,C)
!         Arguments:
!             A(:,:),B(:,:)    double precision,intent(in)
!             C(:,:)           double precision,intent(out),allocatable
! 
!         Performs multiplication of double precision matrices
!         A,B and writes the result on C.
!         Multiplication is performed via intrinsic function MATMUL(A,B).
! 
!         Preallocation is needed only for A,B (input).
! 
!         The subroutine checks the match of size(A,2) and size(B,1)
!         and stops the execution of the program if the check fails.


! *********   program MATTEST1    **************
! This program tests the computational performances of different matrix
! multiplication algorithms, by measuring CPU_TIME().
! The test is performed by multiplying two randomly generated square matrices
! using 1:LOOPMULT, 2:LOOPMULTCOLUMNS and 3:INTRINSICMULT.
! The program allows to either choose the matrix size in advance or perform an 
! automatic test. Computation times are measured increasing the matrix size by
! 100 each step up to maximum size.

! I/O interface:
! The program is called on bash by 
!     $ ./MatTest1.out [filename] [size]
! where both [filename] and [size] are optional arguments. The filename must be
! given without extension, since it is added automatically. [size] is the maximum
! matrix size that will be tested.
! If no arguments are provided:
!     - filename will be asked
!     - test will run up to default max_size (500)
! 
! The output is given on file on four columns:
! 
!     SIZE    T1[s]   T2[s]   T3[s]

module Functions

contains
    subroutine LoopMult(A,B,C)
        ! input variables
        double precision,intent(in), allocatable ::  A(:,:),B(:,:)
        double precision,intent(out),allocatable    ::  C(:,:)
        integer                         ::  nn
        integer                         ::  ii,jj,kk
        ! input check
        if(size(A,dim=2) /= size(B,dim=1)) stop "Input array sizes do not match"
        mm=size(A,dim=1)
        nn=size(B,dim=2)
        allocate(C(mm,nn))
        C=0.
        !$omp parallel do
        do ii=1,mm
            do jj=1,nn
                do kk=1,size(A,dim=2)
                    C(ii,jj)=C(ii,jj)+A(ii,kk)*B(kk,jj)
                end do      
            end do
        end do
    end subroutine

    subroutine LoopMultColumns(A,B,C)
        ! input variables
        double precision,intent(in), allocatable                 ::  A(:,:),B(:,:)
        double precision,intent(out),allocatable    ::  C(:,:)
        integer                         ::  nn
        integer                         ::  ii,jj,kk
        ! input check
        if(size(A,dim=2) /= size(B,dim=1)) stop "Input array sizes do not match"
        mm=size(A,dim=1)
        nn=size(B,dim=2)
        allocate(C(mm,nn))
        C=0.
        !$omp parallel do
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
        double precision,intent(in), allocatable                ::  A(:,:),B(:,:)
        double precision,intent(out),allocatable    ::  C(:,:)
        integer                         ::  mm,nn
        ! input check
        if(size(A,dim=2) /= size(B,dim=1)) stop "Input array sizes do not match"
        mm=size(A,dim=1)
        nn=size(B,dim=2)
        allocate(C(mm,nn))
        C=matmul(A,B)

    end subroutine

end module

module debug

    contains
    function CheckDim(activation,A,B,C,message, print)
        ! Checks dimensions for matrix multiplication
        ! 0 = OK, 1 = ERROR
        integer :: activation, a2, b1, c1, c2
        double precision, dimension(:,:), allocatable :: A,B
        double precision, dimension(:,:), allocatable, optional :: C
        character(*), optional :: message
        integer, optional ::print
        logical :: CheckDim
        if(activation==1) then
            a2=size(A,2)
            b1=size(B,1)
            if(present(C))then
                c1=size(C,1)
                c2=size(C,2)
                CheckDim = .not.((a2==b1).and.(c1==a2).and.(c2==b1))
            else
                CheckDim = .not.((a2==b1))
            end if
            if(present(message)) then
                message = "MESSAGE *** "  // message
                write(*,*) message
            end if
            if(present(print)) then
                write(*,*) "CHK:  Matrix multiplication shape check (double): ", CheckDim
            end if
        else
            return
        end if
    end function CheckDim

    function CheckDimInt(activation,A,B,C, message, print)
        ! Checks dimensions for matrix multiplication
        ! 0 = OK, 1 = ERROR
        integer :: activation, a2, b1, c1, c2
        integer, dimension(:,:), allocatable :: A,B
        integer, dimension(:,:), allocatable, optional :: C
        logical :: CheckDimInt
        character(*), optional :: message
        integer, optional :: print
        if(activation==1)then
            CheckDimInt=.true.
            a2=size(A,2)
            b1=size(B,1)
            if(present(C))then
                c1=size(C,1)
                c2=size(C,2)
                CheckDimInt = .not.((a2==b1).and.(c1==a2).and.(c2==b1))
            else
                CheckDimInt = .not.((a2==b1))
            end if
            if(present(message)) then
                message = "MESSAGE *** "  // message
                write(*,*) message
            end if
            if(present(print)) then
                write(*,*) "CHK:  Matrix multiplication shape check (integer): ", CheckDimInt
            end if
        else 
            return
        end if
    end function CheckDimInt

    function CheckEqInteger(activation,A,B,message,print)
        integer :: activation, a1, a2, b1, b2, ii,jj, err=0
        integer,optional :: print
        integer, dimension(:,:), allocatable :: A,B
        logical :: CheckEqInteger
        character(*), optional :: message
        CheckEqInteger = .false.
        if(activation/=1) then
            return
        else
            a1=size(A,1)
            a2=size(A,2)
            b1=size(B,1)
            b2=size(B,2)
            do jj=1,a2
                do ii=1,a1
                    err=err + abs(A(ii,jj)-B(ii,jj))
                end do
            end do
            if(err/=0) then
                CheckEqInteger = .true.
            end if
            if(present(message))then
                message = "MESSAGE *** " // message
                write(*,*) message
            end if
            write(*,*) "CHK:  Array inequality check (integer):", CheckEqInteger
            if(present(print))then
            write(*,*) "CHK:  Total error detected: ", err
            end if
        end if
    end function CheckEqInteger

    function CheckEqDouble(activation,A,B,message, print)
        integer :: activation, a1, a2, b1, b2, ii,jj
        integer, optional :: print
        double precision err
        double precision, dimension(:,:), allocatable :: A,B
        logical :: CheckEqDouble
        character(*),optional :: message
        err = 0.0
        CheckEqDouble = .false.
        if(activation/=1) then
            return
        else
            a1=size(A,1)
            a2=size(A,2)
            b1=size(B,1)
            b2=size(B,2)
            do jj=1,a2
                do ii=1,a1
                    err=err + abs((A(ii,jj)-B(ii,jj))/B(ii,jj))
                end do
            end do
            if(err>=1e-5*size(A,1)*size(A,2)) then
                CheckEqDouble = .true.
            end if
            if(present(message))then
                message = "MESSAGE *** "  // message
                write(*,*) message
            end if
            write(*,*) "CHK:  Array inequality check (double):", CheckEqDouble
            if(present(print))then
            write(*,*) "CHK:  Total error detected: ", err
            end if
        end if
    end function CheckEqDouble

    function CheckTrace(activation,A, message, print)
        integer :: activation,ii
        logical :: CheckTrace 
        double precision, dimension(:,:) :: A
        double precision :: trace =0
        character(*), optional :: message
        integer, optional ::print
        if(size(A,1)/=size(A,2)) then
            CheckTrace = .true.
            if(present(message)) then
            message = "MESSAGE *** " // message
            write(*,*) message
            end if
            write(*,*) "CHK:  Matrix trace > 0 check: ", CheckTrace
            write(*,*) "CHK:  the matrix isn't square."
            return
        else if(size(A,1)==size(A,2)) then
            CheckTrace = .true.
            do ii=1,size(A,1)
                trace=trace + A(ii,ii)
            end do
            if(trace>=0) then
                CheckTrace = .false.
            end if
            if(present(message)) then
            message = "MESSAGE *** " // message
            write(*,*) message
            end if
            if(present(print)) then
                write(*,*) "CHK:  The trace is positive and equal to: ", trace
            end if
        end if
    end function CheckTrace

end module debug



program MatTest
use Functions
use debug
implicit none
! default max_size here! -> nn
integer :: nn=5, ii=0, jj, kk, aux, ierror, size, num_args 
integer :: act=1,status ! Debug activation variables
logical :: deb1, deb2
integer :: oo=6 ! to suppress screen output set oo=something
                ! this won't work for subroutines
double precision, dimension(:,:), allocatable :: A,B,C,C1,C2
double precision :: start=0, finish=0, sum=0
double precision, dimension(:,:), allocatable :: time
character(len=30), dimension(:), allocatable :: args
character(len=50) :: filename, message
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
open(unit=1,file=filename,action='write',status='unknown',iostat=status)
write(*,'("File opening status =" i3)') status ! Checking correct file opening
write(oo, '("Filename: ",a20," Testing up to size ", i6 )') filename, nn*100
write(oo,*)
do ii=1,nn,1
    size=ii*100
    write(oo,'("Running test on size ", i6, " matrix...")') size ! courtesy message
    allocate(A(size,size),B(size,size),C(size,size),C1(size,size),C2(size,size))
    deb1=CheckDim(act,A,B,C) ! Test matrix dimensions for multiplication
    if(deb1) stop
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
    write(1,'(i5,3G15.5)') size, time(1:3,ii)
end do
close(1,iostat=status)
Write(*,'("File exit status =" i3)') status ! Checking correct file closure
write(oo,*) "End of the program reached, good job!"

end program