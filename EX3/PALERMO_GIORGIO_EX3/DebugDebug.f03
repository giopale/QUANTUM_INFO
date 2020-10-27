! **********************************************************************
! Project           : Quantum information, Ex3
! 
! Program name      : DebugDebug.f03
! 
! Author            : Giorgio Palermo
! 
! Date created      : 20201026
! 
! Purpose           : To develop and test the debug module
! 
! Revision History  :
!
! Date        Author      Ref    Revision (Date in YYYYMMDD format) 
!
! **********************************************************************
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


program DebugDebug
    use Debug
implicit none
integer ii,jj, act
integer, dimension(:,:), allocatable :: A, B, C
double precision, dimension(:,:), allocatable :: A1, B1, C1
character(len=70) :: message
logical :: deb
act=1

allocate(A(2,2),B(2,2),C(3,7))
allocate(A1(2,2),B1(2,2),C1(8,2))
do jj=1,size(A,1)
    do ii=1,size(A,2)
        A(ii,jj)= ii+jj
        B(ii,jj) = A(ii,jj)
    end do
end do

message = "Test on A,B integer, equal:"
deb= CheckDimInt(activation=act,A=A,B=B,message=message, print=1)
deb= CheckEqInteger(activation=act,A=A,B=B,print=1)
message = "Test on A,B integer, different"
B(2,2) = 65
deb= CheckDimInt(activation=act,A=A,B=B,message=message, print=1)
message = "Test on A,C integer, different dim"
deb= CheckDimInt(activation=act,A=A,B=C,message=message, print=1)
write(*,*)

do jj=1,size(A1,1)
    do ii=1,size(A1,2)
        A1(ii,jj)= .1 
        B1(ii,jj) = A1(ii,jj)
    end do
end do

message = "Test on A,B double, equal"
deb=CheckDim(activation=1,A=A1,B=B1, message=message,print=1)
deb=CheckEqDouble(act,A1,B1,print=1)
message = "Test on A,B double, different"
B(1,1)=.666
deb=CheckDim(activation=1,A=A1,B=B1, message=message,print=1)
deb=CheckEqDouble(act,A1,B1,print=1)
message = "Test on A,C double, different dim"
deb=CheckDim(activation=1,A=A1,B=C1, message=message,print=1)

end program DebugDebug