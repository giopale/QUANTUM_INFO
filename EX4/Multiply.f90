module Functions

    contains
        subroutine LoopMult(A,B,C)
            implicit none 
            ! input variables
            double precision,intent(in)                 ::  A(:,:),B(:,:)
            double precision,intent(out),allocatable    ::  C(:,:)
            integer                         ::  mm,nn
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

        subroutine LoopMultColumns(A,B,C)
            implicit none
            ! input variables
            double precision,intent(in)                 ::  A(:,:),B(:,:)
            double precision,intent(out),allocatable    ::  C(:,:)
            integer                         ::  mm,nn
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
            implicit none
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

module MatrixMultip
    contains
    subroutine MatTest(filename, size, verbose, result)
        use Functions
        use Debug
        implicit none
        ! Local scalars
        character(*), intent(in) :: filename
        character(*), intent(in) :: verbose
        character(len=50) :: filename_in, filename_alpha, filename_bravo, message
        integer, intent(in) :: size
        integer :: nn=5, ii=0, status
        integer :: oo=6 ! to suppress screen output set oo=something
                        ! this won't work for subroutines
        double precision :: start=0, finish=0, sum=0
        ! Local Arrays
        double precision, dimension(4) :: result
        double precision, dimension(:,:), allocatable :: A,B,C,C1,C2
        double precision, dimension(1,3) :: time
        character(len=30), dimension(:), allocatable :: args

        ! Select if verbose
        select case(verbose)
        case("y")
            oo=6
        case("n")
            oo=3456
        case default
            oo=3456
        end select

        write(oo,*)
        write(oo,*) "   *** Matrix multiplication test program ***  "
        
        ! filename_in = filename
        filename_in = filename //"_ByRows"// ".txt"
        filename_alpha = filename //"_ByCols" // ".txt"
        filename_bravo = filename // "_Intrinsic" // ".txt"

        ! Testing section 
        open(unit=10,file=filename_in,action='write',position="append",status='unknown',iostat=status)
        open(unit=20,file=filename_alpha,action='write',position="append",status='unknown',iostat=status)
        open(unit=30,file=filename_bravo,action='write',position="append",status='unknown',iostat=status)
        ! write(*,'("File opening status =" i3)') status ! Checking correct file opening
        ! write(oo, '("Filename: ",a20," Testing up to size ", i6 )') filename, nn*100
        ! write(oo,*)
        write(oo,*) "Running test..." ! courtesy message
        write(oo,*) "Size   ", "    ByRows[s]  ","     ByCols[s]   ","    Intrinsic[s]   "
        
        allocate(A(size,size),B(size,size),C(size,size),C1(size,size),C2(size,size))
        call random_number(A)
        call random_number(B)
        
        call cpu_time(start)
        call LoopMult(A,B,C)    ! by rows
        call cpu_time(finish)
        time(1,1)=finish-start

        call cpu_time(start)
        call LoopMultColumns(A,B,C1) ! by columns
        call cpu_time(finish)
        time(1,2)=finish-start

        call cpu_time(start)
        call IntrinsicMult(A,B,C2)  ! Intrinsic
        call cpu_time(finish)
        time(1,3)=finish-start
        deallocate(A,B,C,C1,C2)          
            
        write(oo,'(i5,3G15.5)') size,time(1,1),time(1,2),time(1,3)
        write(10,'(i10, G15.5)') size, time(1,1)
        write(20,'(i10, G15.5)') size, time(1,2)
        write(30,'(i10, G15.5)') size, time(1,3)
        result(1) = dble(size)
        result(2) = time(1,1)
        result(3) = time(1,2)
        result(4) = time(1,3)
        write(oo,*) "Done"
        close(10,iostat=status)
        close(20,iostat=status)
        close(20,iostat=status)

    end subroutine MatTest

    subroutine ReadGrid(filename,grid_dim, grid)
        implicit none
        character(*), intent(in) :: filename
        character(len=100) :: msg
        integer, intent(in) :: grid_dim
        integer, dimension(:), allocatable,intent(out) :: grid
        integer :: ii=1, ios

        open(unit=100,file=filename,iostat=ios,iomsg=msg)
        if(ios/=0) then
            write(*,*) msg
            stop
        end if

        allocate(grid(grid_dim))
        read(100,*,iostat=ios, end=997,iomsg=msg) grid 
        997 print*, ios, msg, grid(:)
        close(100)

    end subroutine ReadGrid

    subroutine WriteGrid(filename,min_dim,max_dim,grid_dim)
        implicit none
        character(*), intent(in) :: filename
        character(len=100) :: msg
        integer, intent(in) :: min_dim, max_dim
        integer, intent(out) :: grid_dim
        integer :: ii=1, max_idx=1000, nn=1, ios
        integer, dimension(:), allocatable :: vect
        real :: exp_factor = 1.5
        open(unit=10,file=filename,status='unknown', action='write', iostat=ios, iomsg=msg )
        if(ios/=0) then
            print*,  msg
            stop
        end if
        allocate(vect(max_idx))
        vect=0
        vect(1)=min_dim
        do while(vect(ii)<=max_dim .and. ii<=max_idx )
            ii=ii+1
            vect(ii)=floor(exp_factor*vect(ii-1))
            ! print*, ii<=max_idx, ii, max_idx
        end do
        ii=1
        do while(vect(ii)/=0 .and. vect(ii)<= max_dim)
            write(10,*) vect(ii)
            ii=ii+1
        end do
        grid_dim=ii-1
        close(10)
    end subroutine WriteGrid

    subroutine EmptyFile(filename)
        character(*), intent(in) :: filename
        open(unit=99,file=filename ,action='write',status='old')
        write(10,'(a,$)')
        close(99)
    end subroutine EmptyFile

end module MatrixMultip

program Multiply
    use MatrixMultip
    implicit none
    double precision, dimension(4) :: result
    character(len=30), dimension(:), allocatable :: args
    integer :: ios, num_args, ii, grid_dim
    integer :: max_dim, min_dim
    integer, dimension(:), allocatable :: grid


    num_args = command_argument_count()
    allocate(args(num_args)) 
    call get_command_argument(1,args(1))
    call get_command_argument(2,args(2))
    read(args(1),'(i10)') min_dim
    read(args(2),'(i10)') max_dim
    
    call WriteGrid("grid.dat", min_dim, max_dim, grid_dim)
    grid_dim=6
    ! call ReadGrid(args(1),grid_dim,grid)
    ! call EmptyFile("time_ByCols.txt")
    ! call EmptyFile("time_ByRows.txt")
    ! call EmptyFile("time_Intrinsic.txt")
    ! do ii=1,grid_dim,1
    !     call MatTest("time", grid(ii), "n", result)
    !     write(*,'(4G15.5)') grid(ii), result(2:4)
    ! end do
    
end program Multiply































