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
    end module Functions

module cosmod
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
    end module cosmod

program cos
    use Functions
    use cosmod
    implicit none

    integer grid_dim
    integer, dimension(:), allocatable :: grid

    grid_dim=7
    allocate(grid(grid_dim))
    

end program cos

