module Functions

    contains
        subroutine LoopMult(A,B,C,mm,nn)
            ! input variables
            double precision,intent(in)                 ::  A(:,:),B(:,:)
            double precision,intent(out), dimension(mm,nn)    ::  C
            integer, intent(in)                        ::  mm, nn
            integer                        ::  ii,jj,kk
            C=0.

            do ii=1,mm
                do jj=1,nn
                    do kk=1,size(A,dim=2)
                        C(ii,jj)=C(ii,jj)+A(ii,kk)*B(kk,jj)
                    end do      
                end do
            end do
        end subroutine

        subroutine LoopMultColumns(A,B,C,mm,nn)
            ! input variables
            double precision,intent(in)                 ::  A(:,:),B(:,:)
            double precision,intent(out), dimension(mm,nn)    ::  C
            integer, intent(in)                        ::  mm, nn
            integer                        ::  ii,jj,kk
            C=0.

            do jj=1,nn
                do kk=1,size(A,dim=2)
                    do ii=1,mm
                        C(ii,jj)=C(ii,jj)+A(ii,kk)*B(kk,jj)
                    end do      
                end do
            end do
        end subroutine

        subroutine IntrinsicMult(A,B,C,mm,nn)
            ! input variables
            double precision,intent(in)                 ::  A(:,:),B(:,:)
            double precision,intent(out), dimension(mm,nn)    ::  C
            integer, intent(in)                        ::  mm, nn
            integer                        ::  ii,jj,kk
            C=0.
            C=matmul(A,B)

        end subroutine

end module


module MatrixMultip
    contains
    subroutine MatTest(filename, mat_dim_in, verbose, result)
        use Functions
        use Debug
        implicit none
        ! Local scalars
        character(*), intent(in) :: filename
        character(*), intent(in) :: verbose
        character(len=50) :: filename_in, filename_alpha, filename_bravo
        integer, intent(in) :: mat_dim
        integer :: status, mm, nn, size
        integer :: oo=6 ! to suppress screen output set oo=something
                        ! this won't work for subroutines
        double precision :: start=0, finish=0
        ! Local Arrays
        double precision, dimension(4) :: result
        double precision, dimension(:,:), allocatable :: A,B,C,C1,C2
        double precision, dimension(1,3) :: time

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
        open(unit=10,file=filename_in,action='write',status='unknown',iostat=status)
        open(unit=20,file=filename_alpha,action='write',status='unknown',iostat=status)
        open(unit=30,file=filename_bravo,action='write',status='unknown',iostat=status)
        ! write(*,'("File opening status =" i3)') status ! Checking correct file opening
        ! write(oo, '("Filename: ",a20," Testing up to size ", i6 )') filename, nn*100
        ! write(oo,*)
        write(oo,*) "Running test..." ! courtesy message
        write(oo,*) "Size   ", "    ByRows[s]  ","     ByCols[s]   ","    Intrinsic[s]   "
        
        allocate(A(mat_dim_in,mat_dim_in),B(mat_dim_in,mat_dim_in))
        allocate(C(mat_dim_in,mat_dim_in),C1(mat_dim_in,mat_dim_in),C2(mat_dim_in,mat_dim_in))
        call random_number(A)
        call random_number(B)
        mm =size(A,1)
        nn =size(B,2)
        
        call cpu_time(start)
        call LoopMult(A,B,C,mm,nn)    ! by rows
        call cpu_time(finish)
        time(1,1)=finish-start

        call cpu_time(start)
        call LoopMultColumns(A,B,C1,mm,nn) ! by columns
        call cpu_time(finish)
        time(1,2)=finish-start

        call cpu_time(start)
        call IntrinsicMult(A,B,C2,mm,nn)  ! Intrinsic
        call cpu_time(finish)
        time(1,3)=finish-start
        deallocate(A,B,C,C1,C2)          
            
        write(oo,'(i5,3G15.5)') mat_dim_in,time(1,1),time(1,2),time(1,3)
        write(10,'(i10, G15.5)') mat_dim_in, time(1,1)
        write(20,'(i10, G15.5)') mat_dim_in, time(1,2)
        write(30,'(i10, G15.5)') mat_dim_in, time(1,3)
        result(1) = dble(mat_dim_in)
        result(2) = time(1,1)
        result(3) = time(1,2)
        result(4) = time(1,3)
        write(oo,*) "Done"
        close(10,iostat=status)
        close(20,iostat=status)
        close(20,iostat=status)

    end subroutine MatTest

end module MatrixMultip