module ModDmat
    type dmatrix
        integer, dimension(2) :: N = (/ 0,0 /)
        double complex, dimension(:,:), allocatable :: elem
        double complex, dimension(:,:), allocatable :: evec
        double precision, dimension(:), allocatable :: eval
        double complex :: Trace
        double complex :: Det
    end type dmatrix

    interface operator (.Adj.)
            module procedure Adj
    end interface 

    interface operator (.Init.)
        module procedure InitUniDCplx
    end interface

    interface operator (.evalh.)
        module procedure DmatHermEv
    end interface

    contains
        
        function InitUniDouble(dmat) 
            ! Initializes a (m,n) complex matrix with 
            ! real and imaginary part taken from [0,1]
            ! uniform distributions
            implicit none
            type(dmatrix), intent(in) :: dmat
            type(dmatrix) :: InitUniDouble
            integer :: jj,sd=4
            integer, dimension(:), allocatable :: seed
            double complex, dimension(dmat%N(1),1) :: X

            if(dmat%N(1)<1 .or. dmat%N(2)<1) then
                ! Check for positive matrix shape
                print*, "*** ERROR in InitUni: matrix shape not defined"
                print*, "Program terminated"
                stop
            else
                call random_seed()
                call random_seed(size = sd)
                allocate(seed(sd))
                call random_seed(get=seed)
                do jj=1,dmat%N(2)
                    ! '1' stands for uniform
                    call clarnv(1, seed, 2*dmat%N(1), dmat%elem(:,jj) )
                end do
                InitUniDouble = dmat
                return
            end if
        end function InitUniDouble

        function InitUniDCplx(herm_1)
            ! Initializes a (n,n) complex hermitian matrix with 
            ! real and imaginary part taken from [0,1]
            ! uniform distributions.
            ! Local scalars
            implicit none
            integer :: Nmat
            integer :: seed_dim =4
            integer :: info_zlaghe=0
            ! Local arrays
            integer, dimension(:), allocatable :: seed !Seed for random matrix generator
            double precision, dimension(:), allocatable :: Diag_Rnd
            double complex, dimension(:), allocatable :: Work
            ! Local custom types
            type(dmatrix), intent(in) :: herm_1
            type(dmatrix):: InitUniDCplx

            Nmat = herm_1%N(1)
            allocate(Diag_Rnd(Nmat))
            allocate(Work(2*Nmat))
            allocate(InitUniDCplx%elem(Nmat,Nmat))

            call random_seed()
            call random_seed(size=seed_dim)
            allocate(seed(seed_dim))
            call random_seed(get=seed)

            call dlarnv(1,seed,Nmat,Diag_Rnd) ! On exit, the seed is updated, but to be shure...
            call random_seed()
            call random_seed(get=seed)
            call zlaghe(Nmat, Nmat-1,Diag_Rnd,herm_1%elem, Nmat, seed, work, info_zlaghe)
            InitUniDCplx = herm_1

        end function InitUniDCplx



        function DmatHermEv(herm_1)
            ! Computes eigenvalues of an hermitian matrix
            ! contained in a dmatrix type.
            ! Note that the shape of the matrix, the matrix
            ! itself and the vector containing the eigenvalues 
            ! must be initialized in advance.
            implicit none
            ! Local scalars
            character(len=1) :: jobz = "N", uplo = "U"
            integer :: Nmat, lda, lwork, rworkdim, info
            ! Local arrays
            double precision, dimension(:), allocatable :: w
            double complex, dimension(:), allocatable :: work
            double precision, dimension(:), allocatable :: rwork
            ! Local custom types
            type(dmatrix), intent(in) :: herm_1
            type(dmatrix) :: DmatHermEv

            Nmat = herm_1%N(1)
            lda = Nmat
            DmatHermEv = herm_1
            allocate(DmatHermEv%eval(Nmat))
            rworkdim = max(1,3*Nmat-2)
            allocate(w(Nmat))
            allocate(work(1))
            allocate(rwork(rworkdim))
            call zheev(jobz,uplo,Nmat,DmatHermEv%elem,lda,w,work,-1,rwork,info) ! call to determine LWork
            lwork=int(work(1))
            deallocate(work)
            allocate(work(lwork))
            call zheev(jobz,uplo,Nmat,DmatHermEv%elem,lda,w,work,lwork,rwork,info) !call to solve the eigenproblem
            
            DmatHermEv%eval=w

        end function DmatHermEv

        subroutine Tr(dmat)
            ! Computes the trace and assigns it
            ! to the "Trace" field of the input object
            ! of type dmatrix
            type(dmatrix) :: dmat
            integer ::ii
            if(dmat%N(1)==dmat%N(2) .and. dmat%N(1)>0 .and. dmat%N(2) >0) then
                do ii=1,dmat%N(1)
                    dmat%Trace=dmat%Trace +dmat%elem(ii,ii)
                end do
                return
            else
                print*, "*** ERROR in Tr: matrix dimensions must be positive and equal"
                return
            end if
        end subroutine Tr

        function Adj(dmat)
            ! computes the adjoint and passes it
            ! as output
            type(dmatrix),intent(in) :: dmat
            type(dmatrix) ::Adj
            Adj%Trace=conjg(dmat%Trace)
            Adj%Det=conjg(dmat%Det)
            Adj%N=(/ dmat%N(2), dmat%N(1) /)
            allocate(Adj%elem(Adj%N(1),Adj%N(2)))
            Adj%elem=dmat%elem
            Adj%elem=conjg(Adj%elem)
            Adj%elem=transpose(Adj%elem)
            return
        end function Adj

        subroutine MatToFile(dmat,filename)
            ! Writes dmatrix type to txt file
            implicit none
            type(dmatrix) :: dmat
            character(*) :: filename
            character(30) :: auxname
            integer :: jj
            auxname = trim(filename // ".txt")
            open(unit=10,file=auxname,action='write',status='unknown')
            do jj=1,dmat%N(1)
                write(10,*) dmat%elem(jj,:)
            end do
            close(10)
        end subroutine MatToFile

end module ModDmat