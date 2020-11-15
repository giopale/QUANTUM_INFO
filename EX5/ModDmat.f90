module ModDmat
    type dmatrix
        integer, dimension(2) :: N = (/ 0,0 /)
        double complex, dimension(:,:), allocatable :: elem
        double complex, dimension(:,:), allocatable :: evec
        double precision, dimension(:), allocatable :: eval
        double complex :: Trace
        double complex :: Det
    end type dmatrix

    type histogram
        integer :: Nbins, entries
        double precision ::lower, upper
        integer, dimension(:), allocatable :: h
        double precision, dimension(:), allocatable :: hnorm
        double precision, dimension(:), allocatable :: bounds, bincenters
    end type histogram

    interface operator (.Adj.)
            module procedure Adj
    end interface 

    interface operator (.Init.)
        module procedure InitHerm
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

        function InitHerm(herm_1,Nmat)
            ! Initializes a (n,n) complex hermitian matrix with 
            ! real and imaginary part taken from [0,1]
            ! uniform distributions.
            ! call with: herm = herm .init. Nmat
            implicit none
            ! Local scalars
            integer, intent(in) :: Nmat
            integer :: seed_dim =4
            integer :: info_zlaghe=0
            ! Local arrays
            integer, dimension(:), allocatable :: seed !Seed for random matrix generator
            double precision, dimension(:), allocatable :: Diag_Rnd
            double complex, dimension(:), allocatable :: Work
            ! Local custom types
            type(dmatrix), intent(in) :: herm_1
            type(dmatrix):: InitHerm

            ! Nmat = herm_1%N(1)
            allocate(Diag_Rnd(Nmat))
            allocate(Work(2*Nmat))
            ! allocate(InitHerm%elem(Nmat,Nmat))

            call random_seed()
            call random_seed(size=seed_dim)
            allocate(seed(seed_dim))
            call random_seed(get=seed)

            call dlarnv(1,seed,Nmat,Diag_Rnd) ! On exit, the seed is updated, but to be shure...
            call random_seed()
            call random_seed(get=seed)
            call zlaghe(Nmat, Nmat-1,Diag_Rnd,herm_1%elem, Nmat, seed, work, info_zlaghe)
            InitHerm = herm_1

        end function InitHerm

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

        function SpAvg(herm_1)
            implicit none
            ! Local custom types
            type(dmatrix) :: herm_1
            ! Local scalars
            integer :: Nmat,sp_size
            double precision :: s_avg, sp_sum
            ! Local arrays
            double precision, dimension(:), allocatable :: values, sp_0, sp
            double precision, dimension(herm_1%N(1)-1) :: SpAvg
            

            allocate(values(herm_1%N(1)))
            allocate(sp_0(herm_1%N(1)))
            allocate(sp(size(sp_0-1)))
            values = herm_1%eval
            sp_0=values
            values=eoshift(values, shift=-1)
            sp_0 = sp_0-values
            sp=sp_0(2:)
            sp_size = size(sp)
            sp_sum = sum(sp)
            s_avg = sp_sum/sp_size
            sp=sp(:)/s_avg

            SpAvg = sp
        end function SpAvg

        function SpAvgDble(diag)
            implicit none
            ! Local custom types

            ! Local scalars
            integer :: Nmat,sp_size
            double precision :: s_avg, sp_sum
            ! Local arrays
            double precision, dimension(:), intent(in) :: diag
            double precision, dimension(:), allocatable :: values, sp_0, sp
            double precision, dimension(size(diag)-1) :: SpAvgDble
            
            Nmat=size(diag)
            allocate(values(Nmat))
            allocate(sp_0(Nmat))
            allocate(sp(Nmat-1))
            values = diag
                ! print*, "-2 values: ", values(1:3)
            sp_0=values
            values=eoshift(values, shift=-1)
                ! print*, "-1 shifted values: ", values(1:3)
            sp_0 = sp_0-values
                ! print*, "0 sp_0: ", sp_0(1:3)
            ! if (isnan(sum(sp_0))) stop '"sp" is a NaN'
            sp=sp_0(2:)
                ! print*, "1 sp: ", sp(1:3)
                ! print*, "max(sp): ", maxval(sp)
            ! if (isnan(sum(sp))) stop '"sp" is a NaN first'
            sp_size = size(sp)
                ! print*, "2 sp_size: ", sp_size
            ! if (size(sp)==0) stop 'size(sp)==0'
            sp_sum = sum(sp)
                ! print*, "3 sp_sum: ", sp_sum
            ! if (isnan(sp_sum)) stop '"sum(sp)" is Nan'
            s_avg = sp_sum/sp_size
                ! print*, "4 s_avg :", s_avg
            ! if (isnan(s_avg)) stop '"s_avg" is Nan'
            sp=sp(:)/s_avg
                ! print*, "5 sp :", sp(1:3)
                ! print*, "6 sum(sp) :", sum(sp)
            ! if (isnan(sum(sp))) stop '"sp" is a NaN'
            ! Write(*,*)

            SpAvgDble = sp
        end function SpAvgDble

        function LoadArray(filename,nrows,ncols)
            implicit none
                ! Local scalars
            character(*), intent(in) :: filename
            character(len=100) :: msg
            integer :: Ncols, Nrows, ios
            ! Local arrays
            double precision, dimension(:,:), allocatable :: indata
            double precision, dimension(Nrows,Ncols) :: LoadArray

            allocate(indata(Nrows,Ncols))
            open(unit=129, file=filename, iostat=ios, iomsg=msg)
            if(ios/=0) print*, msg
            
            read(129,*,iostat=ios,end=997,iomsg=msg) indata
            997 continue
            if (ios==0) write(*,*) 'File "', trim(filename), '" succesfully loaded...'
            close(129)
            LoadArray=indata

        end function LoadArray

        function InitHisto(low,up,Nb)
            implicit none
            ! Local scalars
            double precision, intent(in) ::low,up
            integer, intent(in) :: Nb
            double precision :: step
            integer :: ii, Nentries
            ! Local arrays
            ! double precision, dimension(:), intent(in) :: indata
            double precision, dimension(:), allocatable :: h_input_aux
            integer, dimension(:), allocatable :: h_input

            ! Local custom types
            type(histogram) :: InitHisto

            InitHisto%Nbins = Nb
            InitHisto%lower = low
            InitHisto%upper = up
            InitHisto%entries=0
            allocate(InitHisto%bounds(Nb+1))
            allocate(InitHisto%bincenters(Nb))
            allocate(InitHisto%h(Nb))
            allocate(InitHisto%hnorm(Nb))
            InitHisto%bounds(1)=low
            InitHisto%h=0
            InitHisto%hnorm=0
            step=(up-low)/Nb
            do ii=1,Nb
                InitHisto%bounds(ii+1)=low + step*ii
                InitHisto%bincenters(ii)= InitHisto%bounds(ii) + step/2
            end do
            
        end function InitHisto

        function FillHisto(h,indata,norm)
            implicit none
            ! Local scalars
            integer ::  Nentries
            double precision :: step
            character(*),intent(in) :: norm
            ! Local arrays
            double precision, dimension(:), intent(in) :: indata
            double precision, dimension(:), allocatable :: h_input_aux
            integer, dimension(:), allocatable :: h_input

            ! Local custom types
            type(histogram),intent(in) :: h
            type(histogram) :: FillHisto

            FillHisto = h
            step = (h%upper - h%lower)/h%Nbins
            Nentries=size(indata)
            allocate(h_input_aux(Nentries),h_input(Nentries))
            h_input_aux= indata
            h_input_aux = ceiling(h_input_aux/step)
            h_input=int(h_input_aux)
            FillHisto%h(h_input) = FillHisto%h(h_input)+1
            FillHisto%entries = sum(FillHisto%h)
            if(norm=="y".and. FillHisto%entries/=0) then
                  FillHisto%hnorm=dble(FillHisto%h)/FillHisto%entries*step
            end if
            print*, "Check fill 5"
        end function FillHisto

        subroutine FillH(h1,indata,norm)
            implicit none
            ! Local scalars
            integer ::  Nentries
            double precision :: step
            character(*),intent(in) :: norm
            ! Local arrays
            double precision, dimension(:), intent(in) :: indata
            double precision, dimension(:), allocatable :: h_input_aux
            integer, dimension(:), allocatable :: h_input

            ! Local custom types
            type(histogram),intent(inout) :: h1

            step = (h1%upper - h1%lower)/h1%Nbins
            Nentries=size(indata)
            allocate(h_input_aux(Nentries),h_input(Nentries))
            h_input_aux= indata
            h_input_aux = ceiling(h_input_aux/step)
            h_input=int(h_input_aux)
            h1%h(h_input) = h1%h(h_input)+1
            h1%entries = sum(h1%h)
            if(norm=="y".and. h1%entries/=0) then
                  h1%hnorm=dble(h1%h)/h1%entries*step
            end if
        end subroutine FillH

        integer*2 function compar( a, b )
            INTEGER*4 a, b
            if ( a <= b ) compar = 1
            if ( a == b ) compar = 0
            if ( a >= b ) compar = -1
            return
        end

        subroutine HistoToFile(h1,filename)
            implicit none
            ! Local scalars
            character(*), intent(in) :: filename
            character(len=200) :: msg, dummy
            integer :: ios,ii, upper
            ! Local arrays
            
            ! Local custom types
            type(histogram),intent(inout) :: h1

            upper = h1%Nbins
            print*, "Check 1"

            open(unit=66, file=filename, status='unknown', iostat=ios, iomsg=msg)
            ! if (ios/=0) then
            !     print*, msg
            !     stop
            ! end if
            do ii=1,upper,1
                print*, ii
                write(*,'(2G15.5)') h1%bincenters(ii), h1%hnorm(ii)
            end do
            close(66,iostat=ios, iomsg=msg)
            ! if (ios/=0) then
            !     print*, msg
            !     stop
            ! else if(ios==0) then
            !     dummy= "File """// trim(filename)//""" succesfully created..."
            !     write(*,*) trim(dummy)
            ! end if

        end subroutine HistoToFile

        subroutine Tr(dmat)
            implicit none
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
            implicit none
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