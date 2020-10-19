! **********************************************************************
! Project           : Quantum information, Ex2
! 
! Program name      : DMatrixCODE.f03
! 
! Author            : Giorgio Palermo
! 
! Date created      : 20201014
! 
! Purpose           : Define and manipulate complex matrices
! 
! Revision History  :
!
! Date        Author      Ref    Revision (Date in YYYYMMDD format) 
!
! **********************************************************************

module stuff

    type dmatrix
        integer, dimension(2) :: N = (/ 0,0 /)
        double complex, dimension(:,:), allocatable :: elem
        double complex :: Trace
        double complex :: Det
    end type dmatrix

    interface operator (.Adj.)
            module procedure Adj
    end interface 

    interface operator (.Init.)
        module procedure InitUni
    end interface

    contains
        
        function InitUni(dmat) 
            ! Initializes a (m,n) complex matrix with 
            ! real and imaginary part taken from [0,1]
            ! uniform distributions
            implicit none
            type(dmatrix), intent(in) :: dmat
            type(dmatrix) :: InitUni
            integer :: jj,sd=4
            integer, dimension(:), allocatable :: seed
            double complex, dimension(dmat%N(1),1) :: X

            if(dmat%N(1)<1 .or. dmat%N(2)<1) then
                ! Check for positive matrix shape
                print*, "*** ERROR in InitUni: matrix shape not defined"
                print*, "Program terminated"
                stop
            else
                call random_seed(size = sd)
                allocate(seed(sd))
                call random_seed(get=seed)
                ! seed = (/ 1, 3, 3 ,1 /)
                do jj=1,dmat%N(2)
                    call clarnv(1, seed, 2*dmat%N(1), dmat%elem(:,jj) ) ! '1' stands for uniform
                end do
                InitUni = dmat
                return
            end if
        end function InitUni

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

        subroutine SAdj(dmatIN,dmatOUT)
            !computes the adjoint and copies it on another 
            ! dmatrix type element
            type(dmatrix) :: dmatIN, dmatOUT
            dmatOUT%Trace=conjg(dmatIN%Trace)
            dmatOUT%Det=conjg(dmatIN%Det)
            dmatOUT%N=(/ dmatIN%N(2), dmatIN%N(1) /)
            allocate(dmatOUT%elem(dmatOUT%N(1),dmatOUT%N(2)))
            dmatOUT%elem=dmatIN%elem
            dmatOUT%elem=conjg(dmatOUT%elem)
            dmatOUT%elem=transpose(dmatOUT%elem)
            return
        end subroutine SAdj

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
                print*, "*** ERROR in Tr: matrix dimensions must be equal"
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

                

end module stuff

program DMatrixCODE
    use stuff
    implicit none
    type(dmatrix) :: dmat, dmatOUT, dmatOUT1
    integer, dimension(2) ::shape = (/ 2, 2 /)
    integer :: ii,jj
    double complex ::suca, suca1

    write(*,*)
    write(*,*) "    *** DMatrixCODE.f03 - Complex matrix manipulation ***  "
    write(*,*)
    print*, "Un buono scrittore deve essere affabulatore, maestro e incantatore."
    print*, "Morte agli zar, viva i Soviet!"
    
    dmat%N=shape
    allocate(dmat%elem(dmat%N(1),dmat%N(2)))
    allocate(dmatOUT1%elem(dmat%N(1),dmat%N(2)))

    dmat = .init.dmat
    call SAdj(dmat,dmatOUT)
    call Tr(dmat)
    dmatOUT1 = .Adj.dmat

    call MatToFile(dmat,"matrix")
    call MatToFile(dmatOUT1,"matrix_conj")
    
    deallocate(dmat%elem,dmatOUT%elem)

end program


