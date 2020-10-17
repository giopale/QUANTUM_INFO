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
    contains
        subroutine Tr(dmat)
            ! use type dmatrix
            implicit none
            type(dmatrix) :: dmat
            integer :: m,n,ii
            double complex :: tmp

            m = dmat%N(1)
            n = dmat%N(2)

            if(m/=n) then
                print*, "Error in Tr: matrix dimensions must be equal"
                return
            else 
                do ii=1,n,1
                    tmp=tmp+dmat%elem(ii,ii)
                enddo
                dmat%Trace=tmp
            end if
            return 
        end subroutine Tr

        subroutine InitUni(dmat) 
            ! Initializes a (m,n) complex matrix with 
            ! real and imaginary part taken from [0,1]
            ! uniform distributions
            implicit none
            type(dmatrix) :: dmat
            integer :: jj
            integer, dimension(4) ::ISEED
            double complex, dimension(dmat%N(1),1) :: X

            if(dmat%N(1)<1 .or. dmat%N(2)<1) then
                ! Check for positive matrix shape
                print*, "*** ERROR in InitUni: matrix shape not defined"
                print*, "Program terminated"
                stop
            else
                ISEED = (/ 1, 2, 3 ,4 /)
                do jj=1,dmat%N(2)
                    call clarnv(1, ISEED, 2*dmat%N(1), dmat%elem(:,jj) ) ! '1' stands for uniform
                end do
                return
            end if
        end subroutine InitUni

        subroutine MatToFile(dmat,filename)
            ! Writes dmatrix type to txt file
            implicit none
            type(dmatrix) :: dmat
            character(*) :: filename
            integer :: jj
            filename = filename // ".txt"
            ! open(unit=10,file=filename,action='write',status='unknown')
            do jj=1,dmat%N(2)
                write(*,*) jj!dmat%elem(1,jj)
            end do
            ! close(10)
        end subroutine MatToFile
                

end module stuff

program DMatrixCODE
    use stuff
    implicit none
    type(dmatrix) :: dmat
    integer, dimension(2) ::shape = (/ 6, 2 /)
    integer :: ii,jj

    write(*,*)
    write(*,*) "    *** DMatrixCODE.f03 - Complex matrix init ***  "
    write(*,*)
    print*, "Un buono scrittore deve essere affabulatore, maestro e incantatore."
    print*, "Morte agli zar, viva i Soviet!"
    
    dmat%N=shape
    allocate(dmat%elem(dmat%N(1),dmat%N(2)))
    call InitUni(dmat)

    write(*,*) "Matrix elements:"
    do ii=1,shape(1),1
        write(*,*) dmat%elem(ii,:)
    end do

    call MatToFile(dmat,"matrix")

    ! C'è un errore di accesso in MatToFile: da risolvere
    
end program


