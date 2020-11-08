! External modules:
include "ModDmat.f90"
include "ModDebug.f90"


program Statistics
    use ModDmat
    use ModDebug
    implicit none
    ! Local scalars
    character(len=100) :: msg
    integer :: Nsp, Ncy, ios, Nbins, ii
    double precision :: h_lenght, h_step, h_lower, h_upper
    ! Local arrays
    double precision, dimension(:,:), allocatable :: indata
    double precision, dimension(:), allocatable :: h_bounds, h_input_aux
    integer, dimension(:), allocatable :: h_input, h
    
    ! Local custom types
    type(histogram) :: h1
    write(*,*) "As I'm getting older, chip up on my shoulder..."


    Nsp=999
    Ncy=100

    allocate(indata(Ncy,Nsp))
    indata = LoadArray("Sp_1000_0100.dat",Ncy,Nsp)
    ! write(*,'(999G15.5)') suca1

    h_lower = 0
    h_upper = 10
    Nbins=50

    
    h1=InitHisto(h_lower,h_upper,Nbins)
   
    do ii=1,Ncy
        h1=FillHisto(h1,indata(ii,:))
    end do

    h1%entries=sum(h1%h)
    print*, h1%entries
    print*, h1%h(:)

    ! Problema: per Ncy=13 non funziona print*, h1%h


    write(*,*) "Salad days are gone, missing hippy Jon..."
end program Statistics