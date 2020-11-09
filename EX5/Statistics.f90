! External modules:
include "ModDmat.f90"
include "ModDebug.f90"


program Statistics
    use ModDmat
    use ModDebug
    implicit none
    ! Local scalars
    character(len=100) :: msg, x1, x2, filename
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
    h_upper = 14
    Nbins=50

    
    h1=InitHisto(h_lower,h_upper,Nbins)
   
    do ii=1,Ncy
        h1=FillHisto(h1,indata(ii,:))
    end do

    write(x1,'(i4.4)') Nsp+1
    write(x2,'(i4.4)') Ncy
    filename = "Hist_"//trim(x1)//"_"//trim(x2)//".dat"
    open(unit=55, file=filename, status='unknown', iostat=ios, iomsg=msg)
    do ii=1,h1%Nbins,1
        write(55,'(2G15.5)') h1%bincenters(ii), h1%h(ii)
    end do
    close(55)



    write(*,*) "Salad days are gone, missing hippy Jon..."
end program Statistics