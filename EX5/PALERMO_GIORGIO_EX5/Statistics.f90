! External modules:
include "ModDmat.f90"
include "ModDebug.f90"


program Statistics
    use ModDmat
    use ModDebug
    implicit none
    ! Local scalars
    character(len=100) :: msg, x1, x2, filename, filename_1
    integer :: Nsp, Ncy, ios, Nbins, ii, seme_dim=4
    double precision :: h_lenght, h_step, h_lower, h_upper
    ! Local arrays
    double precision, dimension(:,:), allocatable :: indata
    double precision, dimension(:), allocatable :: diag_rnd, sp_diag
    double precision, dimension(:), allocatable :: h_bounds, h_input_aux
    integer, dimension(:), allocatable :: h_input, h, seme
    
    ! Local custom types
    type(histogram) :: h1, h2
    write(*,*) "As I'm getting older, chip up on my shoulder..."

    write(x1,'(i4.4)') Nsp+1
    write(x2,'(i4.4)') Ncy
    filename = "Hist_"//trim(x1)//"_"//trim(x2)//".dat"
    filename_1 = "Hist_"//trim(x1)//"_"//trim(x2)//"_diag_mat"//".dat"


    Nsp=1999
    Ncy=450

    allocate(indata(Ncy,Nsp),diag_rnd(Nsp+1),sp_diag(Nsp))
    indata = LoadArray("Sp_2000_0450.dat",Ncy,Nsp)

    h_lower = 0
    h_upper = 8
    Nbins=50

    h1=InitHisto(h_lower,h_upper,Nbins)
    h2=InitHisto(h_lower,h_upper,Nbins)
    do ii=1,Ncy
        call FillH(h1,indata(ii,:),"y")
    end do

    call random_seed()
    call random_seed(size=seme_dim)
    allocate(seme(seme_dim))
    call random_seed(get=seme)

    call random_number(diag_rnd)
    call dlasrt('I',Nsp+1,diag_rnd,ios)
    print*, ios

    sp_diag=SpAvgDble(diag_rnd)


    call FillH(h2,sp_diag,'y')


    call HistoToFile(h1,filename)
    call HistoToFile(h2,filename_1)

    write(*,*) "Salad days are gone, missing hippy John..."
end program Statistics