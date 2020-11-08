! External modules:
include "ModDmat.f90"
include "ModDebug.f90"


program Eigenproblem
    use ModDmat
    use ModDebug
    implicit none
    ! Local scalars
    integer, parameter :: Nmat = 6
    integer :: seed_dim =4
    integer :: info_zlaghe=0
    ! Local arrays
    integer, dimension(:), allocatable :: seed !Seed for random matrix generator
    double precision, dimension(Nmat) :: Diag_Rnd
    double complex, dimension(2*Nmat) :: Work
    ! Local custom types
    type(dmatrix) :: herm_1

    seed_dim=4  ! Random seed definition 
    call random_seed()
    call random_seed(size=seed_dim)
    allocate(seed(seed_dim))
    call random_seed(get=seed)

    call dlarnv(1,seed,Nmat,Diag_Rnd) ! On exit, the seed is updated
                                      ! Uniform distr. in [0,1]: should I scale it?
    herm_1%N = (/ Nmat, Nmat /)
    allocate(herm_1%elem(Nmat,Nmat))
    ! write(*,*) "Nmat: ", Nmat
    ! write(*,*) "Diag_Rnd: ", Diag_Rnd(:)
    ! write(*,*) "Seed: ", seed
    ! write(*,*) "herm_1%elem: ", herm_1%elem
    ! write(*,*) "info_zlaghe: ", info_zlaghe
    call random_seed()
    call random_seed(get=seed)
    call zlaghe(Nmat, Nmat-1,Diag_Rnd,herm_1%elem, Nmat, seed, work, info_zlaghe)

    ! write(*,*) "herm_1%elem: ", herm_1%elem
    ! write(*,*) "info_zlaghe: ", info_zlaghe




    




end program Eigenproblem

