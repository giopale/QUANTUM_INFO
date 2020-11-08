! External modules:
include "ModDmat.f90"
include "ModDebug.f90"


program Eigenproblem
    use ModDmat
    use ModDebug
    implicit none
    ! Local scalars
    character(len=1) :: jobz = "N", uplo = "U"
    integer :: Nmat, lda, lwork, rworkdim, info
    ! Local arrays
    double precision, dimension(:), allocatable :: w
    double complex, dimension(:), allocatable :: work
    double precision, dimension(:), allocatable :: rwork
    ! Local custom types
    type(dmatrix) :: herm_1

    Nmat = 10
    lda = Nmat

    herm_1%N = (/ Nmat, Nmat /)
    allocate(herm_1%elem(herm_1%N(1), herm_1%N(1)))
    herm_1=.Init.herm_1

    herm_1 =.evalh.herm_1

    print*, herm_1%eval

    ! rworkdim = max(1,3*Nmat-2)
    ! allocate(w(Nmat))
    ! allocate(work(1))
    ! allocate(rwork(rworkdim))
    ! call zheev(jobz,uplo,Nmat,herm_1%elem,lda,w,work,-1,rwork,info) ! call to determine LWork
    ! lwork=int(work(1))
    ! deallocate(work)
    ! allocate(work(lwork))
    ! call zheev(jobz,uplo,Nmat,herm_1%elem,lda,w,work,lwork,rwork,info) !call to solve the eigenproblem
    ! allocate(herm_1%eval(Nmat))
    ! herm_1%eval = w
    


    




end program Eigenproblem

