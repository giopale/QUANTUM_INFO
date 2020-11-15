! External modules:
include "ModDmat.f90"
include "ModDebug.f90"


program Eigenproblem
    use ModDmat
    use ModDebug
    implicit none
    ! Local scalars
    integer :: Nmat, Nsp, ios, ii,jj, Ncycles
    character(len=100) :: filename, msg, x1,x2, format
    ! Local arrays
    double precision, dimension(:), allocatable :: spacings
    double precision :: perc, start, end, time
    ! Local custom types
    type(dmatrix) :: herm_1

    Nmat = 1000
    Ncycles=100

    Nsp=Nmat-1
    herm_1%N = (/ Nmat, Nmat /)
    write(x1,'(i4.4)') Nmat
    write(x2,'(i4.4)') Ncycles
    allocate(spacings(Nsp))
    filename = "Sp_"//trim(x1)//"_"//trim(x2)//".dat"
    open(unit=55, file=filename, status='unknown', iostat=ios, iomsg=msg)
    ! write(55,'("! Matrix size: ", i4.5)') Nmat
    do jj=1,Ncycles
        call cpu_time(start)
        allocate(herm_1%elem(Nmat,Nmat))
        herm_1=herm_1 .Init. Nmat
        herm_1 =.evalh.herm_1
        if (isnan(sum(herm_1%eval))) go to 129
        spacings=SpAvg(herm_1)

        write(x1,'(i4.4)') Nmat
        do ii=1,Nsp
        write(55,'(g13.6)', advance='no') spacings(ii)
        end do
        129 deallocate(herm_1%elem)
        deallocate(herm_1%eval)
        write(55,*)
        call cpu_time(end)
        time=end-start
        perc=100*jj/Ncycles
        write(*,'("Elapsed time [s]: ",(G9.2),(G9.2)," % done...")') time, floor(perc)
    end do

    close(55)


write(*,*) "Che fa Muuu Muuu..." ! closing message





end program Eigenproblem

