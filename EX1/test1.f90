program test1
implicit none

double precision, dimension(:), allocatable ::vect
allocate(vect(10))
do 10 ii=1,10

enddo 

print *, vect !for example $a^b= \int dx\ \phi(x)$
deallocate(vect)

stop
end program test1


