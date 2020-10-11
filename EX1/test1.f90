program test1
implicit none

real :: num=3.13
CHARACTER(len=255) :: path
  CALL getcwd(path)
  WRITE(*,*) TRIM(path)

open(10,file="prova2.txt",status='new',action='write')

write(10,*) num

close(10)

end program test1


