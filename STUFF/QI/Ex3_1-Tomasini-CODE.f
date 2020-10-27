	module checkpoint_debug
	implicit none
	
	interface checkpoint!INTERFACE
		module procedure check_none, check_int2, 
     $				 check_real4, check_dcomplex,
     $				 check_r4array
	end interface

	contains

	subroutine check_none(debug,stringg)
	!When no variables are passed
	implicit none
	logical debug
	character(:), allocatable :: stringg
	real ttime

	if(debug.eqv..true.) then
		call cpu_time(ttime)
		print*, " "
		print*, "CHECKPOINT-DEBUGGING"
		print*, stringg
		print*, " "
		print*, "Code time: ", ttime
		print*, " "
	end if
	
	end subroutine check_none

	subroutine check_int2(debug,stringg,var,varcheck)
	!check on an integer*2 variable
	implicit none	
	logical debug
	character(:), allocatable :: stringg
	real ttime
	integer*2 var, varcheck

	if(debug.eqv..true.) then
		call cpu_time(ttime)
		print*, " "
		print*, "CHECKPOINT-DEBUGGING"
		print*, stringg
		print*, " "
		if(var==varcheck) then !checking correctness
			print*, "The value of the variable is:", var
		else 
			print*, "The actual value is:", var
			print*, "It should be", varcheck
		end if
		print*, " "
		print*, "Code time: ", ttime!code time printed
		print*, " "
	end if
	
	end subroutine check_int2

	subroutine check_real4(debug,stringg,var,varcheck)
	!check on an real*4 variable
	implicit none	
	logical debug
	character(:), allocatable :: stringg
	real ttime
	real*4 var, varcheck

	if(debug.eqv..true.) then
		call cpu_time(ttime)
		print*, " "
		print*, "CHECKPOINT-DEBUGGING"
		print*, stringg
		print*, " "
		!checking correctness
		if(abs(var-varcheck)/varcheck < 10E-5 ) then
			print*, "The value of the variable is:", var
		else 
			print*, "The actual value is:", var
			print*, "It should be", varcheck
		print*, " "
		end if
		print*, "Code time: ", ttime!printing code time
		print*, " "
	end if
	
	end subroutine check_real4
	
	subroutine check_dcomplex(debug,stringg,var,varcheck)
	!check on an double complex variable
	implicit none	
	logical debug
	character(:), allocatable :: stringg
	real ttime
	double complex var, varcheck

	if(debug.eqv..true.) then
		call cpu_time(ttime)
		print*, " "
		print*, "CHECKPOINT-DEBUGGING"
		print*, stringg
		print*, " "
		!checking correctness
		if(abs(var-varcheck)/abs(varcheck) < 10E-5 ) then
			print*, "The value of the variable is:", var
		else 
			print*, "The actual value is:", var
			print*, "It should be", varcheck
		print*, " "
		print*, "Code time: ", ttime!printing code time
		print*, " "
		end if
	end if
	
	end subroutine check_dcomplex

	subroutine check_r4array(debug,stringg, m,mcheck,nn,mm)
	implicit none
	!checking if the input arrays are equal
	real*4, dimension(nn,mm) :: m, mcheck
	integer*2 :: tt,ss,nn,mm
	integer*4 accum
	logical debug
	character(:), allocatable :: stringg
	real ttime
	
	if(debug.eqv..true.) then
	call cpu_time(ttime)
	print*, " "
	print*, "CHECKPOINT-DEBUGGING"
	print*, stringg
	print*, " "
	accum=0
	do tt=1,nn
	  do ss=1,mm
           if(abs(m(tt,ss)-mcheck(tt,ss))/mcheck(tt,ss) > 10E-5) then
	  	accum=accum+1
	   end if
	end do
	end do
	!how many different entries?
	if(accum>0) then
		print*, "The two arrays have"
     $		, accum, "different entries"
		print*, " "
	else 
		print*, "Same arrays"
		print*, " "
	end if
	end if

	print*, "Code time: ", ttime
	print*, " "

	end subroutine check_r4array

	end module


	program testing
	use checkpoint_debug	

	implicit none
	logical debug
	character(:), allocatable :: stringg
	integer*2 intvar, intvarcheck
	real*4 realvar, realvarcheck
	double complex dcvar, dcvarcheck
	integer*2 :: ii, jj, nn, mm
	real*4, dimension(:,:), allocatable :: arr,arrcheck
	

	debug=.true.
	nn=4
	mm=5
	allocate(arr(nn,mm))
	allocate(arrcheck(nn,mm))
	do ii=1,nn
		do jj=1,mm
			arr(ii,jj)=ii
			arrcheck(ii,jj)=arr(ii,jj)
		end do
	end do
	dcvar=(3,6)
	dcvarcheck=dcvar

	realvar=5.5555
	realvarcheck=5.556 !bug

	intvar=8
	intvarcheck=intvar

	arr(3,3)=5

	stringg="Checking none"
	call checkpoint(debug,stringg)

	stringg="Checking int2"
	call checkpoint(debug,stringg,intvar,intvarcheck)
	
	stringg="Checking real4"
	call checkpoint(debug,stringg,realvar,realvarcheck)

	stringg="Checking double complex"
	call checkpoint(debug,stringg,dcvar,dcvarcheck)

	stringg="Checking array 2-dim, real4"
	call checkpoint(debug,stringg,arr,arrcheck,nn, mm)

	end program testing


	




			
