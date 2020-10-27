
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!FIRST EXERCISE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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


	
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!SECOND EXERCISE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




C> \brief \b Multiplication matrix-matrix and its computational time
C>
C>  =========== DOCUMENTATION ===========
C>  Definition:
C> ===========
C>
C>      SUBROUTINE ROWROWMATMUL (nn,m1,m2,mris1,timetot1)
C>
C>       .. Scalar Arguments ..
C>       INTEGER*2            nn
C>	 REAL		   timetot1
C>       ..
C>      .. Array Arguments ..
C>       INTEGER*2          m1(nn,nn)
C>       INTEGER*2          m2(nn,nn)
C>	 INTEGER*2          mris1(nn,nn)
C>      ..
C>
C>  Purpose
C>  =======
C>
C>\details \b Purpose:
C>\verbatim
C>
C> ROWROWMATMUL does a matrix-matrix multiplication between two matrices
C> given as input: m1 and m2. The resulting matrix is mris1. ROWROWMATMUL 
C> deals with INTEGER*2 squared matrices, with numbers of rows and columns equal C> to nn. In the entry (i,j) of mris1 there is the scalar product between the C> i-th row of m1 and the j-th column of m2. 
C> The  multiplication is row by row because while i is fixed j 
C> goes from 1 to nn. In other words the j cycle is inside the i cycle. 
C> Calling the function CPU_TIME at the beginning of that moltiplication, we assign the initial time to REAL variable start1. The same is done for the ending C> time, in REAL variable finish1. The total time timetot1 is the difference between the finish time and the ending time.

C>\endverbatim
C>
C>  Arguments:
C>  ==========
C>

C> \param[in] nn
C> \verbatim
C>          nn is INTEGER*2
C>          The dimension of the squared array m1, m2 and mris1
C> \endverbatim
C> \param[in] m1
C> \verbatim
C>          m1 is INTEGER*2 ARRAY (nn,nn)
C> \endverbatim
C>
C> \param[in] m2
C> \verbatim
C>          m2 is INTEGER*2 ARRAY (nn,nn)
C> \endverbatim
C>
C> \param[out] m2
C> \verbatim
C>          mris1 is INTEGER*2 ARRAY (nn,nn)
C> \endverbatim
C>
C>
C> \param[out] timetot1
C> \verbatim
C>          timetot1 is REAL
C> \endverbatim
C>
C>  Authors:
C> ========
C>
C> \author Univ. of Padua
C>
C> \date 23 October 2018
C>
C>
C>       SUBROUTINE COLCOLMATMUL (nn,m1,m2,mris2,timetot2)
C>
C>       .. Scalar Arguments ..
C>       INTEGER*2            nn
C>	 REAL		   timetot2
C>       ..
C>      .. Array Arguments ..
C>       INTEGER*2          m1(nn,nn)
C>       INTEGER*2          m2(nn,nn)
C>	 INTEGER*2          mris2(nn,nn)
C>       ..
C>
C>  Purpose
C>  =======
C>
C>\details \b Purpose:
C>\verbatim
C>
C> COLCOLMATMUL does a matrix-matrix multiplication between two matrices
C> given as input: m1 and m2. The resulting matrix is mris2. COLCOLMATMUL 
C> deals with INTEGER*2 squared matrices, with numbers of rows and columns equal C> to nn. In the entry (j,i) of mris1 there is the scalar product between the j-th row of m1 and the i-th column of m2. 
C> The  multiplication is row by row because while i is fixed j 
C> goes from 1 to nn. In other words the j cycle is inside the i cycle. 
C> Calling the function CPU_TIME at the beginning of that moltiplication, we assign the initial time to REAL variable start2. The same is done for the ending time, in REAL variable finish2. The total time timetot2 is the difference between the finish time and the ending time.

C>\endverbatim
C>
C>  Arguments:
C>  ==========
C>

C> \param[in] nn
C> \verbatim
C>          nn is INTEGER*2
C>          The dimension of the squared array m1, m2 and mris2
C> \endverbatim
C> \param[in] m1
C> \verbatim
C>          m1 is INTEGER*2 ARRAY (nn,nn)
C> \endverbatim
C>
C> \param[in] m2
C> \verbatim
C>          m2 is INTEGER*2 ARRAY (nn,nn)
C> \endverbatim
C>
C> \param[out] m2
C> \verbatim
C>          mris2 is INTEGER*2 ARRAY (nn,nn)
C> \endverbatim
C>
C>
C> \param[out] timetot2
C> \verbatim
C>          timetot2 is REAL
C> \endverbatim
C>
C>  Authors:
C>  ========
C>
C> \author Univ. of Padua
C>
C> \date 23 October 2018
C>
C>
C>       SUBROUTINE MATMULINTRINSIC (nn,m1,m2,mris1,timetot1)
C>
C>       .. Scalar Arguments ..
C>       INTEGER*2            nn
C>	 REAL		   timetot3
C>       ..
C>       .. Array Arguments ..
C>       INTEGER*2          m1(nn,nn)
C>       INTEGER*2          m2(nn,nn)
C>	 INTEGER*2          mris3(nn,nn)
C>       ..
C>
C>  Purpose
C>  =======
C>
C>\details \b Purpose:
C>\verbatim
C> 
C> MATMULINTRINSIC does a matrix-matrix multiplication between two matrices
C> given as input: m1 and m2. The resulting matrix is mris1. MATMULINTRINSIC 
C> deals with INTEGER*2 squared matrices, with numbers of rows and columns equal C> to nn.  
C> The  multiplication is done via the intrinsic function mris3=matmul(m1,m2).
C> Calling the function CPU_TIME at the beginning of that moltiplication, we  assign the initial time to REAL variable start3. The same is done for the ending C> time, in REAL variable finish3. The total time timetot3 is the difference between the finish time and the ending time.

C>\endverbatim
C>
C> Arguments:
C>  ==========
C>

C> \param[in] nn
C> \verbatim
C>          nn*2 is INTEGER
C>          The dimension of the squared array m1, m2 and mris3
C> \endverbatim
C> \param[in] m1
C> \verbatim
C>          m1 is INTEGER*2 ARRAY (nn,nn)
C> \endverbatim
C>
C> \param[in] m2
C> \verbatim
C>          m2 is INTEGER*2 ARRAY (nn,nn)
C> \endverbatim
C>
C> \param[out] m2
C> \verbatim
C>          mris3 is INTEGER*2 ARRAY (nn,nn)
C> \endverbatim
C>
C>
C> \param[out] timetot1
C> \verbatim
C>          timetot3 is REAL
C> \endverbatim
C>
C> Authors:
C>  ========
C>
C> \author Univ. of Padua
C>
C> \date 23 October 2018
C>
C>       SUBROUTINE CHECKDIM(nn,q,,debug)
C>
C>       .. Scalar Arguments ..
C>       INTEGER*2            nn
C>	 INTEGER*2            qq
C>	 LOGICAL	      debug
C>       ..
C>       .. Array Arguments ..
C>       INTEGER*2          m(nn,nn)
C>       INTEGER*2          mcheck(nn,nn)
C>       ..
C>  Purpose
C>  =======
C>
C>\details \b Purpose:
C>\verbatim
C> 
C> CHECKDIM checks if the dimension of the matrix are correct. Firstly checks
C> if the the dimension INTEGER*2 is above 10000. If it is so, it prints a WARNING
C> message because it takes too much time. Secondly, it checks if nn is inferior  
C> to 1, printing a WARNING message. Lastly, it checks if nn is actually qq*100
C> as it should be. If it is not true, it prints a WARNING message, stating which 
C> should be the dimension and the actual one. If everything goes well, it prints "okay: right dimensions matrices". If something is wrong, DEBUG turns .FALSE. (its value at the beginning should be .TRUE.).

C>\endverbatim
C>
C> Arguments:
C>  ==========
C>

C> \param[in] nn
C> \verbatim
C>          nn is INTEGER*2
C>          The dimension of the squared arrays
C> \endverbatim
C> \param[in] qq
C> \verbatim
C>          qq is INTEGER*2
C>	    the number of the cycle, it should be nn=qq*100
C> \endverbatim
C>
C> \param[inout] debug
C> \verbatim
C>          debug is LOGICAL
C> \endverbatim
C>
C> Authors:
C>  ========
C>
C> \author Univ. of Padua
C>
C> \date 23 October 2018

C>       SUBROUTINE CHECKMAT(nn,m,mcheck,debug)
C>
C>       .. Scalar Arguments ..
C>       INTEGER*2            nn
C>	 LOGICAL	      debug
C>       ..
C>       .. Array Arguments ..
C>       INTEGER*2          m(nn,nn)
C>       INTEGER*2          mcheck(nn,nn)
C>       ..
C>  Purpose
C>
C>  Purpose
C>  =======
C>
C>\details \b Purpose:
C>\verbatim
C> 
C> CHECKMAT checks if the two imput matrices m and mcheck are equal. This is done
C> via a do cycle which checks the equality between the two matrices element
C> by element. For each entry which is not equal, the variable INTEGER*4 accum
C> increseas by 1 (starting from 0). If accum in the end is bigger than 0,
C> a WARNING message states how much entries between the two matrices are 
C> different. Moreover, If something is wrong, DEBUG turns .FALSE. (its value at the beginning should be .TRUE.). Otherwise, it is printed "Same matrices".
C>\endverbatim
C>
C> Arguments:
C>  ==========
C>

C> \param[in] nn
C> \verbatim
C>          nn is INTEGER*2
C>          The dimension of the squared arrays
C> \endverbatim
C>
C> \param[in] 
C> \verbatim
C>          m is INTEGER*2 ARRAY (nn,nn)
C> \endverbatim
C>
C> \param[in] 
C> \verbatim
C>          mcheck is INTEGER*2 ARRAY (nn,nn)
C> \endverbatim
C>
C> \param[inout] debug
C> \verbatim
C>          debug is LOGICAL
C> \endverbatim
C>
C> Authors:
C>  ========
C>
C> \author Univ. of Padua
C>
C> \date 23 October 2018

C>
C>       PROGRAM test_performance_mulmat

C>  Purpose
C>  =======
C>
C>\details \b Purpose:
C>\verbatim
C>
C>This program implements matrix-matrix multiplication in three different method: C>row by row, column by column and by intrinsic function. For each method an C>apposite subroutine is called: ROWROWMATMUL, COLCOLMATMUL and MATMULINTRINSIC. C>For each method the CPU_TIME is computed. This is done for squared matrix, from C>100x100 to 500x500, with a step=100. The first matrix is such that at the entry (i,j) there is the value i+j, whilst the second matrix is such that at the at the entry (i,j) there is the value i*j. The multiplications are done in a do cycle. In each cycle, at the beginning, the matrices are allocated and in the ending are deallocated. The results are printed on a file, called "Results", at unit=40.
C>
C> DEBUGGING of the program. This is controlled with a CHARACTER*1 variable called CHOICE. The programmer can choose to give the choice to do or not the debugging to the user keeping or changing the value of CHOICE. Keeping its default value "X", the user can choose at the beginning of the program to do the debug, via a printed choice on the terminal. 
C> If he says no, the logical variable DEBUG turns .FALSE. and all the debugging routines do not work. 
C> If he says yes, the logical variable DEBUG turns .FALSE. and therefore all the debugging procedures below are executed. If one of the them reports an error, a warning message is printed.
C> Otherwise, if the programmer changes the value of CHOICE in "n", the user has not this choice and the debugging is not done (DEBUG is .FALSE.). Lastly, changing the value of CHOICE in "y", the user has not the choice and the debugging is done (DEBUG is .TRUE.).
C> The debugging is done printing various checkpoints throughout the program, in specific points: at the beginning of each cycle, initialization of the input matrices, calling the subroutines of the three methods, checking the resulting matrices. 
C> Moreover the subroutines CHECKDIM and CHECKMAT are used a lot , interfaced via the operator CHECK. CHECKDIM is used every time a subroutine is called. CHECKMAT checks after the calling of the three method that the imput matrices are still the same (after each one). Moreover, it checks if the resulting matrices mris1, mris and mris3 are the same at the end of each cycle.
C> Another check is printing the istant cpu_time before and after each calling to a function which exploits one method. This is done in order to check immediately if there are some problems with time computation.  
C> Lastly, error handling commands are added for the functions OPEN, WRITE, ALLOCATE and DEALLOCATE. They print which error the functions reports and its relative error message.
C> If one of the these procedures reports an error, a warning message is printed.
C> At the end of the program, the total time spent on debugging and the total time spent are printed. Furthermore, the percentage of time spent on debugging respect to the total time is printed.
 
C>\endverbatim
C>  Authors:
C>  ========
C>
C> \author Univ. of Padua
C>
C> \date 23 October 2018

C>  =====================================================================

	!MODULE MAT-MULTPLICATIONS
	module matmultiplications
	
	contains
	subroutine rowrowmatmul(nn,m1,m2,mris1,timetot1)
	!first method:row by row
	implicit none
	integer*2 ii, jj, kk, nn
	integer*2, dimension(nn,nn) :: m1
	integer*2, dimension(nn,nn) :: m2
	integer*2, dimension(nn,nn):: mris1	
	real :: finish1, start1,timetot1
	call cpu_time (start1)!initial time	
	do ii=1,nn
	 do jj=1,nn
	  do kk=1,nn 
		mris1(ii,jj)=mris1(ii,jj)+m1(ii,kk)*m2(kk,jj)
		!SCALAR PRODUCT BTW 
		!ii-th ROW OF m1 AND j-th COLUMN OF m2
		!IT GOES IN THE ENTRY (i,j) OF mris1
	  end do
	 end do
	end do
        call cpu_time (finish1)!finish time
	timetot1=finish1-start1
	end subroutine rowrowmatmul



	subroutine colcolmatmul(nn,m1,m2,mris2,timetot2)
	!second method:col by col
	implicit none
	integer*2 ii, jj, kk, nn
	integer*2, dimension(nn,nn) :: m1
	integer*2, dimension(nn,nn) :: m2
	integer*2, dimension(nn,nn) :: mris2	
	real :: finish2, start2,timetot2
	call cpu_time (start2)!initial time		
	do ii=1,nn
	 do jj=1,nn
	  do kk=1,nn
		mris2(jj,ii)=mris2(jj,ii)+m1(jj,kk)*m2(kk,ii)
		!SCALAR PRODUCT BTW 
		!jj-th ROW OF m1 AND ii-th COLUMN OF m2
		!IT GOES IN THE ENTRY (j,i) OF mris2
	  end do
	 end do
	end do
	call cpu_time (finish2)!finish time
	timetot2=finish2-start2
	end subroutine colcolmatmul


	subroutine matmulintrinsic(nn,m1,m2,mris3,timetot3)
	!third method: intrisic matmul
	implicit none
	integer*2 ii, jj, kk, nn
	integer*2, dimension(nn,nn) :: m1
	integer*2, dimension(nn,nn) :: m2
	integer*2, dimension(nn,nn) :: mris3	
	real :: finish3, start3,timetot3

	call cpu_time (start3)!initial time
	mris3= matmul(m1,m2)
	!matmul DOES A MATRIX-MATRIX MULTIPLICATION m1*m2
	call cpu_time (finish3)!finish time
	timetot3= finish3-start3
	end subroutine matmulintrinsic

	end module matmultiplications

	!MODULE DEBUGGING
	module debugging

	interface check
		module procedure checkdim,checkmat
	end interface
	
	contains
	subroutine checkdim(nn,qq,debug)!checking dimensions nn
	implicit none
	integer*2 :: nn, qq
	logical :: debug
	if(debug.eqv..TRUE.) then

	if(nn>10000) then !is the dim too large?
			print*, "WARNING: too large dimension:"
			print*, nn
		else if(nn<1) then !is the dim minor than 1?
			print*, "WARNING: dimension minor than 1"
		else if((nn-qq*100)>0.5) then !is the dim wrong?
			print*, "the dimension is wrong, it should be:"
     $				, qq*100
		else
			print*, "okay: right dimensions matrices",nn
		end if
	end if
	end subroutine checkdim

	subroutine checkmat(nn,m,mcheck,debug)
	implicit none
	!checking if the input matrixes are equal
	integer*2, dimension(nn,nn) :: m, mcheck
	integer*2 :: tt,ss,nn
	integer*4 accum
	logical debug

     if(debug.eqv..TRUE.) then
 
	accum=0
	do tt=1,nn
		do ss=1,nn
			if(m(tt,ss)/=mcheck(tt,ss)) then
				accum=accum+1
			end if
		end do
	end do

	if(accum>0) then
		print*, "The two matrices have"
     $			, accum, "different entries"
	else 
		print*, "Same matrices"
	end if

	end if
	end subroutine checkmat

	end module debugging
	

	program test_performance_mulmat
	!This program implements matrix-matrix multiplication
	!in three different method: column by column, row by row
	!and by intrinsic function. For each method the CPU_TIME 
	!is computed. This is done for squared matrix, from 
	!100x100 to 1000x1000, with a step=100.
	!The results are printed on a file, called "Results", 
	!at unit=40. 
	
	!DECLARATION VARIABLES

	use matmultiplications
	use debugging
	implicit none 

	integer*2, dimension(:,:), allocatable :: m1
	integer*2, dimension(:,:), allocatable :: m2
	integer*2, dimension(:,:), allocatable :: mcheck1, mcheck2
	integer*2, dimension(:,:), allocatable :: mris1, mris2, mris3
	integer*2  ii, jj, nn, qq
	real :: ttime
	real :: timetot1, timetot2, timetot3
	real :: t1,t2, accumtime, t3,t4
	logical :: debug
	character*1 :: choice
	integer :: my_stat
	character (256) :: my_msg 
	call cpu_time(t1)!INITIAL TIME ALL PROGRAM
	
	accumtime=0.0
	choice="X" !Initialization of choice variable
		   !if it is "X", it asks the debug
		   !otherwise yoi can choose btw
		   !doing the debug or not writing "y" or "n"
	
	
	!Do you want to debug?
	if(choice=="X") then
		print *, "Do you want to debug?"
		print *, "y for yes, n for no"
		read*, choice	
	end if

	if(choice=="y") then	
		debug=.TRUE. !if sth goes wrong, it becomes false
	else if(choice=="n") then
		debug=.FALSE.
		print*, "Remind: following no okay are not true"
	else
		print*, "Not understood"
		stop
	end if


	!OPENING THE "Result" FILE
	open(unit = 40, file = "Results", 
     $	status = "unknown", access="append",
     $  iostat=my_stat, iomsg=my_msg)

	if(my_stat /= 0) then
		print*, 'Open -Results- failed with stat = '
     $		         , my_stat, ' msg = '//trim(my_msg)
	end if 

 
	!HERE STARTS THE CYCLE OVER THE MATRIX DIMENSION/100
	do qq= 1, 5

		nn=qq*100!MATRIX DIMENSION

		!ALLOCATION MATRICES
		
		!DEBUG
		call cpu_time(t3)
		if(debug.eqv..TRUE.) then!DEBUG
			print*, " "
			print*, "MATRICES OF ORDER ", nn
			print*, " "
		end if	

	
		if(debug.eqv..TRUE.) then!checkdim
			call check(nn,qq,debug)
		end if
		call cpu_time(t4)
		accumtime=t4-t3	

		allocate(m1(nn,nn), stat=my_stat, errmsg=my_msg)
		!allocating m1
		if(my_stat /= 0) then
			print*, 'Failed to allocate m1 with stat = '
     $		         , my_stat, ' and msg = '//trim(my_msg)
		end if 

	        allocate(m2(nn,nn), stat=my_stat, errmsg=my_msg)
		!allocating m2
		if(my_stat /= 0) then
			print*, 'Failed to allocate m2 with stat = '
     $		         , my_stat, ' and msg = '//trim(my_msg)
		end if 

	        allocate(mris1(nn,nn), stat=my_stat, errmsg=my_msg)
		!allocating mris1
		if(my_stat /= 0) then
			print*, 'Failed to allocate mris1 with stat = '
     $		         , my_stat, ' and msg = '//trim(my_msg)
		end if

		allocate(mris2(nn,nn), stat=my_stat, errmsg=my_msg)
		!allocating mris2
		if(my_stat /= 0) then
			print*, 'Failed to allocate mris2 with stat = '
     $		         , my_stat, ' and msg = '//trim(my_msg)
		end if

		allocate(mris3(nn,nn), stat=my_stat, errmsg=my_msg)
		!allocating mris3
		if(my_stat /= 0) then
			print*, 'Failed to allocate mris3 with stat = '
     $		         , my_stat, ' and msg = '//trim(my_msg)
		end if

		allocate(mcheck1(nn,nn), stat=my_stat, errmsg=my_msg)
		!allocating mcheck1
		if(my_stat /= 0) then
			print*, 'Failed to allocate mcheck1 with stat = '
     $		         , my_stat, ' and msg = '//trim(my_msg)
		end if
		allocate(mcheck2(nn,nn), stat=my_stat, errmsg=my_msg)
		!allocating mcheck2
		if(my_stat /= 0) then
			print*, 'Failed to allocate mcheck2 with stat = '
     $		         , my_stat, ' and msg = '//trim(my_msg)
		end if

		!write on file the order of matrix
		write(40,*) " "
		write(40,*) " "
		write(40,*,iostat=my_stat, iomsg=my_msg)
     $          "Matrix squared, order", nn
		write(40,*) " "

		if(my_stat /= 0) then 
			print*, 'Write -Results- failed with stat = '
     $		         , my_stat, ' msg = '//trim(my_msg)
		end if 

		!INITIALIZATION m1: in the entry (i,j)	
		!the value i+j is assigned
        	do ii=1,nn
		 do jj=1,nn
			m1(ii,jj)=ii+jj
			mcheck1(ii,jj)= m1(ii,jj)
		 end do
		end do

		!DEBUG
		call cpu_time(t3)
		if(debug.eqv..TRUE.) then
			print*, " "
			print*, "INITIALIZATION m1"
		end if
		if(debug.eqv..TRUE.) then!checkdim
			call check(nn,qq,debug)
		end if
		call cpu_time(t4)
		accumtime=t4-t3	

		!INITIALIZATION m2, in the entry (i,j)
		!the value i*j is assigned
		do ii=1,nn
		 do jj=1,nn
			m2(ii,jj)=ii*jj
			mcheck2(ii,jj)= m2(ii,jj)
		 end do
		end do
		
		!DEBUG
		call cpu_time(t3)
		if(debug.eqv..TRUE.) then 
			print*, " "
			print*, "INITIALIZATION m2"
		end if
		if(debug.eqv..TRUE.) then !checkdim
			call check(nn,qq,debug)
		end if

		if(debug.eqv..TRUE.) then !check matrices
			print*, "Imput matrix 1"
			call check(nn,m1,mcheck1,debug)
			print*, "Imput matrix 2"
			call check(nn,m2,mcheck2,debug)
		end if
		call cpu_time(t4)
		accumtime=t4-t3	

		!FIRST METHOD: row by row in resulting matrix

		write(40,*) " "		
		write(40,*,iostat=my_stat, iomsg=my_msg)
     $          "	FIRST METHOD: row by row"

		if(my_stat /= 0) then
			print*, 'Write -Results- failed with stat = '
     $		         , my_stat, ' msg = '//trim(my_msg)
		end if 

		call cpu_time(ttime)!printing time istant before 
		print*, ttime
		call rowrowmatmul(nn,m1,m2,mris1,timetot1)
		call cpu_time(ttime)!printing the time istant after
		print*, ttime

                write(40,*,iostat=my_stat, iomsg=my_msg)
     $          	"Time in seconds= ", timetot1

		if(my_stat /= 0) then
			print*, 'Write -Results- failed with stat = '
     $		         , my_stat, ' msg = '//trim(my_msg)
		end if 
		!printing on file the cpu_time
		
		!DEBUG
		call cpu_time(t3)
		if(debug.eqv..TRUE.) then
			print*, " "
			print*, "FIRST METHOD: row by row"
		end if

		if(debug.eqv..TRUE.) then!checkdim
			call check(nn,qq,debug)
		end if

		if(debug.eqv..TRUE.) then!checkmat
			print*, "Imput matrix 1"
			call check(nn,m1,mcheck1,debug)
			print*, "Imput matrix 2"
			call check(nn,m2,mcheck2,debug)
		end if
		call cpu_time(t4)
		accumtime=t4-t3

		!SECOND METHOD: column by column in resulting matrix	
		
		write(40,*) " "
		write(40,*,iostat=my_stat, iomsg=my_msg)
     $          "	SECOND METHOD: col by col"

		if(my_stat /= 0) then
			print*, 'Write -Results- failed with stat = '
     $		         , my_stat, ' msg = '//trim(my_msg)
		end if 


		call cpu_time(ttime)!printing the time istant before
		print*, ttime
		call colcolmatmul(nn,m1,m2,mris2,timetot2)
		call cpu_time(ttime)!printing the time istant after
		print*, ttime


 		write(40,*,iostat=my_stat, iomsg=my_msg)
     $          	"Time in seconds= ", timetot2

		if(my_stat /= 0) then
			print*, 'Write -Results- failed with stat = '
     $		         , my_stat, ' msg = '//trim(my_msg)
		end if 
		!printing on file the cpu_time

		!DEBUG
		call cpu_time(t3)
		if(debug.eqv..TRUE.) then
			print*, " "
			print*, "SECOND METHOD: col by col"
		end if

		if(debug.eqv..TRUE.) then!checkdim
			call check(nn,qq,debug)
		end if
		if(debug.eqv..TRUE.) then!checkmat
			print*, "Imput matrix 1"
			call check(nn,m1,mcheck1,debug)
			print*, "Imput matrix 2"
			call check(nn,m2,mcheck2,debug)
		end if
		call cpu_time(t4)
		accumtime=t4-t3

		!THIRD METHOD: intrinsic function matmul
		
		write(40,*) " "
		write(40,*,iostat=my_stat, iomsg=my_msg)
     $          "	THIRD METHOD: intrinsic matmul"

		if(my_stat /= 0) then
			print*, 'Write -Results- failed with stat = '
     $		         , my_stat, ' msg = '//trim(my_msg)
		end if


		call cpu_time(ttime)!printing the time istant before
		print*, ttime
		call matmulintrinsic(nn,m1,m2,mris3,timetot3)
		call cpu_time(ttime)!printing the time istant after
		print*, ttime


		write(40,*,iostat=my_stat, iomsg=my_msg)
     $          	"Time in seconds= ", timetot3

		if(my_stat /= 0) then
			print*, 'Write -Results- failed with iostat = '
     $		         , my_stat, ' iomsg = '//trim(my_msg)
		end if 
		!printing on file the cpu_time


		!DEBUG
		call cpu_time(t3)
		if(debug.eqv..TRUE.) then
			print*, " "
			print*, "THIRD METHOD: intrinsic matmul"
		end if
		if(debug.eqv..TRUE.) then!checkdim
			call check(nn,qq,debug)
		end if
		if(debug.eqv..TRUE.) then!checkmat
			print*, "Imput matrix 1"
			call check(nn,m1,mcheck1,debug)
			print*, "Imput matrix 2"
			call check(nn,m2,mcheck2,debug)
		end if

		!the resulting matrices are different?
		!DEBUG
		if(debug.eqv..TRUE.) then
			print*, " "
			print*, "CHECKING RESULTING MATRICES"
		end if
		if(debug.eqv..TRUE.) then!checkdim
			print*, "Checking first and second method"
			call check(nn,mris1,mris2,debug)
			print*, "Checking third and second method"
			call check(nn,mris2,mris3,debug)
			print*, "Checking first and third method"
			call check(nn,mris1,mris3,debug)
		end if
		call cpu_time(t4)
		accumtime=t4-t3
		
		!DEALLOCATION matrices:

		deallocate(m1, stat=my_stat, errmsg=my_msg)
		!deallocating m1
		if(my_stat /= 0) then
			print*, 'Failed to deallocate m1 with stat = '
     $		         , my_stat, ' and msg = '//trim(my_msg)
		end if


		deallocate(m2, stat=my_stat, errmsg=my_msg)
		!deallocating m2
		if(my_stat /= 0) then
			print*, 'Failed to deallocate m2 with stat = '
     $		         , my_stat, ' and msg = '//trim(my_msg)
		end if

		deallocate(mris1, stat=my_stat, errmsg=my_msg)
		!deallocating mris1
		if(my_stat /= 0) then
			print*, 'Failed to deallocate mris1 with stat = '
     $		         , my_stat, ' and msg = '//trim(my_msg)
		end if

		deallocate(mris2, stat=my_stat, errmsg=my_msg)
		!deallocating mris2
		if(my_stat /= 0) then
			print*, 'Failed to deallocate mris2 with stat = '
     $		         , my_stat, ' and msg = '//trim(my_msg)
		end if

		deallocate(mris3, stat=my_stat, errmsg=my_msg)
		!deallocating mris3
		if(my_stat /= 0) then
			print*, 'Failed to deallocate mris3 with stat = '
     $		         , my_stat, ' and msg = '//trim(my_msg)
		end if

		deallocate(mcheck1, stat=my_stat, errmsg=my_msg)
		!deallocating mcheck1
		if(my_stat /= 0) then
			print*, 'Failed to deallocate mcheck1 with stat = '
     $		         , my_stat, ' and msg = '//trim(my_msg)
		end if

		deallocate(mcheck2, stat=my_stat, errmsg=my_msg)
		!deallocating mcheck2
		if(my_stat /= 0) then
			print*, 'Failed to deallocate mcheck2 with stat = '
     $		         , my_stat, ' and msg = '//trim(my_msg)
		end if
	end do
	
	call cpu_time(t2)!ENDING TIME ALL PROGRAM

	if(debug.eqv..TRUE.) then
			ttime=(accumtime*100)/(t2-t1)
			print*, " "
			print*, "Time spent on debugging", accumtime
			print*, "Total time spent", t2-t1
			print*, "Percentage", ttime
			print*, " "
	end if

	stop
	end program test_performance_mulmat

	
