	module matrices
	implicit none	
	type dcm !defining the type: doublecomplex matrix
		integer :: nr !number of rows
		integer :: nc !number of columns
		!the actual matrix
		double complex, dimension (: , :), allocatable :: elem
		!eigenvalues
		double precision, dimension (:), allocatable :: eval
		!trace
		double complex :: m_trace
		!determinant
		double complex :: m_det
	end type dcm

	contains

	subroutine init_hermmat (aa, aacheck, numrow, 
     $					numcol,debug, hh, omega) 
	!this subroutine initializes a type dcm, 
	!given as imput. The matrix becomes the one of
	!harmonic oscillator problem
	!Also the numbers of rows and columns 
	!are given as imput. The subroutine creates random matrices
	!with values in range [0,1].
	!It gives a default value 0 to the determinant, the trace
	!and to eigenvalues. 
	real*8 hh, omega
	integer ii,jj 
        type(dcm) :: aa
	type(dcm) :: aacheck
	integer numrow, numcol, nn
	integer :: my_stat
	character (256) :: my_msg 
	logical debug
	
	
	aa%nr= numrow
	aa%nc= numcol
	aacheck%nr= numrow
	aacheck%nc= numcol
	nn= numrow

	!ERROR HANDLING ALLOCATION VECTORS
	allocate (aa%eval(aa%nr), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate aa%eval with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 
	allocate (aacheck%eval(aa%nr), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate aacheck%eval with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 

	!ERROR HANDLING ALLOCATION MATRICES	
	allocate(aa%elem(aa%nr, aa%nc), stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate aa with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 

        allocate(aacheck%elem(aa%nr, aa%nc), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate aacheck with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 


	!WARNING IF NOT SQUARED
	if(debug .eqv. .true.) then
		if(numrow==numcol) then !checking squareness
			print*, " "
			print*, "INITIALIZATION"
			print*, "Okay, squared"
		else 
		        print*, "WARNING: NOT SQUARED"
			print*, "num rows= ", numrow
			print*, "num columns= ", numcol
		end if
	end if

	!INITIALIZATION
	
	!INITIALIZATION
	do ii= 1, aa%nr
		do jj= 1, aa%nc				
		if(ii==jj) then
			aa%elem(ii,jj)=cmplx(2+omega*omega
     $					*(ii-nn/2+1)*(ii-nn/2+1)*hh**4,0)
			aacheck%elem(ii,jj)= aa%elem(ii,jj)
		else if(jj==(ii+1)) then
			aa%elem(ii,jj)=cmplx(-1,0)
			aacheck%elem(ii,jj)= aa%elem(ii,jj)
		end if

		if(jj==(ii-1)) then
			aa%elem(ii,jj)=cmplx(-1,0)
			aacheck%elem(ii,jj)= aa%elem(ii,jj)
		else 
			aacheck%elem(jj,ii)=(0d0,0d0)
			aacheck%elem(jj,ii)= aa%elem(jj,ii)
		end if
		end do
	end do

c	/hh**2

	aa%m_trace=(0d0,0d0)
	aa%m_det=(0d0,0d0)
	aacheck%m_trace=(0d0,0d0)
	aacheck%m_det=(0d0,0d0)
	
	do ii= 1, numrow
		aa%eval(ii)=(0d0,0d0)
		aacheck%eval(ii)=(0d0,0d0)
	end do
		
	end subroutine init_hermmat


	subroutine print_matrices (AA, namefile) 
	! This subroutine prints a type doublecomplex_matrix: 
	!the numbers of row and columns, the elements (in the proper order)
	!, the trace and the determinant. This is all printed on a file
        type(dcm) :: AA
	character(:), allocatable :: namefile
	integer :: ii,jj
	integer :: my_stat
	character (256) :: my_msg 
	      
	open(unit = 40, file = namefile, status = "unknown", 
     $			iostat=my_stat, iomsg=my_msg)

	if(my_stat /= 0) then
		print*, "Open failed with stat = ",
     $		       my_stat, " msg = "//trim(my_msg)
		end if

	write (40,*) "NUMBER OF ROWS:", AA%nr
	write (40,*) "NUMBER OF COLUMNS:", AA%nc
	write (40,*) "ELEMENTS:"
	do ii=1, AA%nr !do cycle in order to print in the proper order
	     write (40,*) (AA%elem(ii, jj), jj = 1, AA%nc)
      	end do	
        write (40,*) "TRACE:", AA%m_trace
	write (40,*) "DET:", AA%m_det
	write (40,*) "EIGENVALUES:"
	do ii=1, AA%nr !do cycle in order to print in the proper order
	       write (40,*) AA%eval(ii)
      	end do
	close(40)
	      
	end subroutine print_matrices

	subroutine print_vector (vec, nn, namefile) 
	! This subroutine prints a vector on file, 
	! given its length nn
     	character(:), allocatable :: namefile
        integer :: ii,nn
	integer :: my_stat
	character (256) :: my_msg 
	double precision, dimension(:), allocatable :: vec
	
	open(unit = 80, file = namefile, status = "unknown", 
     $			iostat=my_stat, iomsg=my_msg)

		if(my_stat /= 0) then
		print*, "Open failed with stat = "
     $		         , my_stat, " msg = "//trim(my_msg)
		end if

        do ii=1, nn !do cycle in order to print
	        write (80,*) vec(ii)
  	end do
	      
	end subroutine print_vector
	
	end module matrices

c/////////////////////////////////////////////////////////

	!MODULE DEBUGGING
	module debugging

	use matrices

	interface check
		module procedure checkdim,checkmat
	end interface
	
	contains
	subroutine checkdim(nn,nncheck,debug)!checking dimensions nn
	implicit none
	integer :: nn, nncheck
	logical :: debug
	
	if(debug.eqv..TRUE.) then

	if(nn>10000) then !is the dim too large?
			print*, "WARNING: too large dimension:"
			print*, nn
		else if(nn<1) then !is the dim minor than 1?
			print*, "WARNING: dimension minor than 1"
		else if((nn-nncheck)>0.5) then !is the dim wrong?
			print*, "the dimension is wrong, it should be:"
     $				, nncheck
			print*, "but is", nn
			print*, " "
		else
			print*, "okay: right dimensions matrices",nn
		end if
	end if
	end subroutine checkdim

	subroutine checkmat(nn,mm,mmcheck,debug)
	implicit none
	!checking if the input matrixes are equal
	type(dcm) :: mm
	type(dcm) :: mmcheck
	integer :: tt,ss,nn
	integer*4 accum
	logical debug

	if(debug.eqv..TRUE.) then
	accum=0
	do tt=1,nn
		do ss=1,nn
			if(abs(mm%elem(tt,ss)-mmcheck%elem(tt,ss))
     $			/abs(mmcheck%elem(tt,ss))> 10E-10) then
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

c////////////////////////////////////////////////////////


	!PROGRAM	
	program harmosc

	use debugging
	use matrices	


	type(dcm) :: aa
	type(dcm) :: aacheck
	type(dcm) :: ddcheck
	type(dcm) :: dd
	
	integer  ii, jj, nn, nncheck
	logical :: debug
	character*1 :: choice
	integer :: my_stat
	character (256) :: my_msg 
	integer, allocatable :: ipiv(:)
	integer :: iinfo
	
	complex*16, allocatable:: work(:)
	integer lwork
	double precision, allocatable:: rwork(:)
	integer iinfodiag
	character(:), allocatable :: diagfile, matfile

	real*8 :: omega, hh
	real*8 :: xmin, xmax
	integer :: kk
	



	choice="n" !Initialization of choice variable
		   !if it is "X", it asks the debug
		   !otherwise the programmer can choose between
		   !doing the debug or not writing "y" or "n"
	
	!Do you want to debug?
	if(choice=="X") then
		print *, "Do you want to debug?"
		print *, "y for yes, n for no"
		read*, choice	
	end if

	if(choice=="y") then	
		debug=.TRUE. 
	else if(choice=="n") then
		debug=.FALSE.
	else
		print*, "Not understood"
		stop
	end if

	!HERE THE MATRIX DIMENSION IS TAKEN AS IMPUT FROM 
	!FILE "MatDimension.txt"
	!THE MATRIX DIMENSION IS THE NUMBER OF POINTS
	!TAKEN INTO ACCOUNT IN THE DISCRETIZATION
	!(the ones where psi is computed)

	open(unit = 30, file = "MatDimension.txt", 
     $	status = "unknown",
     $  iostat=my_stat, iomsg=my_msg)

	if(my_stat /= 0) then
		print*, 'Open -MatDimension- failed with stat = '
     $		         , my_stat, ' msg = '//trim(my_msg)
	end if 
	
	read(30,*) nn	
	nncheck=nn


	!INITIALIZATION aa and aacheck:
	!harm oscillator problem
	!m=0.5
	!hbar=1 
	![natural units]
	!omega^2=0.5*usual omega -> usual omega= 2*omega^2

	xmax= 5.0
	xmin= -xmax
	!interval
	hh= (xmax-xmin)/float(nn-1)
	!discretization interval (nn-1 is spacing`s number)

	omega= 1.0	

	!diagonal entries_(ii)= (2+ omega*x^2*hh^4)/(hh^2)
	!x^2=(ii-nn/2)^2*h^2
	!(ii,ii +/- 1) entries = (-1)/(hh^2)
	!(others)=0
	!the factor 1/h^2 is omitted in order to avoid
	!dividing per zero

	
	call init_hermmat (aa, aacheck, nn, nn,debug, hh, omega)


	if(debug.eqv..TRUE.) then!print matrix
		matfile="Matrix.txt"
		call print_matrices (aa, matfile) 		
	end if

	!DIAGONALIZATION


	!calling the LAPACK subroutine for eigen values
	!necessary to allocate work and rwork
	!and initialize lwork
	!for information read relative documentation
	lwork= 2*nn

	allocate(work(lwork), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate work with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 

	allocate(rwork(lwork), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate rwork with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 
		
	
	call zheev( "V", "U", nn , aa%elem, nn , aa%eval , work , lwork, 
     $			rwork, iinfo1 )
	!"V"->eigenvalues & eigenvectors
	!"U"-> aa%elem contains the orthonormal eigenvectors of the matrix A.
	!aa%eval will contain the eigen values of aa
	!in ascending order(if INFO==0)

	do ii=1, aa%nr !divide per the neglected factor /(hh*hh)
	    	aa%eval(ii)=aa%eval(ii)/(hh*hh)
      	end do	
	
	do ii=1, aa%nr !divide per the neglected factor /(hh*hh)
		do jj=1, aa%nr
	    		aa%elem(ii,jj)=aa%elem(ii,jj)/sqrt(hh)
      		end do
	end do

	!DEBUG
	if(debug.eqv..TRUE.) then
		print*, " "
		print*, "DIAGONALIZATION"
		print*, " "
	end if
	if(debug.eqv..TRUE.) then!checkdim
		call check(nn,nncheck,debug)
	end if


	if(debug.eqv..TRUE.) then!check diag
		if (iinfo1==0) then
			print*, " "
			print*, "Successful DIAGONALIZATION"
		else if (iinfo1 < 0) then
			print*, " "
			print*, "the", iinfo1, "-th argument 
     $				of ipiv had an illegal value"

		else 
			print*, " "			
			print*,  "the algorithm failed to converge"
		end if
	end if 

	if(debug.eqv..TRUE.) then!print diag
		diagfile="Eigen.txt"
		call print_matrices (aa, diagfile) 		
	end if

	!PRINTING ON FILE FIRST kk EIGENVECTORS

	kk=3
	
	open(unit = 20, file = "first_k_eigenvec.txt", 
     $	status = "unknown",
     $  iostat=my_stat, iomsg=my_msg)

	if(my_stat /= 0) then
		print*, 'Open -first_k_eigenvec- failed with stat = '
     $		         , my_stat, ' msg = '//trim(my_msg)
	end if 
	
	do ii=1, aa%nr !do cycle in order to print in the proper order
	     write (20,*) xmin+(ii-1)*hh, 
     $			(real(aa%elem(ii, jj)), jj = 1, kk)
      	end do	

	!PRINTING ON FILE FIRST kk EIGENVALUES

	
	open(unit = 10, file = "first_k_eigenval.txt", 
     $	status = "unknown",
     $  iostat=my_stat, iomsg=my_msg)

	if(my_stat /= 0) then
		print*, 'Open -first_k_eigenval- failed with stat = '
     $		         , my_stat, ' msg = '//trim(my_msg)
	end if 
	
	do ii=1, kk !do cycle in order to print in the proper order
	     write (10,*) aa%eval(ii)
      	end do	


	stop
	end program harmosc	
