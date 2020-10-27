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

	subroutine init_hermmat (aa, aacheck, numrow, numcol,debug) 
	!this subroutine initializes a type dcm, 
	!given as imput. The matrix becomes Herminian.
	!Also the numbers of rows and columns 
	!are given as imput. The subroutine creates random matrices
	!with values in range [0,1].
	!It gives a default value 0 to the determinant, the trace
	!and to eigenvalues. 
	real*8 yy, xx
	integer ii,jj 
        type(dcm) :: aa
	type(dcm) :: aacheck
	integer numrow, numcol
	integer :: my_stat
	character (256) :: my_msg 
	logical debug
	
	
	aa%nr= numrow
	aa%nc= numcol
	aacheck%nr= numrow
	aacheck%nc= numcol

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
	do ii= 1, aa%nr
		do jj= ii, aa%nc				
		call random_number(xx)
		call random_number(yy)
		aa%elem(ii,jj)=cmplx(xx,yy)!complex number
		aa%elem(jj,ii)=conjg(aa%elem(ii,jj))!hermitianity
		!in order to be hermitian, complex part on diagonale
		!must be null
		if(ii==jj) then
			aa%elem(ii,jj)=cmplx(xx,0)
		end if	
		aacheck%elem(ii,jj)= aa%elem(ii,jj)
		aacheck%elem(jj,ii)= aa%elem(jj,ii)
		end do
	end do
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
	module normalized_spacings
	
	use matrices

	contains

	subroutine normspac(aa,kk,normsp,debug)
	!computes normalized spacings between eigenvalues
	!of a double complex matrix, which is inside type aa
	type(dcm) :: aa
	integer :: kk
	!kk-1 is the number of the spacings to be considered
	!in the normalization computation
	integer :: ii
	integer :: iimax,iimin
	double precision, dimension(:), allocatable :: sp
	double precision, dimension(:), allocatable :: norm
	double precision, dimension(:), allocatable :: normsp
	integer :: my_stat
	character (256) :: my_msg 
	character(:), allocatable :: namefile
	logical debug
	

	allocate(sp(aa%nr-1), stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate normsp with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 

	allocate(norm(aa%nr-1), stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate normsp with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 

	!spacings
	do ii=1, aa%nr-2
		sp(ii)=aa%eval(ii+1)-aa%eval(ii)
	end do
	
	if(debug.eqv..TRUE.) then!check after changing
		namefile="SPACINGS.txt"
		call print_vector(sp, aa%nr-1,namefile)
	end if

	!NORMALIZATION and NORMALIZED SPACINGS
	
	!if cycle to be sure to have (kk-1) in [1,nn-1]
	if((kk-1)>(aa%nr-1)) then
		
		print*, " Too many spaces (kk-1)in the normalization"
		print*, " Rescaled to", (aa%nr-1)
		print*, " It was", kk
		kk=aa%nr
	else if((kk-1)<1) then
		print*, " Not enough spaces (kk-1) in the normalization"
		print*, " Rescaled to 2"
		print*, " It was", kk
		kk=2
	end if
	
	!if cycle to distinguish between odd case and even case	
	do ii=1,aa%nr-1
		!setting the borders of the interval of normalization
		if(mod(kk,2)==1) then
			iimax=ii+(kk-1)/2
			iimin=ii-(kk-1)/2			
		else 
			iimax=ii-1+(kk)/2 
			iimin=ii-(kk)/2 
		end if
		!boundary conditions
		
		if(debug.eqv..TRUE.) then!check before changing
			print*, "BEFORE corrections"
			print*, "Spacing number", ii
			print*, "right border:", iimax
			print*, "left border:", iimin
			print*, " "
		end if
		if(iimax>aa%nr) then!exceeding at right
			iimax=aa%nr
			iimin=aa%nr-(kk-1)
		else if(iimin<1) then!exceeding at left
			iimin=1
			iimax=kk
		end if

		if(debug.eqv..TRUE.) then!check after changing
			print*, "AFTER corrections"
			print*, "Spacing number:", ii
			print*, "right border:", iimax
			print*, "left border:", iimin
			print*, " "
		end if

		!normalization of ii-spacing
		norm(ii)=(aa%eval(iimax)-aa%eval(iimin))/(kk-1)

		if(debug.eqv..TRUE.) then!check after changing
			namefile="NORMS.txt"
			call print_vector(norm, aa%nr-1,namefile)
		end if		

		!normalized ii-spacing
		normsp(ii)=sp(ii)/norm(ii)		
		
	end do
	end subroutine normspac
	
	

			                                                                       
	

	end module normalized_spacings

	

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
	
	

	!PROGRAM	
	program eigenproblem

	use debugging
	use matrices	
	use normalized_spacings

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
	character(:), allocatable :: lufile, diagfile, normfile

	double precision, dimension(:), allocatable :: normsp
	integer kk
	logical diagherm
	type(dcm) :: diag
	type(dcm) :: diagcheck
	real*8 :: temp
	double precision, dimension(:), allocatable :: normspdiag


	diagherm = .true. !if it is .true. then the procedure is done 
			   !only on Hermitian matrices. 
			   !If .false on Diagonal.

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
	!AND THE PARAMETER REGARDING NORMALIZED SPACINGS
	!BETWEEN EIGENVALUES

	open(unit = 30, file = "MatDimension.txt", 
     $	status = "unknown",
     $  iostat=my_stat, iomsg=my_msg)

	if(my_stat /= 0) then
		print*, 'Open -MatDimension- failed with stat = '
     $		         , my_stat, ' msg = '//trim(my_msg)
	end if 
	
	read(30,*) nn, kk	
	nncheck=nn

	if(diagherm.eqv..true.) then

	!INITIALIZATION aa and aacheck:
	!random numbers in range [0,1], flat distribution
	!hermitian matrices
	
	call init_hermmat (aa, aacheck, nn, nn,debug)

	!DEBUG
	if(debug.eqv..TRUE.) then
		print*, " "
		print*, "INITIALIZATION DEBUG"
	end if
	if(debug.eqv..TRUE.) then!checkdim
		call check(nn,nncheck,debug)
	end if

	!INITIALIZATION dd:
	!random numbers in range [0,1]
	!hermitian matrices
	!later they will be overwritten
	!dd used for diagonalization
	!ddcheck useless
	call init_hermmat (dd, ddcheck, nn, nn, debug)



	!DIAGONALIZATION

	!the matrix dd%elem is equal to aa%elem	
	do ii=1,nn
		do jj=1,nn
			dd%elem(ii,jj)=aa%elem(ii,jj)
		end do
	end do	

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
		
	
	call zheev( "N", "U", nn , dd%elem, nn , aa%eval , work , lwork, 
     $			rwork, iinfo1 )
	!"N"->only eigenvalues 
	!"U"->upper triangle of aa is stored, inside dd%elem
	! aa%eval will contain the eigen values of aa
	!in ascending order(if INFO==0)
	
	!DEBUG
	if(debug.eqv..TRUE.) then
		print*, " "
		print*, "DIAGONALIZATION"
		print*, " "
	end if
	if(debug.eqv..TRUE.) then!checkdim
		call check(nn,nncheck,debug)
	end if

	if(debug.eqv..TRUE.) then!checkmat
		print*, " "
		print*,"aa and aacheck"
		call check(nn,aa,aacheck,debug)
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
		diagfile="Eigenvalues.txt"
		call print_matrices (aa, diagfile) 		
	end if

	!NORMALIZED SPACINGS BETWEEN EIGENVALUES
	
	!allocations
	
	allocate(normsp(aa%nr-1), stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate normsp with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 


	!calling the subroutine normspac to compute normalized spacings
	!IMPUT: the type aa, the parameter kk (kk-1 spacings
	!considered in the normalization) taken as imput from file
	!OUTPUT: the vector of normalized spaces, nn-1 entries
	
	call normspac(aa,kk,normsp,debug)
	
	normfile="NormSpacingsHerm.txt"

	call print_vector(normsp,nn-1,normfile)

	!*******************************************************
	else
	!*******************************************************

	!COMPUTING NORMALIZED SPACING for diagonal real matrices

	!INITIALIZATION diag:
	!using the same function used for dcm type
	!just to use the same subroutines
	!then we RE-INTIALIZE to diagonal real matrices
	!same dimension as hermitian ones	
	call init_hermmat (diag, diagcheck, nn, nn, debug)

	do ii=1,nn
		do jj=1,nn
			if(ii==jj) then
			  	call random_number(temp)
			  	diag%elem(ii,ii)=complex(temp, 0d0)
				diagcheck%elem(ii,ii)=diag%elem(ii,ii)
				diag%eval(ii)= temp
			else
			  	diag%elem(ii,jj)=complex(0d0, 0d0)
				diagcheck%elem(ii,jj)=diag%elem(ii,jj)
			end if
		end do
	end do	 

	!SORTING OF EIGENVALUES
	!via the LAPACK subroutine hpsort_eps_epw  
	!"I" for incresing order
	!"D" for decreasing
	
	call dlasrt ("I",diag%nr, diag%eval, iinfo)

	if(debug.eqv..TRUE.) then!check diag
		if (iinfo1==0) then
			print*, " "
			print*, "Successful SORTING"
		else if (iinfo1 < 0) then
			print*, " "
			print*, "the", iinfo1, "-th argument 
     $				had an illegal value"

		end if
	end if 

	if(debug.eqv..TRUE.) then!print diag
		diagfile="Eigenvalues.txt"
		call print_matrices (diag, diagfile) 		
	end if

	!NORMALIZED SPACINGS
	!allocations
	
	allocate(normspdiag(diag%nr-1), stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate normspdiag with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 


	!calling the subroutine normspac to compute normalized spacings
	!IMPUT: the type diag, the parameter kk (kk-1 spacings
	!considered in the normalization) taken as imput from file
	!OUTPUT: the vector of normalized spaces, nn-1 entries
	
	call normspac(diag,kk,normspdiag,debug)
	
	normfile="NormSpacingsDiag.txt"

	call print_vector(normspdiag,nn-1,normfile)

	end if

	stop
	end program eigenproblem









