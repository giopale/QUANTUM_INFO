	module wavefunctions
	
	implicit none
		
	type state 
	!defining the type: 
	!separable pure state for nn subsystems
	!which wavefunction belongs to H^D, i.e.
	!i.e. D dimensional Hilbert space
	
		integer :: nsub 
		!number of subsystems

		integer :: ndim 
		!dimension Hilbert space single psi

		!the coefficients of the total wave function
		double complex, dimension (:), allocatable :: coeff

		logical :: sep
		!if .true. it is separable, if .false., no

		logical :: pur
		!if .true. it is pure, if .false., no

		integer :: len_state
		!length of the state

	end type state 


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

	subroutine init_state (psi, nn, dd, debug, sepp, purr) 
	!this subroutine initializes a type state 
	!with random real numbers between 0 and 1
	!If sepp.eqv .true. then the state is separable
	!hence a nn*dd vector
	!the first coefficient of each wavefunction
	!of one subsystem is set to real
	!each single subsystem wavefunction
	!is normalized to one
	!If sepp.eqv .false. then the state is not separable
	!hence a dd**n vector
	!the first coefficient of the wavefunction
	!of one subsystem is set to real
	!the wavefunction is normalized to one
 
	integer :: ii, jj, kk
	integer :: nn, dd
        type(state) :: psi
	integer :: my_stat
	character (256) :: my_msg 
	logical :: debug, sepp, purr
	double precision :: xx, yy
	double precision :: summ, tempp
	character(:), allocatable :: namefile
	
	
	psi%nsub= nn
	psi%ndim= dd

	psi%sep=sepp
	

	!SEPARABLE
	if(psi%sep .eqv. .true.) then

	psi%len_state=nn*dd

	!WARNING IF NOT POSITIVE
	if(debug.eqv..TRUE.) then
		if(psi%nsub <= 0) then
			print*, "nn <= 0"
			print*, "The actual value is", nn
		end if
	
		if(psi%ndim <= 0) then
			print*, "dd <= 0"
			print*, "The actual value is", dd
		end if
	end if	


	!ALLOCATION VECTOR
	allocate (psi%coeff(nn*dd), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate psi%coeff with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 

	!INITIALIZATION COEFFICIENTS
	!subsystem by subsystem
	do ii= 1, psi%nsub
		
		do jj= 1, psi%ndim
		  call random_number(xx)
		
		  if (jj==1) then
			psi%coeff(jj+(ii-1)*psi%ndim)=cmplx(xx,0)	
			!setting the global phase
			!setting the first coeff real
		  else
			call random_number(yy)
			psi%coeff(jj+(ii-1)*psi%ndim)=cmplx(xx,yy)
		  end if	
		end do
		
		summ=0.0
		!norm computation for each subsystem
		do kk= 1, psi%ndim

		  tempp= abs(psi%coeff(kk+(ii-1)*psi%ndim))**2
		  summ=summ + tempp
		end do
		summ=sqrt(summ)
	
		!normalizing the coefficients 
		!dedicated to each subsystem
		!in such a way that each subsystem is normalized
		!to 1/sqrt(number of subsystems)
		!Hence the total wave function is normalized to 1
		do kk= 1, psi%ndim

		  psi%coeff(kk+(ii-1)*psi%ndim)= 
     $			psi%coeff(kk+(ii-1)*psi%ndim)/(summ)
		end do

		if(debug.eqv..TRUE.) then

		!check on norm
	  	print*, " "
		print*, "NORM of subsystem number: ", ii
		summ=0.0
		!norm computation 
		do kk= 1, psi%ndim*psi%nsub
		 tempp= abs(psi%coeff(kk+(ii-1)*psi%ndim))**2
		 summ=summ + tempp
		end do
		summ=sqrt(summ)
		print*, summ
		print*, " "

		end if
	end do

	if(debug.eqv..TRUE.) then!check
	  namefile="sep_purestate.txt"
	  call print_vector(psi%coeff, nn*dd, namefile, 80)
	  print*, "all subsystems"
	  do kk= 1, psi%ndim*psi%nsub
		 print*, psi%coeff(kk)
	  end do
	  print*, " "

	end if
	
	!------------	
	else !GENERAL
	!------------

	psi%len_state=dd**nn

	!ALLOCATION VECTOR
	allocate (psi%coeff(dd**nn), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate psi%coeff with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 

	!INITIALIZATION COEFFICIENTS
	do ii= 1, psi%ndim**psi%nsub

		  call random_number(xx)
		
		  if (ii==1) then
			psi%coeff(ii)=cmplx(xx,0)	
			!setting the global phase
			!setting the first coeff real
		  else
			call random_number(yy)
			psi%coeff(ii)=cmplx(xx,yy)
		  end if	
	end do
		
	summ=0.0
	!norm computation 
	do kk= 1, psi%ndim**psi%nsub
		  tempp= abs(psi%coeff(kk))**2
		  summ=summ + tempp
	end do
	summ=sqrt(summ)

	!normalizing the coefficients 
	do kk= 1, psi%ndim**psi%nsub
		  psi%coeff(kk)= psi%coeff(kk)/summ
	end do

	if(debug.eqv..TRUE.) then!check
	  namefile="gen_purestate.txt"
	  call print_vector(psi%coeff, dd**nn, namefile, 30)

	  !check on norm
	  print*, "GENERAL PURE STATE"
	  print*, "NORM of general pure state"
	  summ=0.0
	  !norm computation 
	  do kk= 1, psi%ndim**psi%nsub
		  tempp= abs(psi%coeff(kk))**2
		  summ=summ + tempp
	  end do
	  summ=sqrt(summ)
	  print*, summ
	  print*, " "

	  do kk= 1, psi%ndim*psi%nsub
		 print*, psi%coeff(kk)
	  end do


	end if

	end if

	end subroutine init_state

	


	subroutine denmat_pure(psi, denmat, debug)
	!this subroutine computes the
	!density matrix of a general pure state
	!of N subystems, each one belonging
	!to a D-dimensional Hilbert space
	!hence the density matrix is D**N x D**N

	type(state) :: psi
	type(state) :: psi_large
	type(dcm) :: denmat
	integer :: ii,jj, ll, tt , ss
	character(:), allocatable :: namefile
	logical :: debug, sepp, purr
	integer:: nr_denmat, nc_denmat
	integer :: accum

	nr_denmat= psi%len_state
	nc_denmat= psi%len_state

	call init_dcm_mat (denmat, nr_denmat, nc_denmat,debug) 

	!-------------------------
	if(psi%sep .eqv. .true.) then
	!if the psi_is separable, written like a nn*dd vector
	!this subroutine gives the relative dd**nn vector 
	!in the multi-state basis
	call fromsep_togen(psi,psi_large)
	!-------------------------

	!COMPUTATION DENSITY MATRIX SEPARABLE STATE
	do ii= 1, psi_large%len_state
	
		do jj= 1, psi_large%len_state

		denmat%elem(ii,jj)= 
     $		  conjg(psi_large%coeff(ii))*psi_large%coeff(jj)

		end do
	end do

	if(debug.eqv..TRUE.) then!check

		print*, " "
		print*, "DENSITY MATRIX"

		!computation of TRACE
		!if the psi is a pure state
		!then it should be one
		print*, " "

		do ii= 1, psi_large%len_state
		 denmat%m_trace= 
     $			denmat%m_trace +denmat%elem(ii,ii)
		end do

		print*, "TRACE: ", denmat%m_trace
		print*, " "

		!check of hermitianity
		accum=0
		do tt=1,psi%len_state
			do ss=1,psi%len_state

			if(abs(conjg(denmat%elem(tt,ss))-denmat%elem(ss,tt))
     $			/abs(denmat%elem(ss,tt))> 10E-10) then
				accum=accum+1
			end if

			end do
		end do

		if(accum>0) then
			print*, "Hermitianity violated in"
     $			, accum, "different entries"
		else 
			print*, "Hermitianity checked"
		end if
		namefile="sep_denmat.txt"
		call print_matrices(denmat, namefile, 60)
	
	end if

	!-------------------------
	else
	!-------------------------
	

	!COMPUTATION DENSITY MATRIX GENERAL STATE
	do ii= 1, psi%len_state
	
		do jj= 1, psi%len_state

		denmat%elem(ii,jj)= 
     $		  conjg(psi%coeff(ii))*psi%coeff(jj)

		end do
	end do

	if(debug.eqv..TRUE.) then!check
		
		print*, " "
		print*, "DENSITY MATRIX"
		!computation of TRACE
		!if the psi is a pure state
		!then it should be one
		print*, " "

		do ii= 1, psi%len_state
		 denmat%m_trace= 
     $			denmat%m_trace +denmat%elem(ii,ii)
		end do

		print*, "TRACE: ", denmat%m_trace
		print*, " "
		!check of hermitianity
		accum=0
		do tt=1,psi%len_state
			do ss=1,psi%len_state

			if(abs(conjg(denmat%elem(tt,ss))-denmat%elem(ss,tt))
     $			/abs(denmat%elem(ss,tt))> 10E-8) then
				accum=accum+1
			print*, accum
			end if

			end do
		end do

		if(accum>0) then
			print*, "Hermitianity violated in"
     $			, accum, "different entries"
		else 
			print*, "Hermitianity checked"
		end if

		namefile="gen_denmat.txt"
		call print_matrices(denmat, namefile, 50)
	
	end if

	end if

	end subroutine denmat_pure

	subroutine fromsep_togen(psi,psi_large)
	!this subroutine takes a separable state
	!which information is stored in nn*dd vector
	!and gives the dd**nn vector in the multi-state basis

	type(state) :: psi
	type(state) :: psi_large
	logical sepp, purr, debug
	integer :: ii, jj, ll
	integer :: cc, kk
	integer, dimension(:), allocatable :: vec_alpha
	double precision :: summ, tempp
	double complex :: temp

	sepp=.false.
	purr=.true.

	debug=.false.

	call init_state(psi_large, psi%nsub, psi%ndim, debug, sepp, purr)

	allocate(vec_alpha(psi%nsub))
		


	!this cycle allows to find the coefficients
	!of the multistate basis as a function of the 
	!coefficients of the one-state basis
	do ii= 1, psi_large%len_state

		kk=ii-1
		!this cycle finds which indexes (alpha) of the vector nn*dd
		!are associate to an index of the vector dd**nn
		!they are stored in a vector vec_alpha
		do jj= 1, psi%nsub

		cc=psi%nsub-jj
	
		vec_alpha(cc+1)=int((kk-mod(kk,psi%ndim**cc))
     $					/psi%ndim**cc +0.1) +1		

		kk=kk-(vec_alpha(cc+1)-1)*(psi%ndim**cc)


		end do

		temp=(1.,0.)

		do jj= 1, psi%nsub

		temp= temp*
     $				psi%coeff(vec_alpha(jj)+(jj-1)*psi%ndim)

		end do
		psi_large%coeff(ii)=temp

	end do

	if(debug.eqv..TRUE.) then!check

		print*, "Enlarged psi"
		do ii= 1, psi_large%len_state
		  print*, psi_large%coeff(ii)
		end do

		!norm computation 
		summ=0.
	  	do kk= 1, psi_large%len_state
		  tempp= abs(psi_large%coeff(kk))**2
		  summ=summ + tempp
	 	 end do
		 print*, "NORM enlarged psi"
	 	 summ=sqrt(summ)
	 	 print*, summ
	 	 print*, " "



	end if

	end subroutine fromsep_togen


	subroutine red_denmat(denmat_tot, 
     $			denmat_red, dd, nn, kk, debug)
	!this subroutine computes the reduced
	!density matrix of a general density matrix
	!in a dd**nn-dimensional Hilbert space
	!The trace is performed on the kk-th subsystem
	!Hence the reduced density matrix
	!is in a dd**(nn-1)-dimensional Hilbert space

	type(dcm) :: denmat_tot
	type(dcm) :: denmat_red
	integer :: ii, jj, kk
	integer :: aa ,bb
	character(:), allocatable :: namefile
	integer :: nr_red, nc_red
	logical :: debug 
	integer :: nn, dd
	integer :: ind_c, ind_r, ind_kk
	integer :: tt, ss, accum

	nr_red= dd**(nn-1)
	nc_red= dd**(nn-1)

	
	call init_dcm_mat (denmat_red, nr_red, nc_red,debug) 


	!COMPUTATION REDUCED DENSITY MATRIX
	do aa= 1, nr_red

		
		do bb= 1, nc_red
		
		denmat_red%elem(aa,bb)=(0.0, 0.0)	
		 do ind_kk= 1, dd
		  
            ind_r = mod(aa-1,(dd**(kk-1))) + (ind_kk-1)*(dd**(kk-1))
     $	     +(aa-1-mod(aa-1,(dd**(kk-1))))*(dd)+1 

		

            ind_c = mod(bb-1,(dd**(kk-1))) + (ind_kk-1)*(dd**(kk-1))
     $	     +(bb-1-mod(bb-1,(dd**(kk-1))))*(dd)+1		 
		  denmat_red%elem(aa,bb)= denmat_red%elem(aa,bb)+
     $	  	        denmat_tot%elem(ind_r,ind_c)


		 end do

		end do

	end do

	if(debug.eqv..TRUE.) then!check
		
		print*, " "
		print*, "DENSITY MATRIX"
		!computation of TRACE
		!if the psi is a pure state
		!then it should be one
		print*, " "

		do ii= 1, dd**nn
		 denmat_red%m_trace= 
     $			denmat_red%m_trace +denmat_red%elem(ii,ii)
		end do

		print*, "TRACE: ", denmat_red%m_trace
		print*, " "
		!check of hermitianity
		accum=0
		do tt=1,dd**(nn-1)
			do ss=1,dd**(nn-1)
	if(abs(conjg(denmat_red%elem(tt,ss))-denmat_red%elem(ss,tt))
     $			/abs(denmat_red%elem(ss,tt))> 10E-8) then
				accum=accum+1
			print*, accum
			end if

			end do
		end do

		if(accum>0) then
			print*, "Hermitianity violated in"
     $			, accum, "different entries"
		else 
			print*, "Hermitianity checked"
		end if
	
	end if

	if(debug.eqv..TRUE.) then!check
		namefile="denmat_red.txt"
		call print_matrices(denmat_red, namefile, 90)
	end if
	
	end subroutine red_denmat




	subroutine print_matrices (AA, namefile, unitt) 
	! This subroutine prints a type doublecomplex_matrix: 
	!the numbers of row and columns, the elements (in the proper order)
	!, the trace and the determinant. This is all printed on a file
        type(dcm) :: AA
	character(:), allocatable :: namefile
	integer :: ii,jj, unitt
	integer :: my_stat
	character (256) :: my_msg 
	      
	open(unit = unitt, file = namefile, status = "unknown", 
     $			iostat=my_stat, iomsg=my_msg)

	if(my_stat /= 0) then
		print*, "Open failed with stat = ",
     $		       my_stat, " msg = "//trim(my_msg)
		end if

	write (unitt,*) "NUMBER OF ROWS:", AA%nr
	write (unitt,*) "NUMBER OF COLUMNS:", AA%nc
	write (unitt,*) "ELEMENTS:"
	do ii=1, AA%nr !do cycle in order to print in the proper order
	     write (unitt,*) (AA%elem(ii, jj), jj = 1, AA%nc)
      	end do	
        write (unitt,*) "TRACE:", AA%m_trace
	write (unitt,*) "DET:", AA%m_det
	write (unitt,*) "EIGENVALUES:"
	do ii=1, AA%nr !do cycle in order to print in the proper order
	       write (unitt,*) AA%eval(ii)
      	end do
	close(unitt)
	      
	end subroutine print_matrices

	subroutine print_vector (vec, nn, namefile, unitt) 
	! This subroutine prints a vector on file, 
	! given its length nn
     	character(:), allocatable :: namefile
        integer :: ii,nn, unitt
	integer :: my_stat
	character (256) :: my_msg 
	double complex, dimension(:), allocatable :: vec
	
	open(unit = unitt, file = namefile, status = "unknown", 
     $			iostat=my_stat, iomsg=my_msg)

		if(my_stat /= 0) then
		print*, "Open failed with stat = "
     $		         , my_stat, " msg = "//trim(my_msg)
		end if

        do ii=1, nn !do cycle in order to print
	        write (unitt,*) vec(ii)
  	end do
	      
	close(unitt)
	end subroutine print_vector
	

	subroutine init_dcm_mat (aa, numrow, numcol,debug) 
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
	integer numrow, numcol
	integer :: my_stat
	character (256) :: my_msg 
	logical debug
	
	
	aa%nr= numrow
	aa%nc= numcol

	!ERROR HANDLING ALLOCATION VECTORS
	allocate (aa%eval(aa%nr), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate aa%eval with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 

	!ERROR HANDLING ALLOCATION MATRICES	
	allocate(aa%elem(aa%nr, aa%nc), stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate aa with stat = '
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
		do jj= 1, aa%nc				
		call random_number(xx)
		call random_number(yy)
		aa%elem(ii,jj)=cmplx(xx,yy)!complex number
		end do
	end do
	aa%m_trace=(0d0,0d0)
	aa%m_det=(0d0,0d0)
	
	do ii= 1, numrow
		aa%eval(ii)=(0d0,0d0)
	end do
		
	end subroutine init_dcm_mat 

	end module wavefunctions


c       ////////////////////////////////////////////////////
C	PROGRAM

	program test_wf

	use wavefunctions


	logical :: debug, sepp, purr
	integer :: nn, dd
	type(state) :: psi_sep
	type(state) :: psi_gen
	type(dcm) :: denmat_tot
	type(dcm) :: denmat_red

	debug=.True.

	nn=2
	!number of subsystems

	dd=2
	!number of dimension of the Hilbert space
	!of single subsystem

	sepp=.true.
	purr=.true.

c	call init_state (psi_sep, nn, dd, debug, sepp, purr) 
	!initializing separable pure state

	sepp=.false.
	purr=.true.

	call init_state (psi_gen, nn, dd, debug, sepp, purr)
	!initializing general pure state

c	call denmat_pure(psi_sep, denmat_tot, debug)
	!computing density matrix from separable pure state

	call denmat_pure(psi_gen, denmat_tot, debug)
	!computing density matrix from general pure state

	kk=2

	call red_denmat(denmat_tot, denmat_red, 
     $					dd, nn, kk, debug)

	!computing reduced density matrix from general pure state

	stop
	end program test_wf






	
