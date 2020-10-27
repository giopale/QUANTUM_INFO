	module hamiltonian

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


	subroutine en_gen_ising(ha1,hb1,lamb1,dd,debug,vec_e,nnvec)
	!The hamiltonian is the following-->
	!H= sum_{alpha}sum_{i} h^{alpha}_{i} +
	!+ sum_{beta}sum_{i} \lambda_{beta}*
	!*h^{beta}_{i}h^{beta}_{i+1}
	
	!This subroutine accepts as imputs the sets
	!of h^{alpha}, h^{beta} and \lambda_{beta}
	!and computes the energies
	!in Mean Field approximation
	!Another imput is the dimension of the state dd
	!for now, we stick to the qubit case (dd=2)

	!From now on, for alpha-part we mean the external field
	!term, and for beta-part the interaction term

	!translational inviance mean field ansatz


	type(dcm) ::ha1
	type(dcm) :: hb1
	type(state) :: psitemp	
	double complex, dimension (:), allocatable :: vectemp
	integer :: dd,nnpar
	integer :: nndis
	logical :: debug
	double precision, dimension (:), allocatable :: vec_e
	integer :: nnvec
	integer :: ii,jj,kk,hh
	integer :: my_stat
	character (256) :: my_msg
	double complex :: contr_a, contr_b
	double precision :: tempcos_ii, tempsin_ii
	double precision :: tempcos_jj, tempsin_jj
	double precision :: valpar_in, valpar_fin_p, valpar_fin_t
	double precision :: step_par_t, step_par_p
	double precision :: lamb1

	nnpar = 2*dd-2
	!number of parameters for a state

	nndis = 100
	!number points discretization for parameters state

	valpar_in=0.0
	valpar_fin_p=2*acos(-1.0)
	valpar_fin_t=acos(-1.0)/2.
	!initial and finale values of parameters
	
	step_par_t=(valpar_fin_t-valpar_in)/(nndis-1)
	step_par_p=(valpar_fin_p-valpar_in)/(nndis-1)
	!step parameters grid

	nnvec = nndis**nnpar

	!total number states takin into account in
	!the minimization process	
	!allocation vector which will contain energies values
	allocate(vec_e(nnvec), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate vec_e with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 


	
	allocate (vectemp(dd), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate vectemp with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!cycle in order to obtaine the values of the energies


	do ii=1, nndis

	do jj=1, nndis

	!initialization of the 1-body state
	!VALID ONLY FOR dd=2 (qubit)
 	call init_state(psitemp, 1 , dd, debug, .true., .true.)
		
	!one-dimensional state, pure and (separable)
	tempcos_ii = cos(valpar_in_t+(ii-1)*step_par_t)
	tempsin_ii = sin(valpar_in_t+(ii-1)*step_par_t)
	tempcos_jj = cos(valpar_in_p+(jj-1)*step_par_p)
	tempsin_jj = sin(valpar_in_p+(jj-1)*step_par_p)

	
	psitemp%coeff(1)=cmplx(tempcos_ii,0.0)

       	psitemp%coeff(2)=cmplx(tempsin_ii*tempcos_jj,
     $			tempsin_ii*tempsin_jj)



	!contribution by the alpha-part 
	vectemp= matmul(ha1%elem,psitemp%coeff)


	contr_a=cmplx(0.0,0.0)

	do hh=1,dd

	contr_a=contr_a + conjg(psitemp%coeff(hh))*vectemp(hh)  
	 
	end do

	!contribution by the beta-part 
	vectemp= matmul(hb1%elem,psitemp%coeff)


	contr_b=cmplx(0.0,0.0)

	do hh=1,dd
	
	contr_b=contr_b + conjg(psitemp%coeff(hh))*vectemp(hh) 

	end do
	
	contr_b= lamb1*(contr_b**2)
	


	!index in the energy vector
	kk=jj-1+(ii-1)*nndis+1

	!value of energy stored	
	if(real(contr_a+contr_b)>0.) then

		vec_e(kk)=abs(contr_a+contr_b)

	else

		vec_e(kk)= - abs(contr_a+contr_b)

	end if


	deallocate(psitemp%coeff)

	end do

	end do

	
	end subroutine en_gen_ising

	



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
		!call random_number(xx)
		!call random_number(yy)
		aa%elem(ii,jj)=cmplx(0.,0.)!complex number
		end do
	end do
	aa%m_trace=(0d0,0d0)
	aa%m_det=(0d0,0d0)
	
	do ii= 1, numrow
		aa%eval(ii)=0.
	end do
		
	end subroutine init_dcm_mat 

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


	subroutine minim(vec,nnvec,min_val,debug)
	
	double precision, dimension (:), allocatable :: vec
	integer :: nnvec
	double precision :: min_val
	logical :: debug
	integer :: ii
	
	min_val=vec(1)

	do ii=1, nnvec
	
		if(min_val>vec(ii)) then
			
			min_val=vec(ii)

		end if

	end do

	end subroutine minim

	end module hamiltonian


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!PROGRAM
	program ex10

	use hamiltonian

	integer:: dd
	logical :: debug
	type(dcm) :: pauli_x
	type(dcm)  :: pauli_z
	double precision, dimension (:), allocatable :: vec_en
	integer :: nnvec
	double precision :: egs
	double precision :: lamb
	integer :: my_stat
	character (256) :: my_msg

	debug=.false.

	!HERE THE PARAMETER lambda IS TAKEN AS IMPUT FROM 
	!FILE "Parameters.txt"

	open(unit = 30, file = "Parameters.txt", 
     $	status = "unknown",
     $  iostat=my_stat, iomsg=my_msg)

	if(my_stat /= 0) then
		print*, 'Open -Parameters- failed with stat = '
     $		         , my_stat, ' msg = '//trim(my_msg)
	end if 
	
	read(30,*) lamb
	close(30)
	
	dd=2
	!dimension psi

	!allocation Pauli matrices
	allocate (pauli_x%elem(2,2), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate pauli_x with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 

	pauli_x%nr=2
	pauli_x%nc=2

	!allocation Pauli matrices
	allocate (pauli_z%elem(2,2), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate pauli_z with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 

	pauli_z%nr=2
	pauli_z%nc=2

	!initialization Pauli matrices
	!along x
	pauli_x%elem(1,1)=(0.,0.)
	pauli_x%elem(1,2)=(1.,0.)
	pauli_x%elem(2,1)=(1.,0.)
	pauli_x%elem(2,2)=(0.,0.)
	!along z
	pauli_z%elem(1,1)=(1.,0.)
	pauli_z%elem(1,2)=(0.,0.)
	pauli_z%elem(2,1)=(0.,0.)
	pauli_z%elem(2,2)=(-1.,0.)	

	call en_gen_ising(pauli_z,pauli_x, lamb, dd,debug,vec_en,nnvec)
	!computation energies

	call minim(vec_en ,nnvec, egs, debug)
	!searching the min

	open(unit = 20, file = "egs.txt", status = "unknown", 
     $		iostat=my_stat, iomsg=my_msg)

		if(my_stat /= 0) then
		print*, "Open failed with stat = "
     $		         , my_stat, " msg = "//trim(my_msg)
		end if

        !print
        write (20,*) egs   
	close(20)
	
	end program ex10
	
		
	
	
	
