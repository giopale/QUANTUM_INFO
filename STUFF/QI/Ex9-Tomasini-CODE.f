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


	contains

	subroutine ten_prod (m1,m2, mr, debug)
	!This subroutine performs a tensor product
	!between two matrices m1 (at left) and m2
	!(at right), giving as a result the matrix mr

	type(dcm) :: m1
	type(dcm) :: m2
	type(dcm) :: mr
	logical :: debug
	integer :: ii, jj, kk, hh

	call init_dcm_mat (mr, m1%nr*m2%nr, m1%nc*m2%nc,debug)

	do ii=1, m2%nr

		do jj=1, m2%nc
			
			do kk=1, m1%nr

			do hh=1, m1%nc

			ind_r=kk+(ii-1)*m1%nr
			ind_c=hh+(jj-1)*m1%nc
			
			mr%elem(ind_r,ind_c)= 
     $				m2%elem(ii,jj)*m1%elem(kk,hh)


			end do
			
			end do


		end do

	end do



	end subroutine ten_prod


	subroutine ten_prod_vec1(mat_i, mat_f, nn, pos, debug)
	!This subroutine performs the tensor product of 
	!N-1 identies 2x2
	!and one pauli_z located in the position pos 
	!of the chain
	
	type(dcm) :: mat_i
	type(dcm) :: mat_f
	type(dcm) :: idd
	type(dcm) :: temp_mat_i
	type(dcm) :: tempp_matf
	logical :: debug
	integer :: pos
	integer :: ii, jj
	integer :: my_stat
	character (256) :: my_msg

	!allocation iddentity
	allocate (idd%elem(2,2), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate idd with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 

	idd%nr=2
	idd%nc=2


	!initialization iddentity 2x2
	idd%elem(1,1)=(1.,0.)
	idd%elem(1,2)=(0.,0.)
	idd%elem(2,1)=(0.,0.)
	idd%elem(2,2)=(1.,0.)

	
		
	!allocate initial matrix of cycle
	allocate(temp_mat_i%elem(2,2))
	!initialization initial matrix of cycle
	temp_mat_i%nr=2
	temp_mat_i%nc=2
	
	if(pos==1) then

		temp_mat_i%elem= mat_i%elem

	else

		temp_mat_i%elem = idd%elem

	end if


	!cycle of tensor product
	do ii=2, nn

		!if pos==ii then the tensor product with the 
		!initial matrix is performed
		!otherwise with the identity

		if(pos==ii) then
						
			call ten_prod(temp_mat_i,mat_i,tempp_matf, debug)

		else

			call ten_prod(temp_mat_i,idd,tempp_matf, debug)

		end if

		deallocate(temp_mat_i%elem)	

		allocate(temp_mat_i%elem(2**(ii),2**(ii)))

		temp_mat_i%nr=2**(ii)
		temp_mat_i%nc=2**(ii)

		temp_mat_i%elem=tempp_matf%elem
		deallocate(tempp_matf%elem)
		deallocate(tempp_matf%eval)

	end do
	
	allocate(mat_f%elem(2**nn,2**nn), stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate mat_f with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 	

	!final matrix 2**n x 2**n
	mat_f%elem = temp_mat_i%elem

	mat_f%nr = 2**nn
	mat_f%nc = 2**nn


	end subroutine ten_prod_vec1


	subroutine ten_prod_vec2(mat_i, mat_f, nn, pos, debug)
	!This subroutine performs the tensor product of 
	!N-1 identies 2x2
	!and one pauli_z located in the position pos 
	!of the chain
	
	type(dcm) :: mat_i
	type(dcm) :: mat_f
	type(dcm) :: idd
	type(dcm) :: temp_mat_i
	type(dcm) :: tempp_matf
	logical :: debug
	integer :: pos
	integer :: ii, jj
	integer :: my_stat
	character (256) :: my_msg

	!allocation iddentity
	allocate (idd%elem(2,2), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate idd with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 

	idd%nr=2
	idd%nc=2


	!initialization iddentity 2x2
	idd%elem(1,1)=(1.,0.)
	idd%elem(1,2)=(0.,0.)
	idd%elem(2,1)=(0.,0.)
	idd%elem(2,2)=(1.,0.)

	
		
	!allocate initial matrix of cycle
	allocate(temp_mat_i%elem(2,2))
	!initialization initial matrix of cycle
	temp_mat_i%nr=2
	temp_mat_i%nc=2
	
	if(pos==1) then

		temp_mat_i%elem = mat_i%elem

	else

		temp_mat_i%elem = idd%elem

	end if


	!cycle of tensor product
	do ii=2, nn

		!if pos==ii then the tensor product with the 
		!initial matrix is performed
		!otherwise with the identity

		if(pos==ii) then
						
			call ten_prod(temp_mat_i,mat_i,tempp_matf, debug)

		else if ((pos+1)==ii) then

			call ten_prod(temp_mat_i,mat_i,tempp_matf, debug)

		else

			call ten_prod(temp_mat_i,idd,tempp_matf, debug)

		end if

		deallocate(temp_mat_i%elem)	

		allocate(temp_mat_i%elem(2**(ii),2**(ii)))

		temp_mat_i%nr=2**(ii)
		temp_mat_i%nc=2**(ii)

		temp_mat_i%elem=tempp_matf%elem
		deallocate(tempp_matf%elem)
		deallocate(tempp_matf%eval)

	end do
	
	allocate(mat_f%elem(2**nn,2**nn), stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate mat_f with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 	

	!final matrix 2**n x 2**n
	mat_f%elem = temp_mat_i%elem

	mat_f%nr = 2**nn
	mat_f%nc = 2**nn


	end subroutine ten_prod_vec2



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
	double precision, dimension(:), allocatable :: vec
	
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


	end module hamiltonian

	

	program week9

	use hamiltonian
	

	logical :: debug
	integer :: ii, jj, kk, hh
	character(:), allocatable :: namefile
	type(dcm) :: pauli_x
	type(dcm)  :: pauli_z
	type(dcm)  :: mat_f
	type(dcm)  :: ham
	integer :: my_stat
	character (256) :: my_msg
	integer :: nn, pos
	type(dcm) :: m1,m2,mr
	double precision :: lamb
	double precision, allocatable:: tempeval(:)
	double complex, dimension (: , :), allocatable :: tempmat

	complex*16, allocatable:: work(:)
	integer lwork
	double precision, allocatable:: rwork(:)
	integer iinfo
	integer :: aa,bb

c	call init_dcm_mat(m1,3,3,debug)
c	call init_dcm_mat(m2,3,3,debug)


	debug=.false.

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

	!HERE THE NUMBER SITES nn IS TAKEN AS IMPUT FROM 
	!FILE "MatDimension.txt"
	!AND THE PARAMETER lambda

	open(unit = 30, file = "Parameters.txt", 
     $	status = "unknown",
     $  iostat=my_stat, iomsg=my_msg)

	if(my_stat /= 0) then
		print*, 'Open -MatDimension- failed with stat = '
     $		         , my_stat, ' msg = '//trim(my_msg)
	end if 
	
	read(30,*) nn, lamb


	!allocation ham
	call init_dcm_mat(ham,2**nn,2**nn,debug)
	

	!computing the external field part of the hamiltonian	
	do pos=1,nn
	
		!tensor product of N-1 identies and one pauli_z
		!located in the position pos of the chain
		call ten_prod_vec1(pauli_z,mat_f,nn,pos,debug)

		ham%elem=ham%elem+mat_f%elem
		deallocate(mat_f%elem)

	end do

	
	!computing the sites interactions part of the hamiltonian	
	do pos=1,nn-1
	
		!tensor product of N-2 identies and 
		!two closepauli_x, the first is located in
		!the position pos of the chain
		call ten_prod_vec2(pauli_x,mat_f,nn,pos,debug)

		ham%elem=ham%elem+lamb*mat_f%elem
		deallocate(mat_f%elem)

	end do


c	do ii=1,ham%nr
c	do jj=1,ham%nc
c
c	print*, ham%elem(ii,jj)
c
c	end do
c	end do

	
c	print*, "vv",ham%elem(ham%nr,ham%nc)
c	print*, ham%nr, ham%nc

	!calling the LAPACK subroutine for eigen values
	!necessary to allocate work and rwork
	!and initialize lwork
	!for information read relative documentation
	lwork= 2*(2**nn)

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
		
	allocate(tempeval(2**nn))

	allocate(tempmat(2**nn,2**nn))

	tempmat=ham%elem
	
	aa=2**nn
	bb=2**nn

	call zheev( "N", "U", aa , tempmat , bb , tempeval , work , lwork, 
     $			rwork, iinfo )

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


	if(debug.eqv..TRUE.) then!check diag
		if (iinfo==0) then
			print*, " "
			print*, "Successful DIAGONALIZATION"
		else if (iinfo < 0) then
			print*, " "
			print*, "the", iinfo, "-th argument 
     $				of ipiv had an illegal value"

		else 
			print*, " "			
			print*,  "the algorithm failed to converge"
		end if
	end if 

	namefile="ham_eval.txt"
	call print_vector(tempeval, 2**nn, namefile, 10)



	
	end program week9
