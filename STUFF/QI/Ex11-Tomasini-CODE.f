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
	!and one mat_i located in the position pos 
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

	subroutine ten_prod_vec1m(mat_i, mat_f, aa, bb,nn, pos, debug)
	!This subroutine performs the tensor product of 
	!N-1 identies 2x2
	!and one mat_i located in the position pos 
	!of the chain
	!aa and bb number of rows and columns of mat_i
	
	type(dcm) :: mat_i
	type(dcm) :: mat_f
	type(dcm) :: idd
	type(dcm) :: temp_mat_i
	type(dcm) :: tempp_matf
	logical :: debug
	integer :: pos, aa, bb
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
	allocate(temp_mat_i%elem(aa,bb))
	!initialization initial matrix of cycle
	
	if(pos==1) then

		temp_mat_i%elem= mat_i%elem

		temp_mat_i%nr=aa
		temp_mat_i%nc=bb

	else

		temp_mat_i%elem = idd%elem

		temp_mat_i%nr=2
		temp_mat_i%nc=2

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

		if(pos > ii) then

		allocate(temp_mat_i%elem(2**(ii),2**(ii)))

		temp_mat_i%nr=2**(ii)
		temp_mat_i%nc=2**(ii)

		print*, temp_mat_i%nr, ii

		else 

		allocate(temp_mat_i%elem((2**(ii-1))*aa,(2**(ii-1))*bb))

		temp_mat_i%nr=(2**(ii-1))*aa
		temp_mat_i%nc=(2**(ii-1))*bb

		print*, temp_mat_i%nr, ii

		end if


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

	mat_f%nr = (2**(nn-1))*aa
	mat_f%nc = (2**(nn-1))*bb


	end subroutine ten_prod_vec1m


	subroutine ten_prod_vec2(mat_i, mat_f, nn, pos, debug)
	!This subroutine performs the tensor product of 
	!nn-2 identies 2x2
	!and two mat_i located in the position pos 
	!and pos+1 of the chain
	
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
	allocate(temp_mat_i%elem(2,2), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate temp_mat_i%elem with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 


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
	deallocate(temp_mat_i%elem)
	deallocate(idd%elem)

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

	subroutine rowrowmatmul(nn,mm,hh,m1,m2,mris1)
	!first method:row by row
	implicit none
	integer:: ii, jj, kk, nn, mm, hh
	double complex, dimension (:,:), allocatable :: m1
	double complex, dimension (:,:), allocatable :: m2
	double complex, dimension (:,:), allocatable :: mris1	
	do ii=1,nn
	 do jj=1,hh
		mris1(ii,jj)=(0.,0.)
	 end do
	end do	

	do ii=1,nn
	 do jj=1,hh
	  do kk=1,mm 
		mris1(ii,jj)=mris1(ii,jj)+m1(ii,kk)*m2(kk,jj)
		!SCALAR PRODUCT BTW 
		!ii-th ROW OF m1 AND j-th COLUMN OF m2
		!IT GOES IN THE ENTRY (i,j) OF mris1
	  end do
	 end do
	end do

	end subroutine rowrowmatmul

	end module hamiltonian

	

	program week11

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
	double precision :: lamb, truelamb
	double precision, allocatable:: tempeval(:)
	double complex, dimension (: , :), allocatable :: tempmat
	double complex, dimension (: , :), allocatable :: tempmat1

	complex*16, allocatable:: work(:)
	integer lwork
	double precision, allocatable:: rwork(:)
	complex*16, allocatable:: work1(:)
	integer lwork1
	double precision, allocatable:: rwork1(:)
	integer iinfo
	integer :: aa,bb
	type(dcm) :: dham
	type(dcm) :: temp_dmat
	type(dcm) :: temp_dmat1
	type(dcm) :: temp_dmat2
	type(dcm) :: temp_dmat22
	type(dcm) :: idd
	type(dcm) :: half1
	type(dcm) :: half2
	type(dcm) :: temp_int
	double complex, dimension (: , :), allocatable :: pr
	double complex, dimension (: , :), allocatable :: prt
	double complex, dimension (: , :), allocatable :: temp_tr
	integer :: ntr,gg

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

	!number iterations RG algorithm
	
	ntr=100


	!allocation ham
	call init_dcm_mat(ham,2**nn,2**nn,debug)


	!computing the external field part of the hamiltonian	
	do pos=1,nn
	
		!tensor product of N-1 identies and one pauli_z
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


	!duplicating the hamiltonian: from nn to 2nn

	!temporary double matrix
	temp_dmat%nr=2**(2*nn)
	temp_dmat%nc=2**(2*nn)

	temp_dmat1%nr=2**(2*nn)
	temp_dmat1%nc=2**(2*nn)

	!allocation identity 
	allocate (idd%elem(2**(nn),2**(nn)), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate idd with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 

	idd%nr=2**(nn)
	idd%nc=2**(nn)


	!initialization identity
	do ii=1, 2**(nn)

		do jj=1, 2**(nn)

			if(ii==jj) then
				
				idd%elem(ii,ii)=(1.0,0.0)

			else

				idd%elem(ii,jj)=(0.0,0.0)

			end if

		end do

	end do	


	!two halves of the interaction term, first cycle
	call ten_prod_vec1(pauli_x, temp_dmat2, nn, nn, debug)
	call ten_prod_vec1(pauli_x, temp_dmat22, nn, 1, debug)

	!allocation double hamiltonian
	call init_dcm_mat(dham,2**(2*nn),2**(2*nn),debug)

	!necessary to allocate work and rwork
	!and initialize lwork
	!for information read relative documentation
	lwork= 2*(2**(2*nn))
	
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

	!second
	lwork1= 2*(2**(nn))
	
	allocate(work1(lwork1), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate work1 with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 

	allocate(rwork1(lwork1), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate rwork1 with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if






	!temporary matrix
	allocate(temp_tr(2**(2*nn),2**(nn)), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate temp_tr with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 
	
	!projections

	allocate(pr(2**(2*nn),2**(nn)), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate pr with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 

	allocate(prt(2**(nn),2**(2*nn)), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate prt with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if 

	!auxialiary matrices for zheev

	
	allocate(tempeval(2**(2*nn)))

	allocate(tempmat(2**(2*nn),2**(2*nn)))

	allocate(tempmat1(2**(nn),2**(nn)))

	
	!CYCLE OVER TRUNCATIONS

	!in order to keep the numbers low, divide per 2 each cycle
	do ii=1,2**(nn)
	do jj=1,2**(nn)
		ham%elem(ii,jj)=ham%elem(ii,jj)/2.
c		temp_dmat2%elem(ii,jj)=temp_dmat2%elem(ii,jj)/2.
	end do
	end do
	truelamb=lamb
	lamb=lamb/2.

	!BEGINNING CYCLE

	do gg=1, ntr	

	!connecting the two systems			

	!direct product
	
	!first part
	call ten_prod (ham,idd,temp_dmat, debug)

	dham%elem = temp_dmat%elem

	deallocate(temp_dmat%elem)
	deallocate(temp_dmat%eval)


	!second part

	call ten_prod (idd,ham,temp_dmat1, debug)

	dham%elem = dham%elem + temp_dmat1%elem

	deallocate(temp_dmat1%elem)
	deallocate(temp_dmat1%eval)

	!nn interaction among the system

	call ten_prod (temp_dmat2,temp_dmat22,temp_int, debug)
	dham%elem = dham%elem + lamb*temp_int%elem
	
	deallocate(temp_int%elem)
	deallocate(temp_int%eval)

	!calling the LAPACK subroutine for eigen values

	tempmat=dham%elem
	
	aa=2**(2*nn)
	bb=2**(2*nn)


	
	call zheev( "V", "U", aa , tempmat , bb , tempeval , work ,
     $			 lwork, rwork, iinfo )

	!"V"->eigenvalues and eigenvectors, 
	!stored inside tempeval and tempmat
	!"U"->upper triangle of aa is stored, inside dd%elem
	! tempeval will contain the eigen values of aa
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

	!take only first nn eigenvectors and form a projector
	!i.e. the first nn columns of the temporary matrix tempmat

	!projector
	
	do ii=1, 2**(2*nn)

		do jj=1, 2**(nn)

			pr(ii,jj)=tempmat(ii,jj)
			
		end do

	end do


		
	!transpose projector

	do ii=1, 2**(nn)

		do jj=1, 2**(2*nn)

			prt(ii,jj)=conjg(pr(jj,ii))

		end do

	end do


	!computation truncated matrix
	
	
	call rowrowmatmul(2**(2*nn),2**(2*nn),2**(nn),dham%elem,
     $			pr,temp_tr)




	call rowrowmatmul(2**(nn),2**(2*nn),2**(nn),prt,temp_tr,
     $			ham%elem)



	aa=2**(nn)
	bb=2**(nn)

	!in order to keep the numbers low, divide per 2 each cycle
	do ii=1,2**(nn)
	do jj=1,2**(nn)
		tempmat1(ii,jj)=ham%elem(ii,jj)/2.
		ham%elem(ii,jj)=ham%elem(ii,jj)/2.
c		temp_dmat2%elem(ii,jj)=temp_dmat2%elem(ii,jj)/2.
	end do
	end do

	lamb=lamb/2.

	!COMPUTATION EIGENVALUE GS TRUNCATED SYSTEM

	call zheev( "N", "U", aa , tempmat1 , bb , ham%eval , work1,
     $			 lwork1, rwork1, iinfo )

	print*, gg, (2**(gg+1)), ham%eval(1) 

CCCCCCC	!old halves tensor product with identities

	call ten_prod_vec1m(temp_dmat2, half1 , 2**nn, 2**nn, 
     $		nn+1, nn+1, debug)
	call ten_prod_vec1m(temp_dmat22, half2 , 2**nn, 2**nn, 
     $		nn+1, 1, debug)

	

	!computation truncated interaction terms
	
	!first half
	call rowrowmatmul(2**(2*nn),2**(2*nn),2**(nn),half1%elem,
     $			pr,temp_tr)




	call rowrowmatmul(2**(nn),2**(2*nn),2**(nn),prt,temp_tr,
     $			temp_dmat2%elem)

	!second half
	call rowrowmatmul(2**(2*nn),2**(2*nn),2**(nn),half2%elem,
     $			pr,temp_tr)




	call rowrowmatmul(2**(nn),2**(2*nn),2**(nn),prt,temp_tr,
     $			temp_dmat22%elem)

	deallocate(half1%elem)
	deallocate(half2%elem)


	end do

	open(unit = 50, file = "RG_1gs.txt", status = "unknown",
     $		  iostat=my_stat, iomsg=my_msg)

		if(my_stat /= 0) then
		print*, "Open failed with stat = "
     $		         , my_stat, " msg = "//trim(my_msg)
		end if

	write(50,*) truelamb, -abs(ham%eval(1)) 
	

	end program week11
