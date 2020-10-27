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

	subroutine init_mat(aa, aacheck, numrow, numcol,debug, xmin,
     $ 			 xmax, omega, hbar, mm) 
	!this subroutine initializes a type dcm, 
	!given as imput. The matrix becomes the one of
	!harmonic oscillator problem
	!Also the numbers of rows and columns 
	!are given as imput
	!It gives a default value 0 to the determinant, the trace
	!and to eigenvalues.  
	real*8 hh
	integer ii,jj 
        type(dcm) :: aa
	type(dcm) :: aacheck
	integer numrow, numcol, nn
	integer :: my_stat
	character (256) :: my_msg 
	logical debug
	real*8 :: omega, hbar, mm
	real*8 :: xmax, xmin

	
	aa%nr= numrow
	aa%nc= numcol
	aacheck%nr= numrow
	aacheck%nc= numcol
	nn= numrow

	!CHECKS

	if(debug.eqv..TRUE.) then
	if(mm <= 0) then
	print*, "mass <= 0"
	print*, "The actual value is", mm
	end if
	end if	
	
	if(debug.eqv..TRUE.) then
	if(hbar <= 0) then
	print*, "hbar <= 0"
	print*, "The actual value is",hbar
	end if
	end if	
	
	if(debug.eqv..TRUE.) then
	if(omega <= 0) then
		print*, "omega <= 0"
		print*, "The actual value is", omega
	end if
	end if	

	!symmetric interval
		
	if(debug.eqv..TRUE.) then
	if(xmax <= 0) then
		print*, "x_max <= 0"
		print*, "The actual value is", xmax
	end if
	if(xmin >= 0) then
		print*, "x_min => 0"
		print*, "The actual value is", xmin
	end if
	end if	
	

	hh= (xmax-xmin)/float(nn-1)
	!discretization interval (nn-1 is spacing`s number)
		

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
	!diagonal entries_(ii)= 
	!(2*(hbar*hbar)/(2*m)+ omega*omega*x^2*hh^2)/(hh^2)
	!with x^2=(xmin+(ii-1)*hh)**2
	!(ii,ii +/- 1) entries = ((hbar*hbar)/(2*m))*(-1)/(hh^2)
	!(others)=0
	!the factor 1/h^2 is omitted in order to avoid
	!dividing per zero

	do ii= 1, aa%nr
		do jj= 1, aa%nc				
		if(ii==jj) then
			aa%elem(ii,jj)=cmplx(2*((hbar*hbar)/(2*mm))+
     $			omega*omega*(xmin+(ii-1)*hh)*(xmin+(ii-1)*hh)*hh**2,0)
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
		
	end subroutine init_mat


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

	subroutine print_vector (vec, nn, namefile, unitt) 
	! This subroutine prints a vector on file, 
	! given its length nn
     	character(:), allocatable :: namefile
        integer :: ii,nn, unitt
	integer :: my_stat
	character (256) :: my_msg 
	real, dimension(:), allocatable :: vec
	
	open(unit = unitt, file = namefile, status = "unknown", 
     $			iostat=my_stat, iomsg=my_msg)

		if(my_stat /= 0) then
		print*, "Open failed with stat = "
     $		         , my_stat, " msg = "//trim(my_msg)
		end if

        do ii=1, nn !do cycle in order to print
	        write (unitt,*) vec(ii)
  	end do
	      
	end subroutine print_vector

	function pot(xmin,hh,ii,ti,hhtt,jj,tf) result(vv)
	!giving in imput lower boundaries of space interva;
	!and time interval, their steps of discretization
	!and the indexes of the considered point
	!the potential in that point ai that time is given
	integer :: ii,jj
	real*8 :: xmin, ti, hh, hhtt, tf
	real*8 :: vv
	
	vv= 0.5* ((xmin+ii*hh)-(ti+jj*hhtt)/(tf-ti))**2
	
	end function pot

	
	subroutine tfourier(nn, vec,signn, nnkk, tfvec,xmin,xmax,hh)
	!this subroutine performs a fourier transform
	!of vec if sign=-1, an anti fourier transform
	!if sign=1

	double complex, allocatable:: vec(:)
	double complex, allocatable:: tfvec(:)	
	integer :: signn
	integer :: nn, nnkk
	integer :: ii, jj
	double precision :: twopi
	double precision :: xmin,xmax, hh

	twopi = 2*acos(-1.0)
	do ii=-nnkk/2, nnkk/2
		tfvec(ii+nnkk/2+1)=(0.,0.)
		do jj=1, nn
	tfvec(ii+nnkk/2+1)=tfvec(ii+nnkk/2+1) + vec(jj)*cexp(cmplx(0.,
     $	signn*twopi*((xmin+float(jj-1)*hh)/(xmax-xmin))*float(ii)))	
		end do
		tfvec(ii+nnkk/2+1)=tfvec(ii+nnkk/2+1)
     $		/sqrt(float(nn))
	end do
	end subroutine tfourier
	
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
	program harmosctime

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
	real*8 :: hbar, mm

	integer :: nntt, nnttcheck
	real*8:: ti, tf, hhtt

	real*8 :: vv
	double complex, dimension (:,:), allocatable :: psi
	double complex, allocatable:: tfpsi(:)
	double complex, allocatable:: psitemp(:)
	integer :: signn
	double precision :: p, two pi, xxx, nnn, sig
	integer :: cc, bb
	double complex :: temp

	integer :: ist
	character(len=100) :: istlabel, tempch
	character(len=:), allocatable :: filename
	character(len=1000) :: filenametemp
	
	twopi=2*acos(-1.0)


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
	mm=1.0
	hbar=1.0 
	![natural units]

	if(debug.eqv..TRUE.) then
	if(mm <= 0) then
	print*, "mass <= 0"
	print*, "The actual value is", mm
	end if
	end if	
	
	if(debug.eqv..TRUE.) then
	if(hbar <= 0) then
	print*, "hbar <= 0"
	print*, "The actual value is",hbar
	end if
	end if	
	

	xmax= 5.0
	xmin= -xmax
	!symmetric interval, order 10^0
		
	if(debug.eqv..TRUE.) then
	if(xmax <= 0) then
		print*, "x_max <= 0"
		print*, "The actual value is", xmax
	end if
	if(xmin >= 0) then
		print*, "x_min => 0"
		print*, "The actual value is", xmin
	end if
	end if	

	hh= (xmax-xmin)/float(nn-1)
	!discretization interval (nn-1 is spacing`s number)

	omega=1.0/sqrt(2.0)

	!omega

	if(debug.eqv..TRUE.) then
	if(omega <= 0) then
		print*, "omega <= 0"
		print*, "The actual value is", omega
	end if
	end if		

	!diagonal entries_(ii)= 
	!(2*(hbar*hbar)/(2*m)+ omega*omega*x^2*hh^2)/(hh^2)
	!with x^2=(xmin+(ii-1)*hh)**2
	!(ii,ii +/- 1) entries = ((hbar*hbar)/(2*m))*(-1)/(hh^2)
	!(others)=0
	!the factor 1/h^2 is omitted in order to avoid
	!dividing per zero

	
	call init_mat (aa, aacheck, nn, nn,debug, xmin, xmax, omega,
     $							 hbar, mm)



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
	
	do ii=1, aa%nr !divide per the normalization factor
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


c	///////////////////////////////////////////////////

	!TIME EVOLUTION of the first eigenvector,
	! ACCORDING TO THE HAMILTONIAN:
	!H(t)=0.5*p^2 + 0.5*(q-q(t))^2
	!q(t)=t/(tf-ti)
	!hence m=1, hbar=1 and omega=1/sqrt(2)	

	!FAST FOURIER TRANSFORM METHOD 

	!time interval[ti,tf]
	ti=0. 
	tf=5
	
	!HERE THE NUMBER OF ISTANTS IS TAKEN AS IMPUT FROM 
	!FILE "NumberIstants.txt"
	!DISCRETIZATION OF TIME

	open(unit = 50, file = "NumberIstants.txt", 
     $	status = "unknown",
     $  iostat=my_stat, iomsg=my_msg)

	if(my_stat /= 0) then
		print*, 'Open -NumberIstants- failed with stat = '
     $		         , my_stat, ' msg = '//trim(my_msg)
	end if 
	
	read(50,*) nntt	
	nnttcheck=nntt

	hhtt=(tf-ti)/float(nntt-1)
	!discretization interval (nn-1 is spacing`s number)

	!initialize a double complex matrix
	!with each column equal to starting psi
	!number of rows=nn
	!number of col=nntt
	allocate(psi(nn,nntt), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate psi with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if

	do cc=1, nntt  
		do bb=1, nn
	    	 	psi(bb,cc)=aa%elem(bb,1)
		end do
      	end do
		



	!initialize two double complex vectors 
	!useful for fourier
	!length=# of spatial points
	allocate(tfpsi(nn), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate psi with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if
	

	!initialize a double complex vector useful for fourier
	!length=# of temporal points
	allocate(psitemp(nn), 
     $			stat=my_stat, errmsg=my_msg)
	!error handling allocation	
	if(my_stat /= 0) then
		print*, 'Failed to allocate psi with stat = '
     $	         , my_stat, ' and msg = '//trim(my_msg)
	end if
	
	!cycle over the istants
	do jj=1, nntt-1
		!multiplication per factor exp(-0.5*img*V_ll*hhtt) 
		!with V_ll= V(xmin+hh*(ll-1)) 
		do ll=1, nn 
		  vv=pot(xmin,hh,ll-1,ti,hhtt,jj, tf)
	    	  psi(ll,jj+1)=psi(ll,jj)*cexp(cmplx(0., -0.5*vv*hhtt))
		  temp = psi(ll,jj+1)
		  psitemp(ll)= temp
      		end do	
		
		signn=-1!fourier transfom

		call tfourier(nn, psitemp ,signn, 
     $			nn, tfpsi, xmin, xmax, hh)
	
		!multiplication per factor exp(-img*0.5*p*p*hhtt) 
		!p goes from -nn/2 to nn/2
		do ii=-nn/2, -1
			p=ii*twopi/(xmax-xmin)
	    		tfpsi(ii+nn/2+1)=tfpsi((ii+nn/2+1))*cexp(cmplx(0.,
     $				-0.5*p*p*hhtt))
     		end do	
	
		do ii= 1, nn/2
			p=ii*twopi/(xmax-xmin)
	    		tfpsi(ii+nn/2+1)=tfpsi((ii+nn/2+1))*cexp(cmplx(0.,
     $				-0.5*p*p*hhtt))
     		end do	
	
		signn= 1!fourier anti transfom
		call tfourier(nn, tfpsi ,signn, 
     $			nn, psitemp, xmin, xmax, hh)

		!multiplication per factor exp(-0.5*img*V_ll*hhtt) 
		!with V_ll= V(xmin+hh*(ll-1)) 
		do ll=1, nn 
		  vv=pot(xmin,hh,ll-1,ti,hhtt,jj-1, tf)
	    	  psitemp(ll)=psitemp(ll)*cexp(cmplx(0, -0.5*vv*hhtt))
		  temp=psitemp(ll)
		  psi(ll,jj+1)=temp
      		end do	

	end do
	
	!PRINTING psi(t+hhtt*(jj-1)) on file
	!jj goes from 1 to nntt
	!nntt psi
	!for each istant a different file is created
		
		open(unit = 20, file = "psi_time", 
     $			status = "unknown",
     $  		iostat=my_stat, iomsg=my_msg)

		if(my_stat /= 0) then
			print*, 'Open -psi_time_- failed with stat = '
     $		         , my_stat, ' msg = '//trim(my_msg)
		end if 
	
		do ii=1, nn !do cycle in order to print in the proper order
	     		write (20,*) xmin+(ii-1)*hh, 
     $				(abs(psi(ii, ist)), ist = 1, nntt)		
      		end do
		close(20)


		open(unit = 70, file = "averages", 
     $			status = "unknown",
     $  		iostat=my_stat, iomsg=my_msg)

		if(my_stat /= 0) then
			print*, 'Open -average_- failed with stat = '
     $		         , my_stat, ' msg = '//trim(my_msg)
		end if

		!COMPUTATION OF AVERAGE, MEAN SQUARE DISTANCE 
		!FROM AVERAGE AND NORM
		!For each istant
		!on a file, called "averages"
		
		do ist=1, nntt
			xxx=0.0
			nnn=0.0
			sig=0.0
			do ii=1,nn
			 !average
			 xxx=xxx+ (xmin+(ii-1)*hh)*hh*
     $					abs(psi(ii, ist))**2
			 !mean square distance from average
			 sig= sig + ((xmin+(ii-1)*hh)**2)*hh*
     $					abs(psi(ii, ist))**2
			 !norm
 			 nnn=nnn+ hh*abs(psi(ii, ist))**2
			end do
			xxx=xxx/nnn
			sig=sig/nnn
		  write(70,*) ti+(ist-1)*hhtt, xxx , nnn , sqrt(sig)
			
		end do
		print*, "f"
		
			


	stop
	end program harmosctime	
