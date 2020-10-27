C=================================================
C
C       FUNCTION init_state (psi, nn, dd, debug, sepp, purr) 
C
C       .. Scalar Arguments ..
C       integer           nn
C	integer           dd
C	logical		  debug
C       logical		  sepp
C	logical		  purr
C       ..
C       .. Array Arguments ..
C	type(state)	psi
	
C
C
C	  Purpose
C	  =======
C
C	\details \b Purpose:
C	\verbatim
C 
C 	!This subroutine initializes a type state 
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
	!The state and the norms are printed on terminal 
	!if .debug. is .true.
C	
C	 Arguments:
C	  ==========
C


C	 \param[in] nn
C	 \verbatim
C          ii is integer
C	          number subsystems
C	 \endverbatimm
C
C	 \param[in] dd
C	 \verbatim
C          ti is integer
C	          number dimension of Hilbert space
C	 \endverbatim
C
C	 \param[in] debug
C	 \verbatim
C          debug is logical
C	          discretization step time grid
C	 \endverbatim
C
C	 \param[in] jj
C	 \verbatim
C          jj is integer
C	          the jj point is computed (counting ti as 0)
C	 \endverbatimm
C
C	 \param[in] tf
C	 \verbatim
C          tf is real*8
C	          Upper boundary time interval
C	 \endverbatim
C
C  =================================================C=================================================
C
C       SUBROUTINE tfourier(nn, vec,signn, nnkk, tfvec,xmin,xmax,hh)
C
C       .. Scalar Arguments ..
C       INTEGER      nn
C	integer	     signn
C	INTEGER      nnkk
C	real*8 	     xmin
C	real*8 	     xmax
C	real*8 	     hh
C       ..
C       .. Array Arguments ..
C       double complex         vec(nn)
C       double complex         tfvec(nn)
C       ..
C
C
C	  Purpose
C	  =======
C
C	\details \b Purpose:
C	\verbatim
C 
C	This subroutine performs a fourier transform of vec if sign=-1, 
C	an anti fourier transform if sign=1. The basis of Fourier space is 
C	cexp(cmplx(0., signn*twopi*x/(xmax-xmin)*float(ii)))
C	where ii is integer and goes from -nnkk/2, (nnkk-2)/2
C	\endverbatim
C	
C	 Arguments:
C	  ==========
C

C	 \param[in] nn
C	 \verbatim
C          nn is INTEGER
C	          The dimension of the vector to be transformed
C	 \endverbatim
C
C	 \param[in] signn
C	 \verbatim
C          signnn is INTEGER
C	          -1 is forward FT, +1 is backwards FT
C	 \endverbatim
C
C	 \param[in] nnkk
C	 \verbatim
C          nnkk is INTEGER
C	          The number of harmonics used
C	 \endverbatim
C
C	 \param[in] xmin
C	 \verbatim
C          xmin is real*8
C	          Lower boundary space interval
C	 \endverbatim
C
C	 \param[in] xmax
C	 \verbatim
C          xmin is real*8
C	          Upper boundary space interval
C	 \endverbatim
C
C	 \param[in] hh
C	 \verbatim
C          hh is real*8
C	          discretization step space grid
C	 \endverbatim
C
C	 \param[in] vec
C	 \verbatim
C         vec is double complex
C	          the vector which will be transformed
C	 \endverbatim
C
C	 \param[in] tfvec
C	 \verbatim
C         tfvec is double complex
C	          the vector which is transformed
C	 \endverbatim
C
C  =================================================
