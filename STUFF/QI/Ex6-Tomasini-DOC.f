C
C  =========== DOCUMENTATION ===========
C  Definition:
C ===========
C
C*********************************************************
C	 MODULE MATRICES
C*********************************************************
C

C	TYPE dcm
	
C	Arguments:
C
C	INTEGER nr		number of rows
C	INTEGER nr		number of columns
C	double complex	elem	elements of matrix
C	double precision eval	eigenvalues hermitian matrix(real)
C	double complex m_trace	trace
C	double complex m_det	determinat


C	  Purpose
C	  =======
C
C	\details \b Purpose:
C	\verbatim
C	Can be used to contain double complex matrices and their features.
C	The eigenvalues are double precision.
C	\endverbatim

C
C      SUBROUTINE init_hermmat (aa, aacheck, numrow, numcol,debug) 
C
C       .. Scalar Arguments ..
C       INTEGER          numrow
C	INTEGER          numcol	
C	LOGICAL		 debug
C
C      .. Array Arguments ..
C       type(dcm)          aa
C       type(dcm)          aacheck
C	
C
C  	Purpose
C  	=======
C
C	\details \b Purpose:
C	\verbatim
C
C	Purpose: solving Schrodinger equation for harmonic oscillator
C	in one-dimension: H=p^2+omega^2*x^2.
C	The Finite Difference method is exploited.
C	Natural units:
C	m=0.5
C	hbar=1.0
C	xmin,xmax of order 10^0
C	omega=1.0
C	Step of discretization=hh=(xmax-xmin)/(nn-1), with nn the number of points
C	in the interval.
C	diagonal entries_(ii)= (2+ omega*x^2*hh^4)/(hh^2)
C	x^2=(ii-nn/2)^2*h^2
C	(ii,ii +/- 1) entries = (-1)/(hh^2)
C	(others)=0
C	The factor 1/h^2 is omitted in order to avoid dividing per zero (small hh).
C	Then is restored.
C	The dimension of matrix should be equal to numrow and numcol,
C	if different, a warning message is printed. This debug is done 
C	if debug==.true. 
C
C	\endverbatim
C
C  	Arguments:
C 	 ==========
C

C 	\param[in] aa
C 	\verbatim
C          aa is type(dcm)
C          	To be initialized
C 	\endverbatim
C 	\param[in] aacheck
C 	\verbatim
C          aacheck is type(dcm)
C		A copy of the initialized matrix
C 	\endverbatim
C
C 	\param[in] numrow
C 	\verbatim
C          numrow is INTEGER 
C 	\endverbatim
C
C 	\param[in] numcol
C 	\verbatim
C          numcol is INTEGER 
C 	\endverbatim
C
C
C 	\param[in] debug
C 	\verbatim
C          numrow is LOGICAL
C 	\endverbatim
C
C==========================================================
C
C
C       SUBROUTINE print_matrices (AA, namefile)
C
C       .. Scalar Arguments ..
C       character(:)            namefile
C       ..
C      .. Array Arguments ..
C       type(dcm)          AA
C       ..
C
C  	Purpose
C  	=======
C
C	\details \b Purpose:
C	\verbatim
C
C	This subroutine prints a type(dcm) variable on a file, called
C	"namefile", at unit 40. 
C	Number of rows, columns, matrix elements, determinant, trace and eigenvalues
C	are printed. The matrix elements are printed in proper order.
C	\endverbatim
C
C  	Arguments:
C  	==========
C

C 	\param[in] AA
C 	\verbatim
C          AA is type(dcm)
C          	The variable that will be printed
C 	\endverbatim
C
C 	\param[in] namefile
C 	\verbatim
C          	namefile is character(:) 
C 	\endverbatim
C
C =====================================================
C
C
C       SUBROUTINE print_vector (vec, nn, namefile) 
C
C       .. Scalar Arguments ..
C       character(:)            namefile
C	integer			nn
C       ..
C      .. Array Arguments ..
C       double precision         vec(nn)
C       ..
C
C  	Purpose
C  	=======
C
C	\details \b Purpose:
C	\verbatim
C
C	This subroutine prints a double precision vector on a file, called
C	"namefile", at unit 80. 
C	\endverbatim
C
C  	Arguments:
C  	==========
C

C 	\param[in] vec
C 	\verbatim
C          vec is double precision, dimension(:) 
C          	The vector that will be printed
C 	\endverbatim
C 	\param[in] nn
C 	\verbatim
C          nn is INTEGER 
C 	\endverbatim
C 	\param[in] namefile
C 	\verbatim
C          	namefile is character(:) 
C 	\endverbatim
C  ===============================================================
C
C*********************************************************
C	 MODULE debugging
C*********************************************************
C
C       SUBROUTINE CHECKDIM(nn,nncheck,debug)
C
C       .. Scalar Arguments ..
C       INTEGER*2            nn
C	 INTEGER*2            nncheck
C	 LOGICAL	      debug
C       ..
C       .. Array Arguments ..
C       INTEGER*2          m(nn,nn)
C       INTEGER*2          mcheck(nn,nn)
C       ..
C  	Purpose
C 	 =======
C
C	\details \b Purpose:
C	\verbatim
C 
C 	CHECKDIM checks if the dimension of the matrix are correct. Firstly checks
C 	if the the dimension INTEGER*2 is above 10000. If it is so, it prints a WARNING
C 	message because it takes too much time. Secondly, it checks if nn is inferior  
C 	to 1, printing a WARNING message. Lastly, it checks if nn is actually nncheck
C 	as it should be. If it is not true, it prints a WARNING message, stating which 
C 	should be the dimension and the actual one. If everything goes well, it prints
C	"okay: right dimensions matrices", 
C 	and the actual dimension.

C	\endverbatim
C
C 	Arguments:
C  	==========
C

C 	\param[in] nn
C 	\verbatim
C          nn is INTEGER*2
C          The dimension of the squared arrays
C 	\endverbatim
C 	\param[in] qq
C 	\verbatim
C          qq is INTEGER*2
C	    	the number of the cycle, it should be nn=qq*100
C 	\endverbatim
C
C 	\param[inout] debug
C 	\verbatim
C          debug is LOGICAL
C 	\endverbatim
C
C=========================================================
C
C       SUBROUTINE CHECKMAT(nn,m,mcheck,debug)
C
C       .. Scalar Arguments ..
C       INTEGER*2            nn
C	 LOGICAL	      debug
C       ..
C       .. Array Arguments ..
C       INTEGER*2          m(nn,nn)
C       INTEGER*2          mcheck(nn,nn)
C       ..
C
C
C	  Purpose
C	  =======
C
C	\details \b Purpose:
C	\verbatim
C 
C	 CHECKMAT checks if the two imput matrices m and mcheck are equal. This is done
C	 via a do cycle which checks the equality between the two matrices element
C	 by element. For each entry which is not equal, the internal variable INTEGER*4 accum
C	 increseas by 1 (starting from 0). If accum in the end is bigger than 0,
C	 a WARNING message states how much entries between the two matrices are 
C	 different. Otherwise, it is printed "Same matrices".
C	\endverbatim
C	
C	 Arguments:
C	  ==========
C

C	 \param[in] nn
C	 \verbatim
C          nn is INTEGER*2
C	          The dimension of the squared arrays
C	 \endverbatim
C
C	 \param[in] 
C	 \verbatim
C	          m is INTEGER*2 ARRAY (nn,nn)
C	 \endverbatim
C
C	 \param[in] 
C	 \verbatim
C	          mcheck is INTEGER*2 ARRAY (nn,nn)
C	 \endverbatim
C
C	 \param[inout] debug
C	 \verbatim
C	          debug is LOGICAL
C	 \endverbatim
C
C  ===============================================================
C
C*********************************************************
C	 MODULE normalized_spacings
C*********************************************************
C=========================================================
C
C       SUBROUTINE normspac(aa,kk,normsp,debug)
C
C       .. Scalar Arguments ..
C       INTEGER		kk
C	 LOGICAL	debug
C       ..
C       .. Array Arguments ..
C       type(dcm)           aa
C       doubleprecision     normsp(aa%nr)
C       ..
C
C
C	  Purpose
C	  =======
C
C	\details \b Purpose:
C	\verbatim
C 
C	 normspac computes normalized spacings between eigenvalues 
C	 disposed in crescent order of a double complex matrix, which
C	 is inside type aa. The impute integer kk is a "locality parameter",
C	 i.e the number of spacings plus one considered in the computation
C	 of the norm for each spacing.
C	 The first thing computed are the spacings, then the norms.
C	 A if cycle secures to have (kk-1) in [1,nn-1] and correct if it is not,
C	 with a warning message.
C	 In a cycle over the spacings, the boundaries over which the norm is computed
C	 are assigned to each spacing, in a such a way to include kk-1 spacings.
C	 In case of even kk the interval is symmetric, in case of odd kk it is not.
C	 In case those boundaries exceed the global boundaries, they are translated
C	 to consider the interval between the global boundary and another point, such
C	 that to have kk-1 spacings in between. 
C	 Lastly, norms and normalized spacing are computed.
C	 If debug is .true., debugging is done (printing on file and on terminal).
C
C	\endverbatim
C	
C	 Arguments:
C	  ==========
C

C	 \param[in] aa
C	 \verbatim
C          aa is type(dcm)
C	          The type contain the eigenvalues in eval
C	 \endverbatim
C
C	 \param[in] 
C	 \verbatim kk
C	          kk is INTEGERk
C	 \endverbatim
C
C	 \param[out] 
C	 \verbatim
C	          normsp is doubleprecision(aa%nr)
C	 \endverbatim
C
C	 \param[inout] debug
C	 \verbatim
C	          debug is LOGICAL
C	 \endverbatim
C
C==================================================================
C*********************************************************
C	 PROGRAM eigenproblem
C*********************************************************
C=========================================================
C
C	  Purpose
C	  =======
C
C	\details \b Purpose:
C	\verbatim
C
C	Purpose: solving Schrodinger equation for harmonic oscillator
C	in one-dimension: H=p^2+omega^2*x^2.
C	The Finite Difference method is exploited.
C	Natural units:
C	m=0.5
C	hbar=1.0
C	xmin,xmax of order 10^0
C	omega=1.0
C	Step of discretization=hh=(xmax-xmin)/(nn-1), with nn the number of points
C	in the interval.
C	Given the matrix dimension from file "MatDimension.txt", initialize a type(dcm),
C	with inside the matrix of the problem.
C	diagonal entries_(ii)= (2+ omega*x^2*hh^4)/(hh^2)
C	x^2=(ii-nn/2)^2*h^2
C	(ii,ii +/- 1) entries = (-1)/(hh^2)
C	(others)=0
C	The factor 1/h^2 is omitted in order to avoid dividing per zero (small hh).
C	Then is restored.
C	Via the lapack routine zheev, the (double precision) eigenvalues are stored inside
C	the type, in crescent order. Moreover, eigenvalues are placed as columns in the matrix.
C	The first k eigenvalues are printed on file "first_k_eigenval.txt".
C	The first k eigenvectors are printed on file "first_k_eigenvectors.txt".
C	If debug==.true., debugging is done: checking matrix dimension, if matrix is
C	still the same, printing type(dcm) variables on file (including eigenvalues
C	and eigenvectors).
C	Morever, via the output iinfo, we check if zheev worked.
C
C	\endverbatim
C	  Authors:
C	  ========
C	
C	 \author Univ. of Padua
C	
C	 \date 13 November 2018

C  =====================================================================
