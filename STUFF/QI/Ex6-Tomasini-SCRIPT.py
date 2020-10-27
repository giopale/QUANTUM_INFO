import numpy as np
import sys
import os
import subprocess



nmin=1000
#lower bound dimension matrix
nmax=1000+250
#upper bound dimension matrix
stepp=250
#step of the increasing dimension

########################################################################################
#cycle to re-iterate the launch of the program
#that performs the Finite Difference Method
#for the harmonic oscillator with potential omega^2*x^2.
#The first kk eigenvalues are stored in a file
#called "first_k_eigenvalues-nn", with nn the number of
#points in the interval xmin, xmax
#The first kk eigenvalues are stored in a file
#called "first_k_eigenvectors-nn", with nn the number of
#points in the interval xmin, xmax
#kk,xmin and xmax can be modified in the program
#the eigenvectors are plotted and compared
#with the eigenvector of the harminic oscillator
#with potential 0.5*(omega_1)*x^2 where 
#0.5*(omega_1) is equal to omega^2
#######################################################################################	

for ii in range(nmin,nmax,stepp):
	f = open("MatDimension.txt","w")
	#opening the file with the dimension of the matrix
	f.write(str(ii))
	f.close()
	#write the dimension on the opened file
	subprocess.call("./a.out")
	#execute the fortran program

	kvec = open("first_k_eigenvec-%d" %ii ,"w")
	#open file which will collect data per fixed number of points
	
	nsp = open("first_k_eigenvec.txt","r")
	nsplines = nsp.readlines()		
	#reading first k eigenvec from file
		
	for qq in nsplines:
		kvec.write(qq)
		#copying data in file "first_k_eigenvec-%d"
		nsp.close()
		#closing file
	kvec.close()

	kval = open("first_k_eigenval-%d" %ii ,"w")
	#open file which will collect data per fixed number of points
	
	nsp = open("first_k_eigenval.txt","r")
	nsplines = nsp.readlines()		
	#reading first k eigenvec from file
		
	for qq in nsplines:
		kval.write(qq)
		#copying data in file "first_k_eigenvec-%d"
		nsp.close()
		#closing file
	kval.close()

	res = open("first_k_eigenvec-%d" %ii,"r")
	reslines = res.readlines()		
	#reading eigenvector data

	temp = open("resultstemp","w")
	for jj in reslines:
		temp.write(jj)
	#copying data in temporary file
	res.close()
	temp.close()
	#closing files

	os.system("gnuplot plotting1")
	#plot of first  eigenvector vs space

	newnamefit= "First_Eigenvalue.pdf"
	os.rename("tempfit.pdf",newnamefit)
	#renaming the pdf file of the fit

	os.system("gnuplot plotting2")
	#plot of second eigenvector vs space

	newnamefit= "Second_Eigenvalue.pdf"
	os.rename("tempfit.pdf",newnamefit)
	#renaming the pdf file of the fit

	os.system("gnuplot plotting3")
	#plot of third eigenvector vs space

	newnamefit= "Third_Eigenvalue.pdf"
	os.rename("tempfit.pdf",newnamefit)
	#renaming the pdf file of the fit
	
	os.system("gnuplot plotting3inv")
	#plot of third eigenvector vs space
	#the theorical one is inversed

	newnamefit= "Third_Eigenvalue_inverse.pdf"
	os.rename("tempfit.pdf",newnamefit)
	#renaming the pdf file of the fit
exit()
