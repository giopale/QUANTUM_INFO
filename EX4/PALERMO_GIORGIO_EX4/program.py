import os
import numpy


fid=open("grid.dat","w+")

Nmin=10
Nmax=2500
numbers=numpy.logspace(numpy.log10(Nmin),numpy.log10(Nmax),num=10)
numbers=numpy.floor(numbers)

for item in numbers:
	fid.write("%s\n" % item)
fid.close()

opt_flags = ["","-O1", "-O2", "-O3", "-Ofast"]
for item in opt_flags:
	comp_comm="gfortran cos.f90 -o cos.out "+item
	os.system(comp_comm)
	command = "./cos.out " + item
	os.system(command)

os.system("gnuplot plotres.gp")



