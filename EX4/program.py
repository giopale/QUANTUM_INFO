import os
# import subprocess
import numpy


fid=open("grid.dat","w+")
numbers=numpy.logspace(1.,numpy.log10(2500),num=10)
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



