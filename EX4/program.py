import os
import subprocess
import numpy


fid=open("grid.dat","w+")
numbers=numpy.arange(100,800,100)

for item in numbers:
	fid.write("%s\n" % item)

# max_dim=1600
# min_dim=10
# command="./Multiply.out "+ str(min_dim) +" " + str(max_dim)
# os.system(command)
# os.system("gnuplot plotres.gp")



