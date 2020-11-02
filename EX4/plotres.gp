# Set the output to a png file
set terminal png
set output 'out1.png'
set title 'Multiplication Time'
set xlabel 'Matrix size'
set ylabel 'Multiplication time [s]'
set logscale y
set logscale x
plot "time_ByRows.txt" using 1:2 title 'StdLoop' with lines, \
	 "time_ByCols.txt" using 1:2 title 'InvLoop' with lines, \
	 "time_Intrinsic.txt" using 1:2 title 'Intrinsic' with lines

# set output 'timevsize_parallel_O1.png'
# set title 'Multiplication time parallel - O1'
# set xlabel 'Matrix size'
# set ylabel 'Multiplication time [s]'
# plot "time_results_parallel_O1.txt" using 1:2 title 'StdLoop' with lines, \
# 	 "time_results_parallel_O1.txt" using 1:3 title 'InvLoop' with lines, \
# 	 "time_results_parallel_O1.txt" using 1:4 title 'Intrinsic' with lines

# set output 'timevsize_parallel_O2.png'
# set title 'Multiplication time parallel - O2'
# set xlabel 'Matrix size'
# set ylabel 'Multiplication time [s]'
# plot "time_results_parallel_O2.txt" using 1:2 title 'StdLoop' with lines, \
# 	 "time_results_parallel_O2.txt" using 1:3 title 'InvLoop' with lines, \
# 	 "time_results_parallel_O2.txt" using 1:4 title 'Intrinsic' with lines

# set output 'timevsize_parallel_O3.png'
# set title 'Multiplication time parallel - O3'
# set xlabel 'Matrix size'
# set ylabel 'Multiplication time [s]'
# plot "time_results_parallel_O3.txt" using 1:2 title 'StdLoop' with lines, \
# 	 "time_results_parallel_O3.txt" using 1:3 title 'InvLoop' with lines, \
# 	 "time_results_parallel_O3.txt" using 1:4 title 'Intrinsic' with lines