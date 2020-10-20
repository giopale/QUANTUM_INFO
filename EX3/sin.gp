# Set the output to a png file
set terminal png
set output 'timevsizeO0.png'
set title 'Multiplication time - no optimization'
set xlabel 'Matrix size'
set ylabel 'Multiplication time [s]'
set logscale y
plot "time_resultsO0.txt" using 1:2 title 'StdLoop' with lines, \
	 "time_resultsO0.txt" using 1:3 title 'InvLoop' with lines, \
	 "time_resultsO0.txt" using 1:4 title 'Intrinsic' with lines

set output 'timevsizeO1.png'
set title 'Multiplication time - O1'
set xlabel 'Matrix size'
set ylabel 'Multiplication time [s]'
plot "time_resultsO1.txt" using 1:2 title 'StdLoop' with lines, \
	 "time_resultsO1.txt" using 1:3 title 'InvLoop' with lines, \
	 "time_resultsO1.txt" using 1:4 title 'Intrinsic' with lines

set output 'timevsizeO2.png'
set title 'Multiplication time - O2'
set xlabel 'Matrix size'
set ylabel 'Multiplication time [s]'
plot "time_resultsO2.txt" using 1:2 title 'StdLoop' with lines, \
	 "time_resultsO2.txt" using 1:3 title 'InvLoop' with lines, \
	 "time_resultsO2.txt" using 1:4 title 'Intrinsic' with lines

set output 'timevsizeO3.png'
set title 'Multiplication time - O3'
set xlabel 'Matrix size'
set ylabel 'Multiplication time [s]'
plot "time_resultsO3.txt" using 1:2 title 'StdLoop' with lines, \
	 "time_resultsO3.txt" using 1:3 title 'InvLoop' with lines, \
	 "time_resultsO3.txt" using 1:4 title 'Intrinsic' with lines