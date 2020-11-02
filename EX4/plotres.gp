# Set the output to a png file

set style line 1 \
    linecolor rgb '#0060ad' \
    linetype 1 linewidth 1.1 \
    pointtype 7 pointsize 1

set style line 2 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 1.1 \
    pointtype 5 pointsize 1

set style line 3 \
    linecolor rgb '#cc7722' \
    linetype 1 linewidth 1.1 \
    pointtype 3 pointsize 1

set style line 4 \
    linecolor rgb '#000000' \
    linetype 3 linewidth 1.1 dashtype 2 \
    pointtype 3 pointsize 0


set terminal postscript eps color colortext 
set key left top
set output 'O0.eps'
set title 'Multiplication Time - no opt' font "Helvetica, 40"
set xlabel 'Matrix size' font "Helvetica,30"
set ylabel 'Multiplication time [s]'font "Helvetica,30"
set logscale y
set logscale x
set grid
a=1
b=.5
f(x) = a*x**b
g(x) = c*x**d
h(x) = e*x**p

c=a
d=b
e=a
p=b

fit f(x) "result.dat" using 1:2 via a, b
fit g(x) "result.dat" using 1:3 via c,d
# fit h(x) "result.dat" using 1:4 via e,p
plot "result.dat" using 1:2 title 'StdLoop' with points linestyle 1, \
	 "result.dat" using 1:3 title 'InvLoop' with points linestyle 2, \
	 "result.dat" using 1:4 title 'Intrinsic' with points linestyle 3,\
	 f(x) with line linestyle 4,\
	 g(x) with line linestyle 4,\
	 # h(x) with line linestyle 4

set output 'O1.eps'
set title 'Multiplication Time - O1' font "Helvetica,40"
set xlabel 'Matrix size' font "Helvetica,30"
set ylabel 'Multiplication time [s]' font "Helvetica,30"
set logscale y
set logscale x
fit f(x) "result-O1.dat" using 1:2 via a, b
fit g(x) "result-O1.dat" using 1:3 via c,d
# fit h(x) "result-01.dat" using 1:4 via e,p
plot "result-O1.dat" using 1:2 title 'StdLoop' with points linestyle 1, \
	 "result-O1.dat" using 1:3 title 'InvLoop' with points linestyle 2, \
	 "result-O1.dat" using 1:4 title 'Intrinsic' with points linestyle 3,\
	 f(x) with line linestyle 4,\
	 g(x) with line linestyle 4,\
	 # h(x) with line linestyle 4

set output 'O2.eps'
set title 'Multiplication Time - O2' font "Helvetica,40"
set xlabel 'Matrix size' font "Helvetica,30"
set ylabel 'Multiplication time [s]' font "Helvetica,30"
set logscale y
set logscale x
fit f(x) "result-O2.dat" using 1:2 via a, b
fit g(x) "result-O2.dat" using 1:3 via c,d
# fit h(x) "result-02.dat" using 1:4 via e,p
plot "result-O2.dat" using 1:2 title 'StdLoop' with points linestyle 1, \
	 "result-O2.dat" using 1:3 title 'InvLoop' with points linestyle 2, \
	 "result-O2.dat" using 1:4 title 'Intrinsic' with points linestyle 3,\
	 f(x) with line linestyle 4,\
	 g(x) with line linestyle 4,\
	 # h(x) with line linestyle 4

set output 'O3.eps'
set title 'Multiplication Time - O3' font "Helvetica,40"
set xlabel 'Matrix size' font "Helvetica,30"
set ylabel 'Multiplication time [s]' font "Helvetica,30"
set logscale y
set logscale x
fit f(x) "result-O3.dat" using 1:2 via a, b
fit g(x) "result-O3.dat" using 1:3 via c,d
# fit h(x) "result-03.dat" using 1:4 via e,p
plot "result-O3.dat" using 1:2 title 'StdLoop' with points linestyle 1, \
	 "result-O3.dat" using 1:3 title 'InvLoop' with points linestyle 2, \
	 "result-O3.dat" using 1:4 title 'Intrinsic' with points linestyle 3,\
	 f(x) with line linestyle 4,\
	 g(x) with line linestyle 4,\
	 # h(x) with line linestyle 4

set output 'Ofast.eps'
set title 'Multiplication Time - Ofast' font "Helvetica,40"
set xlabel 'Matrix size' font "Helvetica,30"
set ylabel 'Multiplication time [s]' font "Helvetica,30"
set logscale y
set logscale x
fit f(x) "result-Ofast.dat" using 1:2 via a, b
fit g(x) "result-Ofast.dat" using 1:3 via c,d
# fit h(x) "result-0fast.dat" using 1:4 via e,p
plot "result-Ofast.dat" using 1:2 title 'StdLoop' with points linestyle 1, \
	 "result-Ofast.dat" using 1:3 title 'InvLoop' with points linestyle 2, \
	 "result-Ofast.dat" using 1:4 title 'Intrinsic' with points linestyle 3,\
	 f(x) with line linestyle 4,\
	 g(x) with line linestyle 4,\
	 # h(x) with line linestyle 4












