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
set output 'Hist_1000_0100.pdf'
set title 'Normalized eigenvalue spacing' font "Helvetica, 40"
set xlabel 'Subsequent eval spacing' font "Helvetica,30"
set ylabel 'Frequency'font "Helvetica,30"
# set logscale y
# set logscale x
set grid
a=1
b=.5
f(x) = a*x**b


# fit f(x) "result.dat" using 1:2 via a, b

plot "Hist_1000_0100.dat" using 1:2 title 'histo' with points linestyle 1, \
	 # f(x) with line linestyle 4,\





