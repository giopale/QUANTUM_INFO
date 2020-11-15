# Set the output to a png file

set style line 1 \
    linecolor rgb '#0060ad' \
    linetype 1 linewidth 1.1 \
    pointtype 7 pointsize 1.5

set style line 2 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 1.1 \
    pointtype 5 pointsize 1.5

set style line 3 \
    linecolor rgb '#cc7722' \
    linetype 1 linewidth 1.1 \
    pointtype 3 pointsize 1.5

set style line 4 \
    linecolor rgb '#000000' \
    linetype 3 linewidth 1.1 dashtype 2 \
    pointtype 3 pointsize 0


set terminal postscript eps color colortext 
set key right top

set output 'Hist_2000_0450.eps'
set title 'Normalized eigenvalue spacing' font "Helvetica, 30"
set xlabel 's_i' font "Helvetica,30"
set ylabel 'Relative frequency'font "Helvetica,30"
set bmargin 5
# set logscale y
# set logscale x
set grid
a=1
al=-1.3
b=.5
be=1.5
ap=a
alp=al
bp=b
bep=be
f(x) = (a*x**al)*exp(-b*x**be)
g(x) = (ap*x**alp)*exp(-bp*x**bep)


fit f(x) "Hist_2000_0450.dat" using 1:2 via a, al, b, be
fit g(x) "Hist_2000_0450_diag_mat.dat" using 1:2 via ap, alp, bp, bep

plot "Hist_2000_0450.dat" using 1:2 title 'Hermitian' with points linestyle 1, \
     f(x) with line linestyle 4,\
	 "Hist_2000_0450_diag_mat.dat" using 1:2 title 'Real Diagonal' with points linestyle 3, \
	 g(x) with line linestyle 4,\








