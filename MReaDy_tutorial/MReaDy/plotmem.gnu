set terminal postscript eps enhanced color font 'Helvetica,10'
set xrange [1:]
set output 'plotmem.eps'
plot 'fort.85' u 1:2
