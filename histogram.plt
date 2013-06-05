# Grid
set grid

# Draw in a file
set terminal png small color

# Legend
set nokey

set xlabel "Wavelet magnitude"
set ylabel "Point number"
set title "Geometric Detail Histogram"

set output "histogram-00.png"
plot "histogram-00.dat" using ($1):($2) with lines

set output "histogram-01.png"
plot "histogram-01.dat" using ($1):($2) with lines

set output "histogram-02.png"
plot "histogram-02.dat" using ($1):($2) with lines

set output "histogram-03.png"
plot "histogram-03.dat" using ($1):($2) with lines

set output "histogram-04.png"
plot "histogram-04.dat" using ($1):($2) with lines

set output "histogram-05.png"
plot "histogram-05.dat" using ($1):($2) with lines

set output "histogram-06.png"
plot "histogram-06.dat" using ($1):($2) with lines

set output "histogram-07.png"
plot "histogram-07.dat" using ($1):($2) with lines

set output "histogram-08.png"
plot "histogram-08.dat" using ($1):($2) with lines

set output "histogram-09.png"
plot "histogram-09.dat" using ($1):($2) with lines
