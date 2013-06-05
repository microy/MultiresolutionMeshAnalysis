# Grid
set grid

# Draw in a file
set terminal png small color

# Legend
set nokey

set xlabel "Resolution level"
set ylabel "Detail magnitude"
set title "Geometric Details"
set output "GeometricDetails00.png"
plot "coefficient-00.dat" using ($1):($2) with lines
set output "GeometricDetails01.png"
plot "coefficient-01.dat" using ($1):($2) with lines
set output "GeometricDetails02.png"
plot "coefficient-02.dat" using ($1):($2) with lines
set output "GeometricDetails03.png"
plot "coefficient-03.dat" using ($1):($2) with lines
set output "GeometricDetails04.png"
plot "coefficient-04.dat" using ($1):($2) with lines
