# Grid
set grid

# Draw in a file
set terminal png small color

# Legend
set nokey

set xlabel "Resolution level"
set ylabel "Detail magnitude"
set title "Color Details"
set output "ColorDetails.png"
plot "Coefficients.dat" using ($1):($4) with lines
