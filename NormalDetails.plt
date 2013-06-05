# Grid
set grid

# Draw in a file
set terminal png small color

# Legend
set nokey

set xlabel "Resolution level"
set ylabel "Detail magnitude"
set title "Normal Details"
set output "NormalDetails.png"
plot "Coefficients.dat" using ($1):($3) with lines
