#!/usr/bin/gnuplot
set term png
set output "tmp/tmp.png"
set yrange [0:100]
plot "tmp.dat" using 1:3
