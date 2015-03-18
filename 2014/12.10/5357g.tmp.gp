set terminal png enhanced size 2048,768
set output "5357g.png"
set pm3d map
set xl "time (oscillator units)"
set yl "x (oscillator units)"
set xrange [1000:3000]
set yrange [-200:200]
set multiplot layout 1,2
set title 'T = 7'
splot "5357g.tmp.dat" u 1:2:3
splot "5357g.tmp.dat" u 1:2:(cos($4))
