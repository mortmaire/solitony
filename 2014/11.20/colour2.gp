#gnuplot

set terminal png enhanced
set output "out2.png"
set pm3d map
set xl "time (oscillator units)"
set yl "x (oscillator units)"
set cbrange [0:30]
splot "tmp2.dat" u 2:1:3
