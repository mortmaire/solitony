#gnuplot
set terminal png enhanced
set output out
set pm3d map
set xl "time (oscillator units)"
set yl "x (oscillator units)"
splot data u 1:2:5
