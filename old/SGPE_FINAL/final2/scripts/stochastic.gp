#gnuplot
set terminal png enhanced
set output out
set xl "x (oscillator units)"
set yl "psi (oscillator units)"
set yr [0:100]
plot data u 1:4
