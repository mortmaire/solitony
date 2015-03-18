#gnuplot
set terminal png enhanced size 640,480
set output out
set multiplot layout 2,1
set pm3d map
#set zrange [0:200]
set xl "time (oscillator units)"
set yl "x (oscillator units)"
splot data1 u 1:2:5
splot data2 u 1:2:5
