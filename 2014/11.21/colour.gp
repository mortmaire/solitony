#gnuplot

set terminal png enhanced size 1024,768
set output "out.png"
set pm3d map
set xl "time (oscillator units)"
set yl "x (oscillator units)"
set cbrange [0:30]
set xrange [0:120]
set yrange [-50:50]
splot "13081.tmp..dat" u 2:1:3
