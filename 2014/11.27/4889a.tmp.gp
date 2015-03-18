set terminal png enhanced size 2048,768
set output "4889.tmp.png"
set pm3d map interpolate 1,1
set xl "time (oscillator units)"
set yl "x (oscillator units)"
set xrange [0:60]
set yrange [-40:40]
set title 'T = 400K, log(dt) = -3, g = 0.001, {/Symbol m} = 20'
set multiplot layout 1,2
splot "4889a.tmp.dat" u 2:1:3
splot "4889a.tmp.dat2" u 2:1:3
