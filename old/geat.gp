set terminal png enhanced size 1024,768
set output "h1.png"
set pm3d map
set title "Rozkład populacji dla różnych temperatur"
set xlabel "Temperature"
set ylabel "|{/Symbol Y}|"
set xrange [0:8000]
set yrange [0:70]
splot "heat.dat" i 0
set output "h2.png"
set cbrange [0:0.01]
set title "Różnica symetryczna dla różnych temperatur"
splot "heat.dat" i 1
set output "h3.png"
set title "Rozkład populacji dla T<100"
set xrange [1:100]
set cbrange [0:0.05]
splot "heat.dat" i 0
set output "h4.png"
set cbrange [0:0.01]
set title "Różnica symetryczna dla T<100"
splot "heat.dat" i 1

