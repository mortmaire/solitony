set terminal png enhanced size 1024,768
numbers = "50 100 200 500"
set yrange [0:1.6];
set xlabel "Temperatura"
set ylabel "Gęstość solitonów w jednostce długości"
do for[i=1:words(numbers)]{
j=word(numbers,i)+0;
a=sprintf("%03d.png",j);
set output a;
set title sprintf("L = %d",word(numbers,i)*2);
plot "aaa.txt" u 2:($1==word(numbers,i)?$3:1/0) title "";
}
set yrange [0:2]
set output "zbiorczy1.png"
plot for[i=1:words(numbers)] "aaa.txt" u 2:($1==word(numbers,i)?$3:1/0) title sprintf("%d",word(numbers,i)+0);
set output "zbiorczy2.png";
set xrange [0:200];set yrange [0:1];
plot for[i=1:words(numbers)] "aaa.txt" u 2:($1==word(numbers,i)?$3:1/0) title sprintf("%d",word(numbers,i)+0);

