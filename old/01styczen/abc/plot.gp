#LIST = system("ls -1 ?.dat")
#FILES = words(LIST)
#FILE(i)=word(LIST,i)
do for[i=100:2000:100]{
data = "stat2.dat"
czyjest = system("grep ^".i." \"stat2.dat\"")
if(czyjest!=""){
j=i/2
print j
set title "Width L = ".j
set xlabel "Temperature"
set ylabel "Solitons per length unit"
unset key
set terminal png enhanced size 1024,768

set output "T".j.".png"
plot data u 2:($1==i?$3:1/0)
}}
#> gnuplot 4.3 (development version in CVS on SourceForge):
#>
#> LIST = system("ls -1 out*")
#> FILES = words(LIST)
#> FILE(i) = word(LIST,i)
#>
#> plot for [i=1:FILES] FILE(i)
