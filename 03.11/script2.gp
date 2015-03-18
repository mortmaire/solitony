set terminal png enhanced size 1024,768
set output "ciastko.png"

list(i) = word(system('ls -1B ../ss/*.a'),i)
size = system('ls ../ss/*.a | wc -l');
set xlabel "|{/Symbol Y}|"
set ylabel "Liczba zliczeń"
set xzeroaxis lt 1 lc 0
set grid
do for[i=1:size]{
file=list(i)
name = system(sprintf("echo %s | awk -F'[^0-9]*' '$0=$2'",file))
x=0;
x=name+x;
set title sprintf("Rozkład populacji, T = %d",x)
#set output sprintf("../ss/T%04da.png",x)
c="ciastko"
c= sprintf("../ss/T%04d.a",x)
d= sprintf("../ss/T%04d.b",x)
print c
#plot c
set title sprintf("Różnica symetryczna, T = %d",x)
#set output sprintf("../ss/T%04db.png",x)
plot d
set title sprintf("Rozkład populacji, T = %d",x)
set output sprintf("../ss/T%04dc.png",x)
plot c title "Populacja",d title "Różnica"
}
