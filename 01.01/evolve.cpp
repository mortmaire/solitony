void evolve(){
fstream file;
fstream file2;
stringstream ss2;ss2<<"../temp/"<<gpid<<ffrk;
string name1=ss2.str()+".tmp.dat";
string name3=ss2.str()+"S.tmp.dat";

if(ffrk=='a')cerr<<name1<<endl;

file.open(name1.c_str(),ios::out);//czyszczę plik tmp.dat
file2.open(name3.c_str(),ios::out);
file.close();
file2.close();
int i=0;
int I=p.time/p.dt/1e3;
if(ffrk=='a')system("setterm -cursor off");
while(t<p.time){
potential();
kinetic2();
potential();
t+=2*p.dt;
i+=2;
if(i%(I/100)==0){cerr<<setfill('0')<<setw(7)<<fixed<<setprecision(2);if(ffrk=='a')cerr<<t<<"\r";}
if(i>I){i=0;
file.open(name1.c_str(),ios::out|ios::ate|ios::app);
file<<*this;
file.close();
file2.open(name3.c_str(),ios::out|ios::ate|ios::app);
file2>>*this;
file2.close();
}
}
if(ffrk=='a')system("setterm -cursor on ");
stringstream ss;
ss<<"T = "<<p.T;
char title[20];sprintf(title,"T = %.3f",p.T);
char gnuplot[1000];ffrk;

//#define write(x) 

sprintf(gnuplot,"ss=\"%s\"\n\
ss2=ss.\".tmp.dat\"\n\
t1=%f\n\
t2=%f\n\
T=%f\n\
L=%f\n\
set terminal png enhanced size (1024*2),(768)\n\
set output sprintf(\"%%s.png\",ss)\n\
set pm3d map\n\
set xl \"time (oscillator units)\"\n\
set yl \"x (oscillator units)\"\n\
set xrange [t1:t2]\n\
set multiplot layout 1,2 title sprintf(\"T = %%f\",T)\n\
splot ss2 u 1:2:3\n\
splot ss2 u 1:2:(sin($4/2.))\n\
do for [i=2:10]{\n\
unset multiplot\n\
set pm3d map \n\
set pal gray\n\
set xl \"time (oscillator units)\"\n\
set yl \"x (oscillator units)\"\n\
set output sprintf(\"%%s%%02d.png\",ss,i)\n\
set multiplot layout 1,2\n\
set yrange [-L:L]\n\
set cbrange [0:1]\n\
set title sprintf(\"T = %%f, sigma = %%.1f\",T,i/2.)\n\
ii=5+i-2\n\
splot ss2 u 1:2:ii w image\n\
unset pm3d\n\
set xl \"time (oscillator units)\"\n\
set yl \"n (number of solitons)\"\n\
set yrange [0:200]\n\
ii2=5+7+i\n\
unset label\n\
plot ss2 u 1:ii2\n\
unset multiplot}",
ss2.str().c_str(),0.,p.time,p.T,p.L);
string name2=ss2.str()+".gp";
fstream f2;
f2.open(name2.c_str(),ios::out);
f2<<gnuplot;
f2.close();
string systr="gnuplot "+name2;
system(systr.c_str());

int n=0;

}

