#include <iostream>
#include <cstdlib>
#include <complex>
#include <fftw3.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>
#include <cmath>
#include <bmpfile.h>
#define complex complex<double>

using namespace std;
int t2;

string loc;
int gpid=getpid();
char ffrk='a';

class parameters{
public:
double dx;
double L;
double gamma;
double T;
double mu;
double g;
double dt;
double time;
parameters(double dx=0.1,double L=200,double gamma=0.01,double T=4000,double mu=1,double g=0.01,double dt=0.01,double time=1000):dx(dx),L(L),gamma(gamma),T(T),mu(mu),g(g),dt(dt),time(time){}
} p0;

class ffal;
class ffal{
public:
int i;
complex *tab;
parameters p;
double t;

public:
double Rad(int j){
double re=tab[j].real();
double im=tab[j].imag();
return sqrt(re*re+im*im);
}


ffal(parameters p=p0):p(p),t(0){
i=0;
t=0;
int n=2*p.L/p.dx+1;
tab=new complex[n];
for(int i=0;i<n;i++)tab[i]=0;
}
~ffal(){
delete[] tab;
}
complex eta(){
double u1,u2;
u1=rand();
u2=rand();
if(u1==0)u1=RAND_MAX;
u1/=RAND_MAX;
u2/=RAND_MAX;
double r=sqrt(-2*log(u1));
double t=2*M_PI*u2;
complex s=polar(r,t)/sqrt(2);
double delta=1/(p.dx*p.dt);
delta=sqrt(delta);
return delta*s;
}
ffal(const char* name,parameters p=p0):p(p),t(0){
fstream in;
in.open(name,ios::out);
double bufor;
int i=0;
int j=0;
while(in>>bufor){
if(i==0||i==1){i++;continue;}
if(i==2)tab[j].real()=bufor;
if(i==3)tab[j++].imag()=bufor;
i++;
if(i==4)i=0;
}
}
friend ostream& operator<<(ostream& out, ffal& ket);
private:
void shift(int w){
int N=2*p.L/p.dx+1;
fftw_complex *in, *out;
fftw_plan p0;//tutaj przysłaniam globalny parametr p0
in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
for(int i=0;i<N;i++){
in[i][0]=tab[i].real();
in[i][1]=tab[i].imag();
}
p0 = fftw_plan_dft_1d(N, in, out, w, FFTW_ESTIMATE);
fftw_execute(p0);
if(w<0)
for(int i=0;i<N;i++)tab[i]=complex(out[i][0]/N,out[i][1]/N);
else for(int i=0;i<N;i++)tab[i]=complex(out[i][0],out[i][1]);

fftw_destroy_plan(p0);
fftw_free(in);
fftw_free(out);
}

void kinetic2(){//ewolucja kinetyczna, jedna z najważniejszych funkcji
complex i(0,1);
complex A(1,-p.gamma);
shift(-1);//przeskakujemy w przestrzeń pędów...
int L=p.L;
int N=2*L/p.dx+1;
for(int j=0;j<N;j++){
double k;
if(j<N/2)k=j*M_PI/L;
else k=(j-N)*M_PI/L;
tab[j]*=exp(-A*i*k*k*p.dt);//używamy operatora ewolucji kinetycznej...
}
shift(1);//...i zaraz wracamy w przestrzeń położeń.
}
void potential(){//ewolucja potencjalna, druga najważniejsza funkcja
complex i(0,1);
complex A(1,-p.gamma);
int N=2*p.L/p.dx+1;
double x=-p.L;
for(int j=0;j<N;j++){
double f=abs(tab[j]);
tab[j]*=exp(-A*i*(p.g*f*f-p.mu)*p.dt);
tab[j]-=i*sqrt(2*p.gamma*p.T)*eta()*p.dt;
x+=p.dx;
}
}
public:
void evolve(){
fstream file;
stringstream ss2;ss2<<gpid<<ffrk<<".tmp.";
string name1=ss2.str()+"dat";

if(ffrk=='a')cerr<<name1<<endl;

file.open(name1.c_str(),ios::out);//czyszczę plik tmp.dat
file.close();
int i=0;
int I=p.time/p.dt/1e3;
int tt1=0;
while(t<p.time){
potential();
kinetic2();
potential();
t+=2*p.dt;
i+=2;
if(i%(I/100)==0){
cerr<<setfill('0')<<setw(7) <<fixed <<setprecision(2);
if(ffrk=='a')cerr<<t<<"\r";}
if(i>I){i=0;
file.open(name1.c_str(),ios::out|ios::ate|ios::app);
file<<*this;
file.close();
}
}

stringstream ss;
ss<<"T = "<<p.T;
string title=ss.str();
stringstream gnuplot;
gnuplot<<"set terminal png enhanced size 2048,768\n\
set output \""<<gpid<<ffrk<<".png\"\n\
set pm3d map\n\
set xl \"time (oscillator units)\"\n\
set yl \"x (oscillator units)\"\n\
set xrange ["<<t2<<":"<<p.time<<"]\n\
set yrange ["<<-p.L<<":"<<p.L<<"]\n\
set multiplot layout 1,2\n\
set title '"<<title<<"'\n\
splot \""<<name1<<"\" u 1:2:3\n\
splot \""<<name1<<"\" u 1:2:(cos($4))\n\
";
string name2=ss2.str()+"gp";
fstream f2;f2.open(name2.c_str(),ios::out);
f2<<gnuplot.str();
f2.close();
string systr="gnuplot "+name2;
system(systr.c_str());

int n=0;

}
};
ostream& operator<<(ostream& out, ffal& ket){
double L=ket.p.L;
double dx=ket.p.dx;
int n=2*L+1;
complex i0(0,1);
n/=dx;
out<<fixed<<setprecision(6);
int i=0;
stringstream bufor;
for(double x=-L;x<=L;x+=dx)bufor<<ket.t<<"\t"<<x<<"\t"<<abs(ket.tab[i])<<"\t"<<arg(ket.tab[i++])<<"\n";
bufor<<endl;
out<<bufor.str();
return out;
}


int main(int argc, char** argv){
system("setterm -cursor off");
if(argc>1)t2=atoi(argv[1]);
else t2=0;
gpid=getpid();
//int a=fork();
//int b=fork();
//int c=fork();
//ffrk+=(a==0)+2*(b==0)+4*(c==0);

if(argc>2)p0.time=atoi(argv[2]);
else p0.time=1000;
p0.T=((int)(ffrk-'a')+1)*10;
gpid++;
ffal psi;
psi.evolve();
system("setterm -cursor on");
}

