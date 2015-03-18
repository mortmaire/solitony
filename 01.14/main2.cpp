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
#define complex complex<double>
using namespace std;
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
double max;//parametr maksymalnej wartości solitonów
int nn; //parametr kroku rozkładu statystycznego
parameters(double dx=0.1,double L=100,double gamma=0.01,double T=10,double mu=1,double g=0.01,double dt=0.01,double time=1000,double max=20,double nn=100):dx(dx),L(L),gamma(gamma),T(T),mu(mu),g(g),dt(dt),time(time),max(max),nn(nn){}
} p0;
class ffal;
class ssav;
class ssav{
char name[20];
fstream file;
ffal *psi;
char name1[20];
char name2[20];
public:
ssav(ffal* psi,int gpid,char ffrk):psi(psi){
sprintf(name,"../temp/%d%c",gpid,ffrk);
sprintf(name1,"%s.dat",name);
sprintf(name2,"%s.gp",name);
file.open(name1,ios::out|ios::ate);
file.close();
cout<<name1<<endl;
}
~ssav(void){
fstream f2(name2,ios::out);
//tutaj output z instrukcją gnuplota
f2.close();
}
void run(double time);

ssav(void);
};

class ffal{
public:
int i;
complex *tab;
parameters p;
double t;
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
}}
friend ostream& operator<<(ostream& out, ffal& ket);


void evolve(){
double tt=1.;
ssav save(this,gpid,ffrk);
if(ffrk=='a')system("setterm -cursor off");
while(t<p.time){
potential();
kinetic2();
potential();
t+=2*p.dt;
i+=2;
tt-=2*p.dt;
if(tt<=0){
tt=1.;
save.run(t);
cerr<<setfill('0')<<setw(7)<<fixed<<setprecision(2)<<t<<"\r";}
}
if(ffrk=='a')system("setterm -cursor on");
cerr<<endl;
}





};
void ssav::run(double time){
//dwie rzeczy - ma spisać całą tablicę wyników ORAZ rozkład statystyczny
double L=psi->p.L;
double dx=psi->p.dx;
int N=2*L/dx+1;
int M=psi->p.nn;
int tab[M];
double ctab[N];
double max=0;
double min;
for(int i=0;i<N;i++){
ctab[i]=abs(psi->tab[i]);
if(max<ctab[i])max=ctab[i];
if(i==0)min=ctab[i];
else if(min>ctab[i])min=ctab[i];
}//ta pętla służy ustaleniu max wart gęst.

for(int i=0;i<M;i++)tab[i]=0;
for(int i=0;i<N;i++){
int j=(ctab[i]-min)/(max-min)*M;
if(j==M)j--;//to zabezpiecza przed zakwalifikowaniem wart max poza tab
tab[j]++;}//mamy rozkład populacji
int mu=0;int pp=0;double s2=0;
//toto muszę zmodyfikować aby zaczynał od miejsca max rozkładu.
for(int i=0;i<M;i++)if(pp<tab[i]){pp=tab[i];mu=i;}
for(int i=0;i<M;i++)s2+=(i-mu)*(i-mu)*tab[i];
s2/=N;
double A=N;//*(max-min)/M;
//nie ma sensu dawać szerokości funkcji falowych jeśli nigdzie ich nie st.
double ctab2[M];
//teraz czas na yolo
double amp=A/(sqrt(2*M_PI*s2));
#define f(x) amp*exp(-(x-mu)*(x-mu)/(2*s2))
int n=0;
string name1=name;
name1+=".dat";
file.open(name1.c_str(),ios::out|ios::ate|ios::app);
for(double x=-L;x<=L;x+=dx){
stringstream bufor;
bufor<<setw(4)<<setfill(' ')<<time<<"\t"<<x<<"\t"<<ctab[n]<<"\t"<<arg(psi->tab[n]);
if(n<M)bufor<<"\t"<<A*n/N+min<<"\t"<<tab[n]<<"\t"<<f(n);
bufor<<"\n";
file<<bufor.str();
n++;
}file<<"\n\n";
file.close();
}




#define word(x) #x
int main(){
ffal psi;
psi.evolve();
}
