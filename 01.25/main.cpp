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
#define write(...) # __VA_ARGS__
using namespace std;
int gpid=getpid();
char ffrk='a';
double thr=0;

double sol_sr;

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
int nn; //parametr kroku rozkładu statystycznego
parameters(double dx=0.1,double L=100,double gamma=0.01,double T=20,double mu=1,double g=0.01,double dt=0.01,double time=1000,double nn=50):dx(dx),L(L),gamma(gamma),T(T),mu(mu),g(g),dt(dt),time(time),nn(nn){}
} p0;
class ffal;
class ssav;
class ssav{
char name[20];
fstream file;
ffal *psi;
char name1[20];
char name2[20];
parameters p;
double sol_srednia;
double n_sr;
public:
ssav(ffal* psi,int gpid,char ffrk,parameters p=p0):psi(psi),p(p){
sol_srednia=0;
n_sr=0;
sprintf(name,"../temp/%d%c",gpid,ffrk);
sprintf(name1,"%s.dat",name);
sprintf(name2,"%s.gp",name);
file.open(name1,ios::out);
file.close();
//if(ffrk=='a')cerr<<name1<<endl;
}
~ssav(void){
fstream file(name2,ios::out);
stringstream ss;
ss<<write(
s = "%s";
s1=s.".dat";
set title "T = %f";
set pal color;
set terminal png enhanced size 1024,768;
set output s."1.png";
set pm3d map;
splot s1 u 1:2:3;
set output s."2.png";
set palette gray;
splot s1 u 1:2:5;
);
char s2[1000];

sprintf(s2,ss.str().c_str(),name,p.T);
double tt1=p0.T;
file<<s2;
file.close();
string as2="gnuplot ";
as2+=name2;
system(as2.c_str());

sol_sr=sol_srednia/n_sr/(2*p.L);

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
if(ffrk=='a')cerr<<setfill('0')<<setw(7)<<fixed<<setprecision(2)<<t<<"\r";}
}
if(ffrk=='a')system("setterm -cursor on");
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
double mu=0;int pp=0;double s2=0;
//toto muszę zmodyfikować aby zaczynał od miejsca max rozkładu.
for(int i=0;i<M;i++)if(pp<tab[i]){pp=tab[i];mu=i;}

//for(int i=0;i<M;i++){mu+=i*tab[i];}mu/=N;

int n0=0;
for(int i=mu;i<M;i++){s2+=(i-mu)*(i-mu)*tab[i];n0+=tab[i];}
s2/=n0;
double A=n0*2;//*(max-min)/M;
//nie ma sensu dawać szerokości funkcji falowych jeśli nigdzie ich nie st.
double ctab2[M];
//teraz czas na yolo
double amp=A/(sqrt(2*M_PI*s2));
#define f(x) amp*exp(-(x-mu)*(x-mu)/(2*s2))
double delta=0;int k=0;
for(int i=0;i<mu;i++){//zbadaj, gdzie jest maks delty
if(tab[i]-f(i)>delta){k=i;delta=tab[i]-f(i);}
}
if(k>0){
double thr2=k*(max-min)/M+min;
if(thr<thr2)thr+=0.1;
if(thr>thr2)thr-=0.1;
}



string name1=name;
name1+=".dat";
file.open(name1.c_str(),ios::out|ios::ate|ios::app);
int sol_count=0;bool flag=1;



for(int n=0;n<N;n++){
stringstream bufor;
bool soliton=(ctab[n]>thr?1:0);
bufor<<setw(4)<<setfill(' ')<<time<<"\t"<<n*dx-L<<"\t"<<ctab[n]<<"\t"<<arg(psi->tab[n])<<"\t"<<soliton;

if(!soliton&&flag){sol_count++;flag=0;}
if(soliton)flag=1;
if(n==N-1){bufor<<"\t"<<sol_count<<"\n";
if(time>p.time/2.){sol_srednia+=sol_count;n_sr+=1.;}

}
bufor<<"\n";
file<<bufor.str();
}
//file<<endl;
file.close();
}





int main(int argc,char** argv){
srand(time(0));
int a=fork();
int b=fork();
int c=fork();
ffrk+=(a==0)+2*(b==0)+4*(c==0);
string name;name+=ffrk;name+=".dat";
fstream file(name.c_str(),ios::out|ios::app|ios::ate);
p0.L=0;
while(p0.L<500){
if(p0.L==0)p0.L=50;
else if(p0.L==50)p0.L=100;
else if(p0.L==100)p0.L=200;
else p0.L=500;
p0.T=0;
p0.time=2000;
int i=0;
while(p0.T<5000){
i%=8;
if(p0.T==0)p0.T=1;
else if(p0.T==1)p0.T=5;
else if(p0.T==5)p0.T=10;
else p0.T+=10;
if(i++!=ffrk-'a')continue;
ffal psi;
psi.evolve();
file<<2*p0.L<<"\t"<<p0.T<<"\t"<<sol_sr<<endl;
cout<<endl<<p0.L<<"\t"<<p0.T<<endl;
gpid++;
}
}}


