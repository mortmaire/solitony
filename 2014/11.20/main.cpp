#include <iostream>
#include <cstdlib>
#include <complex>
#include <fftw3.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#define complex complex<double>
using namespace std;

class parameters{
public:
double dx;
double L;
double gamma;
double T;
double mu;
double g;
double dt;
parameters(double dx=0.1,double L=20,double gamma=0.01,double T=20,double mu=22.41,double g=0.1,double dt=0.01):dx(dx),L(L),gamma(gamma),T(T),mu(mu),g(g),dt(dt){}
} p0;

class ffal{
complex *tab;
parameters p;
double t;
public:
ffal(parameters p=p0):p(p){
t=0;
int n=2*p.L/p.dx+1;
tab=new complex[n];
for(int i=0;i<n;i++){tab[i]=0;}
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
ffal(const char* name,parameters p=p0){
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
tab[j]*=exp(-A*i*(x*x/2+p.g*f*f-p.mu)*p.dt);
tab[j]-=i*sqrt(2*p.gamma*p.T)*eta()*p.dt;
x+=p.dx;
}
}
public:
void evolve(double time){
double t0=t+time;
fstream file;
file.open("tmp2.dat",ios::out);//czyszczę plik tmp.dat
file.close();
int i=0;
int I=1/p.dt/1000;
cout<<I<<endl;
while(t<t0){
potential();
kinetic2();
potential();
t+=2*p.dt;
i+=2;
if(i>=I){
i=0;
cout<<p.T<<"\t"<<t<<"\t"<<p.dt<<endl;
file.open("tmp2.dat",ios::out|ios::ate|ios::app);
file<<*this;
file.close();}
}

}
};
ostream& operator<<(ostream& out, ffal& ket){
double L=ket.p.L;
double dx=ket.p.dx;
out<<fixed<<setprecision(6);
int i=0;
stringstream bufor;
for(double x=-L;x<=L;x+=dx)bufor<<x<<"\t"<<ket.t<<"\t"<<abs(ket.tab[i++])<<"\n";
bufor<<endl;
out<<bufor.str();
return out;
}



int main(){
p0.dt=0.001;
p0.T=300;
ffal psi(p0);
psi.evolve(300);

system("gnuplot colour2.gp");
stringstream ss;
ss<<"mv out.png T="<<setw(3)<<setfill('0')<<p0.T<<","<<setw(1)<<"p_dt="<<-log10(p0.dt)<<".png";
system(ss.str().c_str());

}
