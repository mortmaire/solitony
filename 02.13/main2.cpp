#include <iostream>
#include <cstdlib>
#include <complex>
#include <fftw3.h>
using namespace std;
#define complex complex<double>
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
parameters(double dx=0.1,double L=204.8,double gamma=0.01,double T=20,
double mu=1,double g=0.01,double dt=0.01,double time=1000):
dx(dx),L(L),gamma(gamma),T(T),mu(mu),g(g),dt(dt),time(time){}}p0;
class ffal{
parameters p;
double **tab;
int tt;
int N;
public:
ffal(parameters p=p0):p(p){
tt=0;
N=p.L/p.dx;
tab=new double*[N];
for(int i=0;i<N;i++){
tab[i]=new double[2];
tab[i][0]=0;
tab[i][1]=0;
}
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

fftw_complex *in, *out;
fftw_plan p0;//tutaj przysłaniam globalny parametr p0
in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
for(int i=0;i<N;i++){
in[i][0]=tab[i][0];
in[i][1]=tab[i][1];
}
p0 = fftw_plan_dft_1d(N, in, out, w, FFTW_ESTIMATE);
fftw_execute(p0);
if(w<0)
for(int i=0;i<N;i++){
tab[i][0]=out[i][0]/N;
tab[i][1]=out[i][1]/N;


}
else for(int i=0;i<N;i++){
tab[i][0]=out[i][0];
tab[i][1]=out[i][1];}

fftw_destroy_plan(p0);
fftw_free(in);
fftw_free(out);
}

void kinetic2(){//ewolucja kinetyczna, jedna z najważniejszych funkcji
complex i(0,1);
complex A(1,-p.gamma);
shift(-1);//przeskakujemy w przestrzeń pędów...
int L=p.L;
for(int j=0;j<N;j++){
double k;
if(j<N/2)k=2*j*M_PI/L;
else k=(j-N)*2*M_PI/L;
complex e=exp(-A*i*k*k*p.dt)*complex(tab[j][0],tab[j][1]);//używamy operatora ewolucji kinetycznej...
tab[j][0]=e.real();
tab[j][1]=e.imag();
}
shift(1);//...i zaraz wracamy w przestrzeń położeń.
}
void potential(){//ewolucja potencjalna, druga najważniejsza funkcja

complex i(0,1);
complex A(1,-p.gamma);
for(int j=0;j<N;j++){
double f=tab[j][0]*tab[j][0]+tab[j][1]*tab[j][1];
complex e=exp(-A*i*(p.g*f-p.mu)*p.dt)*complex(tab[j][0],tab[j][1]);
e-=i*sqrt(2*p.gamma*p.T)*eta()*p.dt;
tab[j][0]=e.real();
tab[j][1]=e.imag();
}}

void evolve(){

while(tt<p.time){
for(int i=0;i<50;i++){
potential();
kinetic2();
potential();
}
tt++;
//for(int i=0;i<N;i++)
//printf("%d\t%f\t%f\n",tt,i*p.dx,tab[i][0]*tab[i][0]+tab[i][1]*tab[i][1]);
//printf("\n");
}}

};







int main(){
double tt=time(0);
ffal psi(p0);
psi.evolve();
cout<<time(0)-tt<<endl;

//wyniki:
//133
//65 -O2
//65 -O2
//136

}
