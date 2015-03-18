#include <iostream>
#include <cstdlib>
#include <complex>
#include <cmath>
#include <stdint.h>
#include <fftw3.h>
using namespace std;
#define write(...) # __VA_ARGS__
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
class ffal;
class ffal{
parameters p;
double ***przebieg;
int tt;//krok czasowy
int N;
public:
ffal(parameters p=p0):p(p){
N=p.L/p.dx;
tt=0;
int T=p.time;
int L=p.L/p.dx;
przebieg=new double**[T];
for(int i=0;i<T;i++)przebieg[i]=new double*[L];
for(int i=0;i<T;i++)for(int j=0;j<L;j++){przebieg[i][j]=new double[2];
przebieg[i][j][0]=0;przebieg[i][j][1]=0;
}
}
~ffal(){
int T=p.time;
int L=p.L/p.dx;

for(int i=0;i<T;i++)for(int j=0;j<L;j++)delete[] przebieg[i][j];
for(int i=0;i<T;i++)delete[] przebieg[i];
delete[] przebieg;
}
private:

void shift(double tab[][2],int w=-1){
/*
for(int i=0;i<n;i++){
out[i][0]=0;
out[i][1]=0;
for(int j=0;j<n;j++){
complex U=exp(complex(0.,s*2.*M_PI/n*j*i));
complex kk=complex(in[j][0],in[j][1]);
out[i][0]+=kk.real();
out[j][0]+=kk.imag();
}
if(s==1){out[i][0]/=n;out[i][1]/=n;}
}
*/

//in=(double(*)[2])in2;
//out=(double(*)[2])out2;
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
tab[i][1]=out[i][1];

}

fftw_destroy_plan(p0);
fftw_free(in);
fftw_free(out);
}

double rand(){
static uint64_t *state=new uint64_t(time(0));
uint32_t c=(*state)>>32, x=(*state)&0xFFFFFFFF;
*state = x*((uint64_t)4294883355U) + c;
uint32_t a = x^c;
if(a!=0)return a/(double)0xFFFFFFFF;
else return 1.;
}

complex eta(){
double u1=rand(),u2=rand();
double r=sqrt(-2*log(u1));
double t=2*M_PI*u2;
complex s=polar(r,t)/sqrt(2);
double delta=1/(p.dx*p.dt);
delta=sqrt(delta);
return delta*s;
}

inline double ab2(double *tab){
return tab[0]*tab[0]+tab[1]*tab[1];
}

inline complex cmpx(double *tab){
return complex(tab[0],tab[1]);
}

void step(){
for(int i=0;i<N;i++)for(int j=0;j<2;j++)
przebieg[tt+1][i][j]=przebieg[tt][i][j];
tt++;
double **tab;
tab=new double*[N];
for(int i=0;i<N;i++){
tab[i]=new double[2];
tab[i][0]=0;
tab[i][1]=0;
}



complex I(0,1);
complex A(1,-p.gamma);
for(int i=0;i<50;i++){
for(int j=0;j<N;j++){
double *pp=przebieg[tt][j];
double f=ab2(pp);

complex e=exp(-A*I*(p.g*f-p.mu)*p.dt)*cmpx(pp);
e*=cmpx(pp);e-=I*sqrt(2*p.gamma*p.T)*eta()*p.dt;
pp[0]=e.real();
pp[1]=e.imag();

}
double tab[N][2];
for(int j=0;j<N;j++){
tab[j][0]=przebieg[tt][j][0];
tab[j][1]=przebieg[tt][j][1];
}


shift(tab);
for(int j=0;j<N;j++){
double k;
if(j<N/2)k=j*2*M_PI/p.L;
else k=(j-N)*2*M_PI/p.L;
complex e=exp(-A*I*k*k*p.dt)*cmpx(tab[j]);
tab[j][0]=e.real();
tab[j][1]=e.imag();
}
shift(tab,1);

for(int j=0;j<N;j++){
przebieg[tt][j][0]=tab[j][0];
przebieg[tt][j][1]=tab[j][1];
}

for(int j=0;j<N;j++){
double *pp=przebieg[tt][j];
double f=ab2(pp);
complex e=exp(-A*I*(p.g*f-p.mu)*p.dt)*cmpx(pp);
e-=I*sqrt(2*p.gamma*p.T)*eta()*p.dt;
pp[0]=e.real();
pp[1]=e.imag();
}
}

}
public:
void evolve(){
while(tt<p.time){
step();
for(int i=0;i<N;i++)
printf("%d\t%f\t%f\n",tt,i*p.L/N,ab2(przebieg[tt][i]));
printf("\n");
}


}



};


int main(){
p0.time=1000;
ffal psi(p0);
psi.evolve();
}
