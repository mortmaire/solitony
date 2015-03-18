#include <iostream>
#include <cstdlib>
#include <fftw3.h>
#include <complex>
#include <cmath>
using namespace std;
#define complex complex<double>

int n=1024;
void fftw(complex *tab,int w){
fftw_complex *in, *out;
fftw_plan p0;
in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
for(int i=0;i<n;i++){
in[i][0]=tab[i].real();
in[i][1]=tab[i].imag();
}
p0 = fftw_plan_dft_1d(n, in, out, w, FFTW_ESTIMATE);
fftw_execute(p0);
if(w==-1)for(int i=0;i<n;i++)
tab[i]=complex(out[i][0]/n,out[i][1]/n);
else
for(int i=0;i<n;i++)
tab[i]=complex(out[i][0],out[i][1]);

fftw_destroy_plan(p0);
fftw_free(in);
fftw_free(out);


}

void naive(complex* tab,int w){
complex out[n];
for(int i=0;i<n;i++){
out[i]=0;
for(int j=0;j<n;j++){
out[i]+=tab[j]*exp(complex(0,w*2*M_PI*i*j/(double)n));
}
if(w==-1)out[i]/=n;
}
for(int i=0;i<n;i++)tab[i]=out[i];

}


void hybrid1(complex *tab,int w){
complex out[n];
for(int i=0;i<n;i++){
out[i]=0;
complex e=exp(complex(0,w*2*M_PI*i/(double)n));
for(int j=0;j<n/2;j++){
out[i]+=exp(complex(0,w*4*M_PI*i*j/(double)n))*(tab[2*j]+e*tab[2*j+1]);
}
if(w==-1)out[i]/=n;
}
for(int i=0;i<n;i++)tab[i]=out[i];
}
int ii=0;
void hybrid2(complex *tab,int w){
complex out[n];
int N=32;
for(int i=0;i<n/N;i++){
complex sum[N];
for(int j=0;j<N;j++)sum[j]=0;
for(int j=0;j<n/N;j++){
for(int k=0;k<N;k++){
sum[k]+=tab[N*j+k]*exp(complex(0,w*2*M_PI*N*i*j/(double)n));
}
}
for(int j=0;j<N;j++){
for(int k=0;k<N;k++){
out[n/N*j+i]+=sum[k]*exp(complex(0,w*2*M_PI*k*(n/N*j+i)/n));
}
}
}
for(int i=0;i<n;i++){
if(w==-1)tab[i]=out[i]/(double)n;
else tab[i]=out[i];
}
}





int main(){


complex tab[n];
for(int i=0;i<n;i++)tab[i]=complex(i,n-i);
int N=n*100;
int tt=time(0);
while(tt==time(0));
tt=time(0);int i=0;
while(tt==time(0)){
fftw(tab,-1);
fftw(tab,1);
i++;
}
cout<<"FFTW: "<<i<<" per second"<<endl;
/*
tt=time(0);
N=100;
for(int i=0;i<N;i++){
//naive(tab,-1);
//naive(tab,1);
}
cout<<"Naive: "<<(double)N/(time(0)-tt)<<" per second"<<endl;
*/
cout<<tab[500];
hybrid2(tab,-1);
cout<<tab[500];
hybrid2(tab,1);
cout<<tab[500]<<endl;
tt=time(0);
while(tt==time(0));
tt=time(0);i=0;
while(tt==time(0)){
hybrid2(tab,-1);
hybrid2(tab,1);
i++;
}
cout<<"Hybrid2: "<<i<<" per second"<<endl;
}
