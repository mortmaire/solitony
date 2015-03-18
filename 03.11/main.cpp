#include <iostream>
#include <cstdlib>
#include <complex>
#include <fftw3.h>
#include <stdint.h>
#include <fstream>
#include <sys/types.h>
#include <unistd.h>
#include <iomanip>
#define complex complex<double>
#include <string.h>
#define READ(...) # __VA_ARGS__
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
parameters(double dx=0.1,double L=204.8,double gamma=0.01,double T=5,
double mu=1,double g=0.01,double dt=0.01,double time=1000):
dx(dx),L(L),gamma(gamma),T(T),mu(mu),g(g),dt(dt),time(time){}}p0;

class ffal{
double ***tab;
parameters p;
int tt;
int N;
int T;
int S;
int *stat;

#include "core.cpp"
double delta;
int s_min;
double max;


void statystyka(const char* name=0){
S=100;
stat=new int[S];
for(int i=0;i<S;i++)stat[i]=0;
max=0;
for(int i=T/2;i<T;i++)
for(int j=0;j<N;j++){
tab[i][j][0]=abs(complex(tab[i][j][0],tab[i][j][1]));
if(tab[i][j][0]>max)max=tab[i][j][0];
tab[i][j][1]=0;
}
//wstępna statystyka
for(int i=T/2;i<T;i++)
for(int j=0;j<N;j++){
int k=tab[i][j][0]*S/max;
if(k==S)k=S-1;
stat[k]++;
}
int n=T/2*N;
s_min=0;
int s_max=S-1;
while(1){
if(stat[s_min]>n/S/S)break;//jeśli jest większe od 1%
s_min++;
}
while(1){
if(stat[s_max]>n/S/S)break;
s_max--;
}
double min=max/S*s_min;
max=max/S*s_max;
for(int i=0;i<S;i++)stat[i]=0;//resetujemy statystykę, jeszcze raz
n=0;
delta=(max-min)/S;
for(int i=T/2;i<T;i++)
for(int j=0;j<N;j++){
if(tab[i][j][0]<min||tab[i][j][0]>max)continue;
int k=(tab[i][j][0]-min)/delta;
if(k==S)k--;
stat[k]++;
n++;}
if(name==0)
for(int i=0;i<S;i++)cout<<(s_min+i)*max/S<<"\t"<<stat[i]<<endl;
else{
fstream file(name,ios::out);
for(int i=0;i<S;i++)file<<(s_min+i)*max/S<<"\t"<<stat[i]<<endl;
file.close();
}
}

double antisym(const char* name=0){
int n=0;
double w=0;
double maxi=0;
double smax=0;
for(int i=0;i<S;i++){
if(smax<stat[i]){maxi=i;smax=stat[i];}
n+=stat[i];
w+=stat[i]*i;
}
w/=n;
w=(int)(w+0.5);

for(int i=0;i<w;i++){
if(2*w-i>=S){stat[i]=0;continue;}
stat[i]-=stat[(int)(2*w-i)];

}

for(int i=w;i<S;i++)stat[i]=0;
if(name==0)for(int i=0;i<S;i++)cout<<(s_min+i)*max/S<<" "<<stat[i]<<endl;
else{
fstream file(name,ios::out);
for(int i=0;i<S;i++)file<<(s_min+i)*max/S<<" "<<stat[i]<<endl;
file.close();
}
return (maxi-w)*delta;
}

};

int choice(int i){
int r_min=1;
int r_max=8192;
int r_size=r_max-r_min+1;
int step=r_size/2;
if(i==0)return r_min;
if(i==1)return r_max;
int j=2;
while(step>=1){
int x=step;
while(x<r_max){
if(i==j)return x;
x+=2*step;
j++;
}
step/=2;
}
}

char* name(int gpid,char a){
stringstream ss;
ss<<"./temp/"<<gpid<<a<<".dat";
string s=ss.str();
int n=s.length();
char* ret=new char[n];
strcpy(ret,s.c_str());
return ret;
}

int main(){
p0.time=10000;
int a=fork();
//int b=fork();
//int c=fork();
ffrk+=(a==0);//+2*(b==0)//+4*(c==0);//+8*(d==0);

p0.T=1;
if(ffrk=='a'){p0.dx=0.01;p0.L=163.84;}
else {p0.dx=0.001;p0.L=131.072;}
cout<<name(gpid,ffrk)<<endl;

ffal psi;
psi.evolve();
psi.save(name(gpid,ffrk));

}
