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

int **stat;

public:
ffal(parameters p=p0):p(p){
tt=0;
N=p.L/p.dx;
T=p.time+1;
tab=new double**[T];
stat=new int*[T];
for(int i=0;i<T;i++){
stat[i]=new int[100];
for(int j=0;j<100;j++)stat[i][j]=0;
tab[i]=new double*[N];
for(int j=0;j<N;j++){
tab[i][j]=new double[2];
tab[i][j][0]=0;
tab[i][j][1]=0;

}}
//system("setterm -cursor off");
}

~ffal(){
for(int i=0;i<T;i++){
for(int j=0;j<N;j++)delete[] tab[i][j];
delete[] stat[i];
}
for(int i=0;i<T;i++)delete[] tab[i];
delete[] tab;
delete[] stat;
//system("setterm -cursor on");
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
double u1,u2;
u1=rand();
u2=rand();
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
in[i][0]=tab[tt][i][0];
in[i][1]=tab[tt][i][1];
}



p0 = fftw_plan_dft_1d(N, in, out, w, FFTW_ESTIMATE);
fftw_execute(p0);
if(w<0)
for(int i=0;i<N;i++){
tab[tt][i][0]=out[i][0]/N;
tab[tt][i][1]=out[i][1]/N;


}
else for(int i=0;i<N;i++){
tab[tt][i][0]=out[i][0];
tab[tt][i][1]=out[i][1];}

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
complex e=exp(-A*i*k*k*p.dt)*complex(tab[tt][j][0],tab[tt][j][1]);//używamy operatora ewolucji kinetycznej...
tab[tt][j][0]=e.real();
tab[tt][j][1]=e.imag();
}
shift(1);//...i zaraz wracamy w przestrzeń położeń.
}
void potential(){//ewolucja potencjalna, druga najważniejsza funkcja

complex i(0,1);
complex A(1,-p.gamma);
for(int j=0;j<N;j++){
double f=tab[tt][j][0]*tab[tt][j][0]+tab[tt][j][1]*tab[tt][j][1];

complex e=exp(-A*i*(p.g*f-p.mu)*p.dt)*complex(tab[tt][j][0],tab[tt][j][1]);
e-=i*sqrt(2*p.gamma*p.T)*eta()*p.dt;
tab[tt][j][0]=e.real();
tab[tt][j][1]=e.imag();
}}

double statystyka(){
//po pierwsze, chcę zamienić wszystkie wartości (x,y) na (r,t)
for(int i=0;i<T;i++)
for(int j=0;j<N;j++){
double *t2=tab[i][j];
complex z(t2[0],t2[1]);
t2[0]=abs(z);
}
double tab2[N];
double max=0;
for(int i=0;i<N;i++)tab2[i]=0;
for(int i=0;i<T;i++)
for(int j=0;j<N;j++){
//tab2[j]+=tab[i][j][0];
if(tab[i][j][0]>max)max=tab[i][j][0];
}

for(int i=0;i<T;i++)
for(int j=0;j<N;j++){
int k=tab[i][j][0]*100./max;

if(k<100&&k>=0)stat[i][k]++;
}

double srednia[100];

for(int i=0;i<100;i++)srednia[i]=0;
for(int i=T/2;i<T;i++)
for(int j=0;j<100;j++)srednia[j]+=stat[i][j];
for(int i=0;i<100;i++)srednia[i]/=(T/2*p.L/204.8);//normowanie do standardowych
//danych szerokości 204.8

double A=0;
int mu=0;
double s2=0;
double n0=0;

for(int i=0;i<100;i++)if(A<srednia[i]){A=srednia[i];mu=i;}
for(int i=mu;i<100;i++){
s2+=(i-mu)*(i-mu)*srednia[i];n0+=srednia[i];
}
s2/=n0;

#define f(x) A*exp(-(x-mu)*(x-mu)/(2*s2))

//teraz dwie wersje wydarzeń: ALBO szukam maksimum różnicy rozkładów
//normalnego i rzeczywistego, albo sprawdzam, gdzie się przecinają

//zacznijmy od tego drugiego. Chciałbym, aby rozkład rzeczywisty był
//nie większy od normalnego niż o 5%.

//albo i nie. W wyższych temperaturach taka strategia niezbyt dobrze działa

double roznica[100];
double mmax=0;double trigger=0;

for(int i=0;i<100;i++){
roznica[i]=srednia[i]-f(i);
if(roznica[i]>mmax){mmax=roznica[i];trigger=i;}
}

trigger*=max/100;

for(int i=0;i<T;i++)for(int j=0;j<N;j++){
if(tab[i][j][0]>trigger)tab[i][j][1]=1;
else tab[i][j][1]=0;//tab[i][j][0]/trigger;

}

//clean();



char name[30];
char name2[30];

sprintf(name2,"../temp/%d%cS.dat",gpid,ffrk);
fstream file;
fstream file2(name2,ios::out);
if(ffrk=='a')cerr<<gpid<<endl;
//tutaj jest taki mnożnik 200. to w celu ładnego wykresu w gnuplocie

double mnoznik=((N/1000)*100.)/N;

double sum=0;
bool flag=1;

for(int t1=0;t1<p.time;t1++){
if(t1%2000==0){
int aa=t1/2000;
if(aa>0)file.close();
sprintf(name,"../temp/%d%c%d.dat",gpid,ffrk,aa);
file.open(name,ios::out|ios::app);
}
for(int i=0;i<N;i++){

file<<t1<<"\t"<<i*mnoznik<<"\t"<<tab[t1][i][0]*tab[t1][i][0]<<endl;
//"\t"<<tab[t1][i][1]<<endl;
//if(t1>T/2&&tab[t1][i][1]<1.&&flag==1){sum++;flag=0;}
//if(t1>T/2&&tab[t1][i][1]==1.&&flag==0)flag=1;

}
file<<endl;



//for(int i=0;i<100;i++)file2<<i/100.*max<<"\t"<<stat[t1][i]<<endl;
//file2<<endl<<endl;
}
file.close();
for(int i=0;i<100;i++){
file2<<p.T<<"\t"<<i/100.*max<<"\t"<<srednia[i]<<"\t"<<f(i)<<"\t"<<roznica[i]<<endl;
}
file2.close();

//chciałbym wygenerować skrypt gnuplotowy do rysowania rzeczy

gnuplot(name,name2);

return 0;//sum/p.L/(T/2);
}

void gnuplot(char* n1,char* n2){
char script[1024];

sprintf(script,READ(
set terminal png size 1024,768;
unset key;

set pm3d map;
set title "T = %.01f, L = %d";
set xlabel "time";
set ylabel "space position";
n0="../temp/%d%c";
do for[i=0:0]{
set pal color;
name = n0.i;
n2 = name.".dat";
set output name.".png";
set xrange [i*2000:(i+1)*2000];
print name;
splot n2 u 1:2:3;
//set pal gray;
//set output "../temp/".name.i."g".".png";
//set xrange [i*2000:(i+1)*2000];
//splot n2 u 1:2:4;
}
set autoscale;
unset pm3d;
nn2 = n0."S.dat";
set output n0."S.png";
set key;
plot nn2 u 2:3 title "real distribution", nn2 u 2:4 title "Gauss distribution", nn2 u 2:5 title "Real - Gauss";
),p.T,(int)p.L/100*100,gpid,ffrk);

char name[30];
sprintf(name,"../temp/%d%c.gp",gpid,ffrk);
fstream file(name,ios::out);
file<<script;
file.close();
char comm[30];
sprintf(comm,"gnuplot %s",name);
//int bb=system(comm);
}

void evolve(){



int ii=0;int t0=time(0);
int i0=0;
int isum=0;
int ile=0;
double isr;
cout<<setprecision(2)<<fixed;
while(tt<p.time){

if(ffrk=='a')cerr<<setprecision(2)<<fixed<<"\t\t\t\r"<<tt+1<<"\t"<<i0/50.<<" cycles per second, ETA: "<<(p.time-tt)/isr<<" s";
for(int i=0;i<50;i++){
ii++;
if(t0!=time(0)){i0=ii;t0=time(0);ii=0;ile++;isum+=i0/50.;
isr=isum/(double)ile;}

potential();
kinetic2();
potential();
}
for(int i=0;i<N;i++){
tab[tt+1][i][0]=tab[tt][i][0];
tab[tt+1][i][1]=tab[tt][i][1];
}
tt++;
}
if(ffrk=='a')cerr<<endl;
}



};






int main(){

p0.time=10000;
cout<<gpid<<endl;

int a=fork();
int b=fork();
int c=fork();

ffrk+=(a==0)+2*(b==0)+4*(c==0);

int s=system("setterm -cursor off");

int i=0;
int I=ffrk-'a';


for(p0.T=21;p0.T<100;p0.T+=1){

if(i++%8!=I)continue;

ffal psi(p0);
psi.evolve();
psi.statystyka();
gpid++;
}

s=system("setterm -cursor on");

//dla FFTW na -O2 mamy 66 sekund/1000 obrotów, ~820 kroków na sekundę
//generowanie pliku tekstowego zajmuje ~40 sekund

}
