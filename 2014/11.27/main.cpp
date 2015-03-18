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
parameters(double dx=0.1,double L=50,double gamma=0.01,double T=500,double mu=22.41,double g=0.1,double dt=0.001,double time=60):dx(dx),L(L),gamma(gamma),T(T),mu(mu),g(g),dt(dt),time(time){}
} p0;

class ffal;
class memory{
double **tab;
double *mean;
parameters p;
int n;
int T;
int i;
public:
memory(parameters p=p0):p(p),i(0){
n=2*p.L/p.dx+1;
T=1e4;//'stała' globalna, mniej więcej tyle danych chcę mieć w pliku - wystarczy dla rozsądnej rozdzielczości.
tab=new double*[T];
mean=new double[n];
for(int j=0;j<n;j++)mean[j]=0;
for(int i=0;i<T;i++)tab[i]=new double[n];
}



void write(ffal *psi);


void draw(void){

fstream file;
stringstream ss2;ss2<<gpid<<ffrk<<".tmp.";
string name1=ss2.str()+"dat2";
file.open(name1.c_str(),ios::out);//czyszczę plik tmp.dat
for(int t0=0;t0<T;t0++){
int i0=0;
for(double x=-p.L;x<=p.L;x+=p.dx)file<<x<<"\t"<<t0*p.time/T<<"\t"<<abs(tab[t0][i0++])<<"\n";
file<<"\n";
}
file.close();

}



} *mem;

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
tab[j]*=exp(-A*i*(x*x/2+p.g*f*f-p.mu)*p.dt);
tab[j]-=i*sqrt(2*p.gamma*p.T)*eta()*p.dt;
x+=p.dx;
}
}
public:
void evolve(){
fstream file;
stringstream ss2;ss2<<gpid<<ffrk<<".tmp.";
string name1=ss2.str()+"dat";
file.open(name1.c_str(),ios::out);//czyszczę plik tmp.dat
file.close();
int i=0;
int I=p.time/p.dt/1e4;
//cout<<I<<endl;

#define ok cout<<"raz"<<endl;

while(t<p.time){
potential();
kinetic2();
potential();
t+=2*p.dt;
i+=2;
if(i>=I){

i=0;
cout<<t<<endl;
file.open(name1.c_str(),ios::out|ios::ate|ios::app);
file<<*this;

mem->write(this);//cancer

file.close();}
}
cout<<"Zakończono parametry g="<<p.g<<" mu="<<p.mu<<"."<<endl;
cout<<"Tworzę negatyw..."<<endl;
mem->draw();

stringstream ss;
ss<<"T = "<<p.T<<"K, "<<setw(1)<<"log(dt) = "<<log10(p.dt)<<", g = "<<p.g<<", {/Symbol m} = "<<p.mu;
string title=ss.str();
stringstream gnuplot;
gnuplot<<"set terminal png enhanced size 2048,768\n\
set output \""<<gpid<<".tmp.png\"\n\
set pm3d map interpolate 1,1\n\
set xl \"time (oscillator units)\"\n\
set yl \"x (oscillator units)\"\n\
set xrange [0:"<<p.time<<"]\n\
set yrange [-"<<p.L<<":"<<p.L<<"]\n\
set title '"<<title<<"'\n\
set multiplot layout 1,2\n\
splot \""<<name1<<"\" u 2:1:3\n\
splot \""<<name1<<"2\" u 2:1:3\n\
";
string name2=ss2.str()+"gp";
fstream f2;f2.open(name2.c_str(),ios::out);
f2<<gnuplot.str();
f2.close();
string systr="gnuplot "+name2;
system(systr.c_str());
/*
//cout<<gpid<<endl;
int n=0;
DIR *dir;

dir=opendir("./negatives/");
ok;
struct dirent *ent;
ok;
while ((ent = readdir (dir)) != NULL)n++;//cancer
ok;
closedir (dir);
stringstream ss3;ss3<<"mv ";
ss3<<gpid<<".tmp.png ./negatives/"<<setw(6)<<setfill('0')<<n<<ffrk<<".png";
system(ss3.str().c_str());
//cout<<ss3.str()<<endl;
*/
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

void memory::write(ffal *psi){
for(int j=0;j<n;j++){
tab[i][j]=psi->Rad(j);}
i++;
if(i==T){
cout<<"Uśrednianie wyników..."<<endl;
for(int j=9e3;j<1e4;j++)
for(int k=0;k<n;k++)mean[k]+=tab[j][k];
for(int k=0;k<n;k++)mean[k]/=1e3;
for(int i0=0;i0<T;i0++)for(int j=0;j<n;j++)tab[i0][j]-=mean[j];
}

}


int main(int argc, char** argv){
/*int a=fork();
int b=fork();
int c=fork();
gpid=getpid();

ffrk+=(a==0)+2*(b==0)+4*(c==0);
/*
for(int T=1000;T<3000;T+=100)
for(int i=0;i<120;i+=16){
*/
system("mv *tmp* ~/thrash");
gpid=getpid();
srand(time(0)+getpid());
p0.mu=20;
p0.g=0.001;
p0.L=40;
p0.T=400;
p0.dt=1e-3;
p0.time=60;
//for(int i=0;i<1000;i++){
memory mm;
mem=&mm;
ffal psi;
psi.evolve();
//}
cout<<gpid<<endl;
stringstream ss;
ss<<"eog "<<gpid<<".tmp.png";
system(ss.str().c_str());
}
