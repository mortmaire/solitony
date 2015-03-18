public:
ffal(parameters p=p0):p(p){
tt=0;
N=p.L/p.dx;
T=p.time+1;
tab=new double**[T];

for(int i=0;i<T;i++){
tab[i]=new double*[N];
for(int j=0;j<N;j++){
tab[i][j]=new double[2];
tab[i][j][0]=0;
tab[i][j][1]=0;
}}}

~ffal(){
for(int i=0;i<T;i++){
for(int j=0;j<N;j++)delete[] tab[i][j];}
for(int i=0;i<T;i++)delete[] tab[i];
delete[] tab;
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
complex e=exp(-A*i*k*k*p.dt)*complex(tab[tt][j][0],tab[tt][j][1]);
//używamy operatora ewolucji kinetycznej...
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

void save(const char* name){
fstream file(name,ios::out);
for(int i=0;i<T;i++){
for(int j=0;j<N;j++)
file<<i<<" "<<j<<" "<<tab[i][j][0]<<endl;
file<<endl;
}
file.close();
}
void load(const char* name){
string buffer;
fstream file(name,ios::in);
while(getline(file,buffer)){
if(buffer=="")continue;
int t,k;
double r;
stringstream bff;bff<<buffer;
bff>>t>>k>>r;
tab[t][k][0]=r;
}
}
