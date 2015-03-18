#ifndef PSI_BODY_H
#define PSI_BODY_H



void parameters::show(void){
cerr<<"PARAMETRY:"<<endl;
cerr<<"dx = "<<dx<<endl;
cerr<<"L  = "<<L<<endl;
cerr<<"n  = "<<n<<endl;
cerr<<"I  = "<<I<<endl;
cerr<<"N  = "<<N<<endl;
cerr<<"g0 = "<<g<<endl;
cerr<<"tau= "<<tau<<endl;
cerr<<"t0  = "<<t<<endl;
cerr<<"gamma="<<gamma<<endl;
cerr<<"mu = "<<mu<<endl;
cerr<<"t1 = "<<t1<<endl;
}

bool parameters::variable(string in, double &var, string str){

int pos=in.find(str);
if(pos!=-1){
int len=str.length();
string s=in.substr(pos);
string sub=s.substr(len);
len=sub.find(" ");
sub=sub.substr(1,len);
double a=atof(sub.c_str());
if(a!=-1)var=a;
else{cerr<<"Błąd parametrów! Czyżby za duży odstęp?"<<endl;exit(1);}
}
}

bool parameters::variable(string in, int &var, string str){
int pos=in.find(str);
if(pos!=-1){
string s=in.substr(pos);
int len=str.length();
string sub=s.substr(len);
len=sub.find(" ");
sub=sub.substr(1,len);
double a=atoi(sub.c_str());
if(a!=-1)var=a;
else{cerr<<"Błąd parametrów! Czyżby za duży odstęp?"<<endl;exit(1);}
}
}


bool parameters::fill(int argc, char** argv){

stringstream ss;
for(int i=1;i<argc;i++)ss<<argv[i]<<" ";
string s0=ss.str();
replaceAll(s0, ",", ".");
replaceAll(s0, "  ", " ");
replaceAll(s0, " = ", "=");

//uzupełniam tutaj wszystkie parametry na podstawie inputu

variable(s0,dx,_(dx));//krok iteracyjny
variable(s0,L,_(L));//szerokość sandboxu (2xL)
variable(s0,n,_(n)); //numer stanu wzbudzonego
variable(s0,N,_(N)); //ilość cząstek
variable(s0,g,"g0"); //coupling constant
variable(s0,tau,_(tau)); //parametr czasu urojonego
variable(s0,I,_(I)); //ilość iteracji
variable(s0,t,"t0"); //czas wykonywania SGPE
variable(s0, dt,_(dt)); //krok czasowy
variable(s0,gamma,_(gamma)); //parametr gamma
variable(s0,T,_(T)); //temperatura
variable(s0,mu,_(mu)); //potencjał chemiczny
variable(s0,t1,_(t1)); //czas wykonywania równań hybrydowych

show();

}


























bool psi::hybrid(psi& con){
complex A(1.,-p.gamma);
int N = show_NL();//wielkość tablicy
complex tab[2][N]; //mamy dwie funkcje
double L=p.L;
double x0=-L;
double dx=p.dx;
double dt=p.dt;
double t=p.t1;
for(int i=0;i<N;i++){
complex a=val(x0);
tab[0][i]=a;
tab[1][i]=con(x0);//W KOŃCU UŻYŁEM TEGO OPERATORA
x0+=dx;
}
double dk=dx*M_PI/L;
complex i(0,1);



for(double t0=0;t0<=t;t0+=dt){

shift(tab[0],-1);
shift(tab[1],-1);
for(int rh=0;rh<2;rh++)
for(int j=0;j<N;j++){
double k;
if(j<N/2)k=j*M_PI/L;
else k=(j-N)*M_PI/L;
tab[rh][j]*=exp(-A*i*k*k*p.dt/2);
}


shift(tab[0],1);
shift(tab[1],1);


double x=-L;
for(int j=0;j<N;j++){
complex f[2];
f[0]=tab[0][j]*conj(tab[1][j]);
f[1]=tab[1][j]*conj(tab[0][j]);
complex xi=eta();
complex e[2];
e[0]=exp(-A*i*(x*x/2+p.g*f[0]-p.mu+sqrt(2*i*p.g)*xi.real())*p.dt);
e[1]=exp(-A*i*(x*x/2+p.g*f[1]-p.mu+sqrt(2*i*p.g)*xi.imag())*p.dt);
tab[0][j]=e[0]*tab[0][j];
tab[1][j]=e[1]*tab[1][j];
x+=dx;
}

x=-L;
int N=show_NL();
for(int j=0;j<N;j++){
complex eta0=eta();
tab[0][j]-=i*sqrt(2*p.gamma*p.T)*eta0*p.dt;
tab[1][j]-=i*sqrt(2*p.gamma*p.T)*eta0*p.dt;
x+=p.dx;
}


if(isfinite(tab[0][show_NL()/2].real())==0)return 0;

}





*this=psi(tab[0],p);
con=psi(tab[1],p);
return 1;
}



psi psi::create_conjugate(void){
int N=show_NL();
complex tab[N];
double x=-p.L;
for(int i=0;i<N;i++){
tab[i]=val(x);
x+=p.dx;
}
psi fun2(tab,p);
return fun2;
}



psi::psi(parameters p,double v):p(p){
for(double x=-p.L;x<=p.L;x+=p.dx)add(x,v);
}

double f(double x,double L);

double psi::iterate(){
int I=p.I;
int n=p.n;
int m;
if(n==0)m=1;
else m=n;
psi a[m+1];complex b[m+1];
//caveat - może się wysypać na innych kompilatorach
for(int i0=0;i0<=n;i0++){
cerr<<"Stan energetyczny "<<i0<<endl;
p.n=i0;
psi fn(f,p);
for(int i=0;i<I;i++){
//cerr<<i<<endl;
for(int i1=0;i1<i0;i1++)b[i1]=scalar_product(a[i1]);
complex p0=0;
complex p1=0;
complex p2=0;
double dx=p.dx;
double L=p.L;
for(double x=-L-dx;x<L;x+=dx){
p0=p1;
p1=p2;
p2=fn.val(x+dx);
for(int i1=0;i1<i0;i1++)p2-=b[i1]*a[i1].val(x+dx);
complex K=p2+p0-2*p1;
K/=-2*dx*dx;
complex V=x*x*p1;
V/=2;
double g=p.g;
double tau=p.tau;
complex G=p1*p1*conj(p1)*g;
complex M=-p.mu;
complex l=p1-(tau)*(K+V+G+M);
fn.add(x,l);
}
fn.norm();
a[i0]=fn;
(*this)=fn;
}
}
return energy_product(*this);

}




double psi::entropy(void){
complex f0,f1,f2;
f1=0;
double entropia=0;
double vmax=0;
for(double x=-p.L+p.dx;x<p.L;x+=p.dx){
if(vmax<abs(val(x)))vmax=abs(val(x));
}
if(vmax==0)return 0;
for(double x=-p.L+p.dx;x<p.L;x+=p.dx){
f0=f1;
f1=f2;
f2=val(x)/vmax;
complex d2=f2-f0;
complex d1=f1-f0;
double v1v2=abs(d1*d2)+2;
double modules=(abs(d2*d2)+4)*(abs(d1*d1)+1);
modules=sqrt(modules);
double p_i=v1v2/modules;

entropia+=-p_i*log(p_i);

}
return entropia;
}


double psi::hybrid(int powt){
double t0=0;
for(int i=0;i<powt;i++){
t0+=hybrid();}
return t0/powt;
}







bool psi::ifex(){//funkcja sprawdzająca czy plik o danej
//nazwie istnieje w folderze


fstream fs;
fs.open(name().c_str(),ios::in);
bool b= fs.is_open();
fs.close();
return b;
}

string psi::name(void){
stringstream ss;
ss<<"./baza/gamma"<<p.gamma<<"g"<<p.g<<"mu"<<p.mu<<".txt";
//cout<<p.g<<endl;
return ss.str();
}




double psi::hybrid(){
complex A(1.,-p.gamma);
int N = show_NL();//wielkość tablicy
complex tab[2][N]; //mamy dwie funkcje
double L=p.L;
double x0=-L;
double dx=p.dx;
double dt=p.dt;
double t=p.t1;
for(int i=0;i<N;i++){
complex a=val(x0);
tab[0][i]=a;
tab[1][i]=a;
x0+=dx;
}
double dk=dx*M_PI/L;
complex i(0,1);



for(double t0=0;t0<=t;t0+=dt){

shift(tab[0],-1);
shift(tab[1],-1);
for(int rh=0;rh<2;rh++)
for(int j=0;j<N;j++){
double k;
if(j<N/2)k=j*M_PI/L;
else k=(j-N)*M_PI/L;
tab[rh][j]*=exp(-A*i*k*k*p.dt/2);
}


shift(tab[0],1);
shift(tab[1],1);


double x=-L;
for(int j=0;j<N;j++){
complex f[2];
f[0]=tab[0][j]*conj(tab[1][j]);
f[1]=tab[1][j]*conj(tab[0][j]);
complex xi=eta();
complex e[2];
e[0]=exp(-A*i*(x*x/2+p.g*f[0]-p.mu+sqrt(2*i*p.g)*xi.real())*p.dt);
e[1]=exp(-A*i*(x*x/2+p.g*f[1]-p.mu+sqrt(2*i*p.g)*xi.imag())*p.dt);
tab[0][j]=e[0]*tab[0][j];
tab[1][j]=e[1]*tab[1][j];
x+=dx;
}

x=-L;
int N=show_NL();
for(int j=0;j<N;j++){
complex eta0=eta();
tab[0][j]-=i*sqrt(2*p.gamma*p.T)*eta0*p.dt;
tab[1][j]-=i*sqrt(2*p.gamma*p.T)*eta0*p.dt;
x+=p.dx;
}


if(isfinite(tab[0][show_NL()/2].real())==0)return t0;

}





*this=psi(tab[0],p);
add_conj(tab[1]);
return t;
}

void psi::add_conj(complex *tab){
double x=-p.L;
int N=show_NL();
for(int i=0;i<N;i++){
value_c[x]=tab[i];
x+=p.dx;
}
}

complex psi::eta(){
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


void psi::fourier_iterate_stochastic(){
int i1=0;
complex A(1,-p.gamma);
//cout<<A<<endl;
int N = show_NL();
complex tab[N];
double L=p.L;
double x0=-L;
double dx=p.dx;
double dt=p.dt;
double t=p.t;
for(int i=0;i<N;i++){
tab[i]=val(x0);
x0+=dx;
}
double dk=dx*M_PI/L;
complex i(0,1);
for(double t0=0;t0<=t;t0+=dt){
//cerr<<t0<<ta<<p.dt<<endl;
shift(tab,-1); //otrzymuję tablicę wartości w przestrzeni pędów

for(int j=0;j<N;j++){
double k;
if(j<N/2)k=j*M_PI/L;
else k=(j-N)*M_PI/L;
tab[j]*=exp(-A*i*k*k*p.dt/2);
}

shift(tab,1); //otrzymuję tablicę wartości w przestrzeni położeń
 
double x=-L;
for(int j=0;j<N;j++){
double f=abs(tab[j]);
f*=f;
tab[j]*=exp(-A*i*(x*x/2+p.g*f-p.mu)*p.dt);
tab[j]-=i*sqrt(2*p.gamma*p.T)*eta()*p.dt;
x+=dx;
}

}

*this=psi(tab,p);

}


void psi::shift(complex *tab,int w){
int N=show_NL();
fftw_complex *in, *out;
fftw_plan p0;
in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
for(int i=0;i<N;i++){
in[i][0]=tab[i].real();
in[i][1]=tab[i].imag();
}

p0 = fftw_plan_dft_1d(N, in, out, w, FFTW_ESTIMATE);
fftw_execute(p0);
double x=-p.L;
if(w<0)
for(int i=0;i<N;i++){
tab[i]=complex(out[i][0]/N,out[i][1]/N);
x+=p.dx;
}
else for(int i=0;i<N;i++){
tab[i]=complex(out[i][0],out[i][1]);
x+=p.dx;
}

fftw_destroy_plan(p0);
fftw_free(in); fftw_free(out);
}


psi::psi(complex *tab, parameters p):p(p){
int N=show_NL();
double x=-p.L;
for(int i=0;i<N;i++){
add(x,tab[i]);
x+=p.dx;
}
}



void psi::prnt(string name){
ofstream fs(name.c_str(),ofstream::out | ofstream::trunc);
fs<<(*this);
fs.close();
}




psi::psi(string name, parameters p):p(p){

ifstream fs(name.c_str());
double x,re,im;
string line;
while(getline(fs,line)){
istringstream iss(line);
iss>>x>>re>>im;
complex z(re,im);
value[x]=z;
}
fs.close();
}






complex psi::scalar_product(psi &ket){
complex integral=0;
double L=p.L;
double dx=p.dx;
for(double x=-L;x<L;x+=dx){
complex con=conj(ket(x));
con*=operator()(x);
integral+=con*dx;
}
return integral;
}

double psi::energy_product(psi &ket2){
psi &ket1 = (*this);
double L=ket1.show_L();
double dx=ket1.show_dx();
complex integral=0;
for(double x=-L;x<=L;x+=dx){
complex kinetic=ket2(x+dx)+ket2(x-dx)-ket2(x)-ket2(x);
kinetic/=-dx*dx*2;
complex potential=ket2(x)*x*x/2.;
complex con=conj(ket1(x));
integral+=(con*kinetic+con*potential)*dx;
}
return abs(integral);
}

ostream& operator<<(ostream& os, psi& ket)
{
double L=ket.show_L();
double dx=ket.show_dx();
for(double x=-L;x<=L;x+=dx)os<<fixed<<setprecision(6)<<x<<"\t"<<ket(x).real()<<"\t"<<ket(x).imag()<<"\t"<< abs(ket(x))<<endl;
return os;
}


void psi::add(double x, complex f){//funkcja dodawania kolejnych wartości
value[x]=f;
}

psi::psi(double (*f)(double,double),parameters p):p(p){
int n=p.n;
int m=n00;
n00=n;
double L=p.L;
double dx=p.dx;
for(double x=-L;x<=L;x+=dx){
complex z((*f)((n+1)*x-n*L,L),0);
add(x,z);
}
norm();
n00=m;//fragment wrzucania odpowiedniego rozw. analitycznego
}

void psi::norm(){//normalizator. Use wisely.
int N=p.N;
double L=p.L;
double dx=p.dx;
double mod=0;
for(double x=-L;x<=L;x+=dx){
double ap=abs(val(x));
mod+=ap*ap*dx;}
double M=sqrt(mod/N);
//cerr<<N<<endl;
for (map<double,complex>::iterator it=value.begin(); it!=value.end(); ++it)
it->second/=M;
}

complex psi::operator()(const double x){

map<double,complex>::iterator it;
it=value.find(x);
if(it!=value.end())return it->second;
else{
for(it=value.begin();(it->first<x)&&it!=value.end();it++);
if(it!=value.begin()&&it!=value.end()){
double x1,x2;//będę tutaj aproksymował wartość pomiędzy
complex z1,z2;//dwoma wartościami
x2=it->first;
z2=it->second;
it--;
x1=it->first;//tutaj muszę założyć że mam dużo punktów
z1=it->second;//w małych odstępach między sobą
complex a = (z2-z1)/(x2-x1);//bardzo niebezpieczny
//kawałek - dzielenie małego przez małe.
complex z=a*(x-(x1+x2)/2)+(z1+z2)/2;//ten kawałek działa całkiem nieźle
return z;}
if(x<value.begin()->first){
complex A = 0;
return A;
}
else{it=value.end();it--;
complex A = 0;//jeśli wyszedłem poza zakres danych - zwraca zero
return A;}}}

void psi::prnt(void){
double L=p.L;
double dx=p.dx;
for(double x=-L;x<=L;x+=dx)
cout<<x<<"\t"<<val(x).real()<<"\t"<<val(x).imag()<<"\t"<<abs(val(x))*abs(val(x))<<endl;
}

void psi::prnt(double time){
double L=p.L;
double dx=p.dx;
for(double x=-L;x<=L;x+=dx)cout<<time<<"\t"<<x<<"\t"<<abs(val(x))<<endl;
cout<<endl;
}

#endif
