#ifndef FUNCT_H
#define FUNCT_H


#define complex complex<double> //to jedyny rodzaj liczb zespolonych z jakich korzystam
#define _(x) #x

void replaceAll(string& str, const string& from, const string& to) {
    if(from.empty())
        return;
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}

void system(stringstream &ss){
string s=ss.str();
system(s.c_str());
}

int fact(int);//działa do n=16!
complex operator/(complex,double);
complex operator*(double,complex);
complex operator*(int,complex);
int knmk(int,int);
double an(double x,double L);
double f(double,double);
int n00;//globalna zmienna dla funkcji analitycznej. Bóg mi wybaczy.
double mod=1;

bool operator==(complex& c1, complex& c2){
double x1,x2,y1,y2;
x1=c1.real();
x2=c2.real();
y1=c1.imag();
y2=c2.imag();
if((x1==x2)&&(y1==y2))return 1;
else return 0;
}


double an(double x,double L){//analityczne rozwiązanie osc. harm.
int n=n00;
//cerr<<n<<endl;
x+=n*L;
x/=n+1;
//x+=n*L/(n+1);
double p=exp(-x*x/2);
p/=sqrt(pow(2,n)*fact(n))*pow(M_PI,0.25);
double h=0;

for(int m=0;m<=n/2;m++){
h+=pow(-1,m)/(fact(m)*fact(n-2*m))*pow(2*x,n-2*m);
}
h*=fact(n);
//h=1;
return mod*p*h;
}

double f(double x,double L){
return cos(M_PI*(x)/(2*L));
//return sin(M_PI/L*(x+L));
}

int fact(int n){
if(n<1)return 1;
int m=1;
for(int i=1;i<=n;i++)m*=i;
return m;
}

complex operator/(complex z0,double n){
double re=z0.real()/n;
double im=z0.imag()/n;
complex z(re,im);
return z;
}

complex operator*(double a,complex z0){
double re=z0.real()*a;
double im=z0.imag()*a;
complex z(re,im);
return z;
}

complex operator*(int a,complex z0){
double re=z0.real()*a;
double im=z0.imag()*a;
complex z(re,im);
return z;
}

int knmk(int k,int n){
if(k>n){int a=k;k=n;n=a;}
int a=1;
for(int i=1;i<=k;i++)a*=i;
for(int i=1;i<=n-k;i++)a*=i;
return a;
}

#endif
