#include "psi.h"

int main(int argc, char ** argv){


cerr<<"Podaj nazwę pliku z wynikami: ";
string wyniki;cin>>wyniki;



srand(time(0));

parameters p;

p.dx=1e-1;
p.L=40;
p.t=60;

p.g=0.01;
p.t1=100;
p.dt=0.01;
p.gamma=0.01;
p.mu=20;
system("mkdir baza 2> /dev/null");
int a,b,c;
a=fork();
b=fork();
c=fork();

int z=0;
if(a)z+=1;
if(b)z+=2;
if(c)z+=4;

int i=0;


while(p.gamma<0.4&&p.mu<500&&p.g>1e-5){


if(i==z){p.dt=0.05;
psi funkcja(p);
for(double x=-p.L;x<p.L;x+=p.dx)funkcja.add(x,0);
double entr;

do{
int odwr=1/p.dt;
while(odwr>2)odwr/=10;
if(odwr==2)p.dt/=5;
else p.dt/=2;
funkcja.p.dt=p.dt;
funkcja.fourier_iterate_stochastic();
entr=funkcja.entropy();
}
while(entr>0.1);
funkcja.prnt(funkcja.name());

stringstream s5;s5<<"try"<<z;
funkcja.prnt(s5.str());
stringstream s6;s6<<"gnuplot ./scripts/try"<<z<<".gp";
system(s6.str().c_str());
stringstream s1;s1<<"cp try"<<z<<".png "<<funkcja.name()<<".png";
stringstream s3;s3<<"rm try"<<z<<" try"<<z<<".png";
system(s3.str().c_str());
system(s1.str().c_str());
double t5=funkcja.hybrid();
stringstream ss1;
ss1<<p.gamma<<ta<<p.g<<ta<<p.mu<<ta<<t5<<endl;
cerr<<ss1.str();
if(t5==p.t1){

stringstream nazwa;
nazwa<<funkcja.name()<<".calc"<<p.dt<<".txt";
funkcja.prnt(nazwa.str());


}
stringstream outp;
outp<<"echo '";
outp<<p.gamma<<ta<<p.g<<ta<<p.mu<<ta<<t5;
outp<<"' >> "<<wyniki;
system(outp.str().c_str());
}


if(p.gamma<=0.099)p.gamma+=0.01;
else p.gamma+=0.1;
if(p.gamma==0.4){
p.gamma=0.01;
if(p.mu<100)p.mu+=10;
else p.mu+=100;
if(p.mu==300){
p.mu=20;
int odwr=1/p.g;
while(odwr>5)odwr/=10;
if(odwr==2){p.g/=5;p.g*=2;}
else p.g/=2;
}}
i++;i%=8;
}
}
