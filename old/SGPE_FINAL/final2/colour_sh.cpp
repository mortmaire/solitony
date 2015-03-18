#include "psi.h"

//program obliczający ewolucję czasową rozwiązania dla stochastycznego równania Grossa Pitajewskiego (SGPE)



int main(int argc, char ** argv){

parameters p;

p.dx=1e-1;
p.L=40;
p.g=0.001;
p.gamma=0.1;
p.mu=50;
p.dt=0.0001;

if(argc<1){
cerr<<"Program generujący kolorowy przebieg ewolucji układu SGPE i \n";
return 0;
}
fstream tem;

p.fill(argc,argv);

cerr<<"Podaj nazwę wyjściowych .txt: ";
string name;cin>>name;

stringstream s1,s2;s1<<name<<1;
s2<<name<<2;
string name1=s1.str();
string name2=s2.str();

cerr<<"Podaj nazwę wyjściowego .png: ";
string n2;cin>>n2;

fstream out1, out2;
out1.open(name1.c_str(),fstream::out|fstream::trunc);
out2.open(name2.c_str(),fstream::out|fstream::trunc);
p.I=p.t/0.1;
p.t=0.1;
cerr<<p.I<<endl;
psi fun(p,0); //tworzy stan próżni kwantowej
for(int i=0;i<p.I;i++){
cerr<<i*p.t<<endl;
fun.fourier_iterate_stochastic();
tem.open("temp2.txt",fstream::out|fstream::trunc);
tem<<fun;
tem.close();
string ciastko;
tem.open("temp2.txt",fstream::in);
while(getline(tem,ciastko)){
out1<<p.t*i<<ta<<ciastko<<endl;
out2<<p.t*i<<ta<<ciastko<<endl;}
out1<<endl;
out2<<endl;
tem.close();
}
int I0=p.I;
p.I=p.t1/0.1;
p.t1=0.1;

psi f1,f2;
fun.p=p;
f1=fun;
f2=fun;

for(int i=0;i<p.I;i++){
cerr<<(i+I0)*p.t1<<endl;
bool f=f1.hybrid(f2);
if(f==0)break;


tem.open("temp2.txt",fstream::out|fstream::trunc);
tem<<f1;
tem.close();
string ciastko;
tem.open("temp2.txt",fstream::in);
while(getline(tem,ciastko))out1<<p.t1*(i+I0)<<ta<<ciastko<<endl;
out1<<endl;
tem.close();

tem.open("temp2.txt",fstream::out|fstream::trunc);
tem<<f2;
tem.close();
ciastko;
tem.open("temp2.txt",fstream::in);
while(getline(tem,ciastko))out2<<p.t1*(i+I0)<<ta<<ciastko<<endl;
out2<<endl;
tem.close();



}


out1.close();
out2.close();
system("rm temp2.txt");


stringstream ss;
ss<<"gnuplot -e \"data1='"<<name1<<"'; data2='"<<name2<<"'; out='";
ss<<n2<<"'\" ./scripts/colour_sh.gp";

system(ss.str().c_str());

}
