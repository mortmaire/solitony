#include "psi.h"

//program obliczający ewolucję czasową rozwiązania dla stochastycznego równania Grossa Pitajewskiego (SGPE)



int main(int argc, char ** argv){

parameters p;

p.dx=1e-1;
p.L=40;
p.g=0.1;
p.gamma=0.01;
p.mu=20;
p.dt=0.0001;

if(argc<1){
cerr<<"Program generujący kolorowy przebieg ewolucji układu SGPE\n";
return 0;
}
fstream tem;

p.fill(argc,argv);

cerr<<"Podaj nazwę wyjściowego .txt: ";
string name;cin>>name;
cerr<<"Podaj nazwę wyjściowego .png: ";
string n2;cin>>n2;
fstream out;out.open(name.c_str(),fstream::out|fstream::trunc);
p.I=p.t/0.1;
p.t=0.1;
cerr<<p.I<<endl;
psi fun(p,0); //tworzy stan próżni kwantowej
for(int i=0;i<p.I;i++){
cerr<<i*p.t<<endl;
fun.fourier_iterate_stochastic();
tem.open("temp.txt",fstream::out|fstream::trunc);
tem<<fun;
tem.close();
string ciastko;
tem.open("temp.txt",fstream::in);
while(getline(tem,ciastko))out<<p.t*i<<ta<<ciastko<<endl;
out<<endl;
tem.close();
}
out.close();
system("rm temp.txt");

stringstream ss;
ss<<"gnuplot -e \"data='"<<name<<"'; out='";
ss<<n2<<"'\" ./scripts/colour.gp";
system(ss.str().c_str());

}
