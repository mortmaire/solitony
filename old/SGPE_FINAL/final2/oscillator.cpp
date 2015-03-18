#include "psi.h"

//program obliczający funkcję falową dla podstawowego oscylatora
//kwantowego



int main(int argc, char ** argv){

parameters p;

p.dx=1e-1;
p.L=40;
p.g=0.1;
p.gamma=0;
p.mu=22.41;
p.tau=0.001;
p.n=0;
p.N=1;

if(argc<1){
cerr<<"Program do tworzenia rozwiązań oscylatora harmonicznego oraz "<<"elementarnych rozwiązań równania Grossa-Pitajewskiego.\n";

return 0;
}
p.fill(argc,argv);

psi fun(f,p);//wrzucam funkcję trygonometryczną na start
cerr<<"Energia funkcji: "<<fun.iterate()<<endl;
cout<<fun;

}
