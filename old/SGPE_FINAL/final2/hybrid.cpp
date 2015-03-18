#include "psi.h"

//program obliczający ewolucję czasową rozwiązań dla równań hybrydowych



int main(int argc, char ** argv){

parameters p;

p.dx=1e-1;
p.L=40;
p.g=0.1;
p.gamma=0.01;
p.mu=20;
p.dt=0.0001;

if(argc<1){
cerr<<"Program do tworzenia ewolucji czasowej rozwiązań równań hybrydowych Grossa-Pitajewskiego\n";
cerr<<"USAGE: time / dt / gamma / g / mu / "<<endl;
return 0;
}
p.fill(argc,argv);

psi fun(p,0); //tworzy stan próżni kwantowej
fun.fourier_iterate_stochastic();
cerr<<fun.entropy()<<endl;
fun.hybrid();
psi fun_c=fun.create_conjugate();
fun.prnt("try1.txt");
fun_c.prnt("try2.txt");

}
