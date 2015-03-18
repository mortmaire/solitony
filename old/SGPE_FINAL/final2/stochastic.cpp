#include "psi.h"

//program obliczający ewolucję czasową rozwiązania dla stochastycznego równania Grossa Pitajewskiego (SGPE)



int main(int argc, char ** argv){

parameters p;

p.dx=1e-1;
p.L=40;
p.g=0.1;
p.gamma=0.01;
p.mu=20;
p.dt=0.001;


cerr<<"Program do tworzenia ewolucji czasowej rozwiązania stochastycznego równania Grossa-Pitajewskiego\n";

cerr<<"Podaj nazwę pliku .avi: ";
string s;cin>>s;

p.fill(argc,argv);

system("mkdir tmp");

p.I=p.t/p.dt/100;
p.t=p.dt*100;
psi fun(p,0); //tworzy stan próżni kwantowej
for(int i=0;i<p.I;i++){

cerr<<i*p.dt*100<<endl;
fun.fourier_iterate_stochastic();

stringstream ss2;ss2<<s<<".txt";
fstream out;
out.open(ss2.str().c_str(),fstream::out|fstream::trunc);

out<<fun;
out.close();
stringstream ss3;
ss3<<"gnuplot -e \"data='"<<ss2.str()<<"'; out='";
ss3<<s<<".png"<<"'\" ./scripts/stochastic.gp";

system(ss3.str().c_str());
stringstream ss;

ss<<"mv "<<s<<".png ./tmp/"<<setw(6)<<setfill('0')<<i<<".png";

system(ss.str().c_str());
}
cerr<<"Tworzenie pliku video..."<<endl;
stringstream ss;
ss<<"mencoder \"mf://tmp/*.png\" -mf fps=10 -o "<<s;
ss<<" -ovc lavc -lavcopts vcodec=mpeg4";
system(ss.str().c_str());
system("rm -r tmp");
}
