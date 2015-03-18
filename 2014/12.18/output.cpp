#include "activator.cpp"


ostream& operator<<(ostream& out, ffal& ket){
double L=ket.p.L;
double dx=ket.p.dx;
int n=2*L+1;
complex i0(0,1);
n/=dx;
out<<fixed<<setprecision(6);

stringstream bufor;
int tab[100];for(int i=0;i<100;i++)tab[i]=0;
double mu=0;double s2=0;
for(int i=0;i<n;i++){
double rad=abs(ket.tab[i]);
if(rad<ket.p.max)tab[(int)(rad*100/ket.p.max)]++;
}
int x_max,max;x_max=0;max=0;
for(int i=0;i<100;i++)if(max<tab[i]){x_max=i;max=tab[i];}
for(int i=0;i<100;i++)mu+=tab[i]/4001.*i;
for(int i=0;i<100;i++)s2+=tab[i]/4001.*(mu-i)*(mu-i);
s2=sqrt(s2);
//cerr<<ket.t<<"\t"<<mu<<"\t"<<s2<<endl;


activator *act[4][5];

for(int j=0;j<5;j++)for(int i=0;i<4;i++){
double thr=(mu-(2+i/2.*s2))*ket.p.max/100.;
act[i][j]=new activator(&bufor,thr,j);
}
int i=0;
for(double x=-L;x<=L;x+=dx){
double rad=abs(ket.tab[i]);
bufor<<"\n";
bufor<<ket.t<<"\t"<<x<<"\t"<<rad<<"\t"<<arg(ket.tab[i++]);
for(int j=0;j<5;j++)for(int i=0;i<4;i++)act[i][j]->run(rad);
}
for(int j=0;j<5;j++)for(int i=0;i<4;i++)act[i][j]->how_many();
bufor<<endl;
out<<bufor.str();
return out;
}

