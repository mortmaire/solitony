#include "activator.cpp"


ostream& operator<<(ostream& out, ffal& ket){
double L=ket.p.L;
double dx=ket.p.dx;
int n=2*L/dx+1;
complex i0(0,1);
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

int n0=9;

activator *act[n0];
int j=0;
for(double ii=1;ii<5.5;ii+=0.5){
double thr=(mu-ii*s2)*ket.p.max/100.;
act[j++]=new activator(&bufor,thr,1);
}
int i=0;
for(double x=-L;x<=L;x+=dx){
double rad=abs(ket.tab[i]);
bufor<<"\n";
bufor<<ket.t<<"\t"<<x<<"\t"<<rad<<"\t"<<arg(ket.tab[i++]);
for(int i=0;i<9;i++)act[i]->run(rad);
}
for(int i=0;i<9;i++)act[i]->how_many();
bufor<<endl;
out<<bufor.str();
return out;
}

ostream& operator>>(ostream& out, ffal& ket){//to tylko inny output
double L=ket.p.L;
double dx=ket.p.dx;
int n=2*L/dx+1;
complex i0(0,1);
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
//double thr=(mu-ii*s2)*ket.p.max/100.;
//act[j++]=new activator(&bufor,thr,1);

//int i=0;
for(int xx=0;xx<100;xx++){
//double rad=abs(ket.tab[i]);
bufor<<"\n";
int y=2*x_max-xx;
if(y<0||y>=100)y=99;
bufor<<ket.t<<"\t"<<(xx)/100.*ket.p.max<<"\t"<<tab[xx]<<"\t"<<tab[xx]-tab[y];}
bufor<<"\t"<<mu<<"\t"<<s2<<endl;
out<<bufor.str();
return out;
}

