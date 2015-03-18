class aktywator{
int z,n;
double thr;
stringstream *bufor;
int total;
public:
aktywator(stringstream *bufor,double thr,int n):bufor(bufor),thr(thr),n(n),z(0),total(0){}
void run(double &val){
if(val<thr&&z<n){z++;(*bufor)<<" "<<1;}
else if(z>=n){
}
};


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
int i=0;

int s_count[8];
bool flag[8];
bool flag2[8];
double thr[8];
for(int i=0;i<1;i++){s_count[i]=0;flag[i]=0;flag2[i]=0;
thr[i]=(mu-s2)*ket.p.max/100.;
//cout<<thr[i]<<endl;
}

for(double x=-L;x<=L;x+=dx){
double rad=abs(ket.tab[i]);
bufor<<"\n";
bufor<<ket.t<<"\t"<<x<<"\t"<<rad<<"\t"<<arg(ket.tab[i++]);
for(int i=0;i<1;i++){


if(flag[i]==1&&flag2[i]==0){flag2[i]=1;s_count[i]++;}
if(rad<thr[i]&&flag[i]==0)flag[i]=1;
if(rad>thr[i]&&flag[i]==1){flag[i]=0;flag2[i]=0;}
bufor<<" "<< ((flag2[i])?0:1);
}
}
for(int i=0;i<1;i++)bufor<<" "<<s_count[i];
bufor<<endl;
out<<bufor.str();
return out;
}

