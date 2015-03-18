#include "soliton.h"

int main(int argc, char** argv){
//if(argc>1)t2=atoi(argv[1]);
else t2=0;
gpid=getpid();
//s_mod=2.5;

int a=fork();
int b=fork();
int c=fork();
ffrk+=(a==0)+2*(b==0)+4*(c==0);

if(ffrk-'a'==8)return 0;
p0.time=2000;
cout<<p0.time<<endl;
//if(argc>2)p0.time=atoi(argv[2]);
p0.T=(int)(ffrk-'a')+1;
ffal psi(p0);
psi.evolve();
}

