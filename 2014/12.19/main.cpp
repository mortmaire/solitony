#include "soliton.cpp"

int main(int argc, char** argv){
if(argc>1)t2=atoi(argv[1]);
else t2=0;
gpid=time(0);
sigma = 2;
int a=fork();
int b=fork();
int c=fork();
ffrk+=(a==0)+2*(b==0)+4*(c==0);

if(ffrk-'a'==8)return 0;
p0.time=1200;
if(argc>2)p0.time=atoi(argv[2]);
for(int i=0;i<100;i++){
p0.T=(int)(ffrk-'a')+8*i+20;
ffal psi;
psi.evolve();
gpid=time(0);}
}

