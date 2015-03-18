#include "soliton.cpp"


int main(int argc, char** argv){

if(argc>1)t2=atoi(argv[1]);
else t2=0;
sigma = 2;
//gpid=getpid();
gpid=0;
//int a=fork();
//int b=fork();
//int c=fork();
//ffrk+=(a==0)+2*(b==0)+4*(c==0);
p0.time=1200;
p0.L=500;
if(argc>2)p0.time=atoi(argv[2]);
p0.T=(int)(ffrk-'a')*5;
if(p0.T==0)p0.T=1;
ffal psi;
psi.evolve();
}

