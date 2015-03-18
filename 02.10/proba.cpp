#include <iostream>
#include <ctime>
#include <cstdlib>
using namespace std;

void ciastko(int a){
static int i=0;
if(a>0)i++;
cout<<i<<endl;

}

int main(){
/*
int **p;
int N=1e4;
p=new int*[N];
for(int i=0;i<N;i++)p[i]=new int[N];

double t;

t=time(0);
for(int I=0;I<10;I++)
for(int j=0;j<N;j++)for(int i=0;i<N;i++)p[i][j]=4;
t=time(0)-t;
cout<<t<<endl;

t=time(0);
for(int I=0;I<10;I++)
for(int i=0;i<N;i++)for(int j=0;j<N;j++)p[i][j]=2;
t=time(0)-t;
cout<<t<<endl;

for(int i=0;i<N;i++)delete[] p[i];
*/
ciastko(0);
ciastko(0);
ciastko(1);
ciastko(1);
ciastko(0);

}
