#include <iostream>

using namespace std;

int main(){
double sum=0;
double sum2=0;
double buff[4];
double s1,s2;
s1=10;s2=10;
int i=0;
double tab[4001];
while(cin>>buff[0]>>buff[1]>>buff[2]>>buff[3]){
sum+=buff[2];
sum2+=buff[2]*buff[2];
if(s1>buff[2])s1=buff[2];
if(s2<buff[2])s2=buff[2];
tab[i++]=buff[2];
}
sum/=i;
sum2/=i;
cerr<<sum<<" "<<sum2<<" "<<s1<<" "<<s2<<" "<<i<<endl;
double szer=s2-s1;
szer/=100;
int rnorm[100];
for(int ii=0;ii<100;ii++)rnorm[ii]=0;
for(int ii=0;ii<i;ii++){
rnorm[(int)(tab[ii]/szer)]++;
}
for(int ii=0;ii<100;ii++)cout<<szer*ii<<"\t"<<rnorm[ii]<<endl;
}

