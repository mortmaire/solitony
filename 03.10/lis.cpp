#include <iostream>
using namespace std;

int main(){
int r_min=1;
int r_max=8192;
int r_size=r_max-r_min+1;
int step=r_size/2;
int i=1;
cout<<r_min<<"\t"<<i<<endl;
cout<<r_max<<"\t"<<i<<endl;
while(step>=1){
int x=step;
while(x<r_max){
if(i++>100)return 0;
cout<<x<<"\t"<<1<<endl;
x+=2*step;
}
step/=2;
}
}
