#include <iostream>
#include <cstdlib>
#include <fstream>
using namespace std;

int number(int i=0,int=1,int=8192);
int main(){

int y=0;
fstream file("ciastko",ios::in);
double x;
while(file>>y){
string buffer;
getline(file,buffer,'\n');
//cout<<x<<endl;


}

file.close();

y=0;
for(int i=0;i<5;i++){
//file.open("ciastko",ios::out|ios::ate|ios::app);
y=number(y);
cout<<y<<endl;
//file.close();
}

}

int number(int i,int r_min,int r_max){
int r_size=r_max-r_min+1;
int step=r_size/2;
if(i==0)return r_min;
if(i==r_min)return r_max;
bool f=0;
if(i==r_max)f=1;
while(step>=1){
int x=step;
while(x<r_max){
if(f)return x;
if(x==i)f=1;
x+=2*step;

}
step/=2;
}
}
