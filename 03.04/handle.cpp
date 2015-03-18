#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
using namespace std;
int main(){
int n=0;
for(int i=1;i<8193;i++){
stringstream ss;
ss<<"../ss/T"<<setw(4)<<setfill('0')<<i<<".a";

fstream file(ss.str().c_str(),ios::in);
if(file.is_open())n++;
file.close();
}
int tab[n];
int j=0;
double max=0;
for(int i=1;i<8193;i++){
stringstream ss;
ss<<"../ss/T"<<setw(4)<<setfill('0')<<i<<".a";
fstream file(ss.str().c_str(),ios::in);
if(file.is_open()){
tab[j++]=i;
string bufor;
stringstream bb;
double x,y;
while(file>>x>>y){
if(max<x)max=x;
}
}
file.close();
}
cout<<max<<endl;
double tab3[n][1000];
for(int i=0;i<n;i++){
stringstream ss;
int T=tab[i];
ss<<"../ss/T"<<setw(4)<<setfill('0')<<T<<".a";
fstream file(ss.str().c_str(),ios::in);
string bufor;
double tab2[100][2];
int j=0;
stringstream bb;
while(getline(file,bufor)){
bb<<bufor;
bb>>tab2[j][0]>>tab2[j++][1];
}
double minn=tab2[0][0];
int maxx=tab2[100][0];
for(int j=0;j<1000;j++){
double val=j*max/1000;
if(val<minn)tab[i][j]=0;
else if(val>maxx)tab[i][j]=0;
else{
if(val==minn)tab[i][j]=minn;
else if(val==maxx)tab[i][j]=maxx;
else tab[i][j]=//muszę tutaj walnąć funkcję interpolacyjną
}

}
}

}
