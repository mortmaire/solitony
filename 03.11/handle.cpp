#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
using namespace std;

class funkcja{
double *tab;
int n;
double min;
double max;
double delta;
public:
funkcja(double* tab,int n,double min,double max):tab(tab),n(n),min(min),max(max){
if(n>1)
delta=(max-min)/(n-1);
else n=1;
}
double operator()(double val){
if(val==min)return tab[0];
if(val==max)return tab[n-1];
if(val>max)return 0;
if(val<min)return 0;
int w=(val-min)/delta;
double p1=(val-min)/delta-w;
double p2=(w+1)+(min-val)/delta;
return p1*tab[w+1]+p2*tab[w];
}
double norm(){
double calka=0;
for(int i=0;i<n;i++){
calka+=tab[i];
}
for(int i=0;i<n;i++)tab[i]/=calka;
return calka;
}
double norm(int S){
for(int i=0;i<n;i++)tab[i]/=S;

}
};


int main(){



int n=0;
for(int i=1;i<8193;i++){
stringstream ss;
ss<<"../ss/T"<<setw(4)<<setfill('0')<<i<<".a";

fstream file(ss.str().c_str(),ios::in);
if(file.is_open())n++;
file.close();
}

int tab[2][n];
int j=0;
double max=0;
for(int i=1;i<8193;i++){
stringstream ss;
ss<<"../ss/T"<<setw(4)<<setfill('0')<<i<<".a";
fstream file(ss.str().c_str(),ios::in);
if(file.is_open()){
tab[0][j++]=i;
string bufor;
stringstream bb;
double x,y;
while(file>>x>>y){
if(max<x)max=x;
}
}
file.close();
}



int N=1000;




for(int i=0;i<n;i++){

stringstream ss;
int T=tab[0][i];
ss<<"../ss/T"<<setw(4)<<setfill('0')<<T<<".a";
fstream file(ss.str().c_str(),ios::in);
string bufor;
double tab2[2][100];
int j=0;

while(getline(file,bufor)){
stringstream bb;
bb<<bufor;
bb>>tab2[0][j]>>tab2[1][j];
j++;
}

//tutaj wrzucić normowanie tablicy

funkcja f(tab2[1],100,tab2[0][0],tab2[0][99]);
tab[1][i]=f.norm();//nie jestem pewien czy moja norma jest odpowiednia
//może lepiej wrzucić norm do max val = 1
for(int k=0;k<N;k++){
cout<<tab[0][i]<<" "<<k/(double)N*max<<" "<<f(k/(double)N*max)<<endl;
}
cout<<endl;
cerr<<tab[0][i]<<endl;

}
cout<<endl;


for(int i=0;i<n;i++){
stringstream ss;
int T=tab[0][i];
ss<<"../ss/T"<<setw(4)<<setfill('0')<<T<<".b";
fstream file(ss.str().c_str(),ios::in);
string bufor;
double tab2[2][100];
int j=0;
while(getline(file,bufor)){
stringstream bb;
bb<<bufor;
bb>>tab2[0][j]>>tab2[1][j];
j++;
}

//tutaj wrzucić normowanie tablicy

funkcja f(tab2[1],100,tab2[0][0],tab2[0][99]);
f.norm(tab[1][i]);
for(int k=0;k<N;k++){
cout<<tab[0][i]<<" "<<k/(double)N*max<<" "<<
f(k/(double)N*max)<<endl;
}
cout<<endl;
cerr<<tab[0][i]<<endl;


}




}
