#include <iostream>
#include <fftw3.h>
#include <cstdlib>
#include <complex>
using namespace std;
#define complex complex<double>
inline complex cmpx(double *tab){
return complex(tab[0],tab[1]);
}

void czytacz(double *tab[2]){
for(int i=0;i<10;i++)cout<<tab[i][0]<<tab[i][1]<<endl;
}


int main(){

int n=1024;

double **tab=new double*[n];
for(int i=0;i<n;i++){
tab[i]=new double[2];
tab[i][0]=i;
tab[i][1]=0;
}
czytacz(tab);

}
