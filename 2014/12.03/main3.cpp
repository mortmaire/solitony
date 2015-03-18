#include <iostream>
#include <fftw3.h>
#include <fstream>
#include "../bmp/EasyBMP.h"

using namespace std;

int main(){

BMP Image;
Image.ReadFromFile("1.bmp");
int W=Image.TellWidth();
int H=Image.TellHeight();
double tab[W];
for(int i=0;i<W;i++)tab[i]=0;//zerowanie

for(int i=0;i<W;i++)
for(int j=0;j<H;j++)
if(Image(i,j)->Red==0)tab[i]+=(double)j;

fftw_complex *in, *out;
fftw_plan p0;//tutaj przysłaniam globalny parametr p0
in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * W);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * W);

//for(int i=0;i<W;i++)cerr<<i<<"\t"<<tab[i]<<endl;

for(int i=0;i<W;i++){
in[i][0]=tab[i];
in[i][1]=0;
}

p0 = fftw_plan_dft_1d(W, in, out, 1, FFTW_ESTIMATE);
fftw_execute(p0);
for(int i=0;i<W;i++)tab[i]=out[i][0];

fftw_destroy_plan(p0);
fftw_free(in);
fftw_free(out);

for(int i=0;i<W;i++)cerr<<i<<"\t"<<tab[i]<<endl;




}
