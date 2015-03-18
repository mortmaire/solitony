void evolve(){
fstream file;
fstream file2;
stringstream ss2;ss2<<"../temp/"<<gpid<<ffrk;
string name1=ss2.str()+".tmp.dat";
string name3=ss2.str()+"S.tmp.dat";

if(ffrk=='a')cerr<<name1<<endl;

file.open(name1.c_str(),ios::out);//czyszczę plik tmp.dat
file2.open(name3.c_str(),ios::out);
file.close();
file2.close();
int i=0;
int I=p.time/p.dt/1e3;
if(ffrk=='a')system("setterm -cursor off");
while(t<p.time){
potential();
kinetic2();
potential();
t+=2*p.dt;
i+=2;
if(i%(I/100)==0){cerr<<setfill('0')<<setw(7)<<fixed<<setprecision(2);if(ffrk=='a')cerr<<t<<"\r";}
if(i>I){i=0;
file.open(name1.c_str(),ios::out|ios::ate|ios::app);
file<<*this;
file.close();
file2.open(name3.c_str(),ios::out|ios::ate|ios::app);
file2>>*this;
file2.close();
}
}
if(ffrk=='a')system("setterm -cursor on ");
stringstream ss;
ss<<"T = "<<p.T;
char title[20];sprintf(title,"T = %.3f",p.T);
char gnuplot[1000];ffrk;

#define write(x) #x 

//sprintf(gnuplot,ss2.str().c_str(),0.,p.time,p.T,p.L);
string name2=ss2.str()+".gp";
fstream f2;
f2.open(name2.c_str(),ios::out);
f2<<gnuplot;
f2.close();
string systr="gnuplot "+name2;
system(systr.c_str());

int n=0;

}

