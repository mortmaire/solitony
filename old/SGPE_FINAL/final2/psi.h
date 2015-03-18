#ifndef PSI_H
#define PSI_H
#define ta "\t"




//131649
//131427
//131324
//131218
using namespace std;

class psi;
#include <iostream>
#include <cmath>
#include <complex>
#include <map>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <fftw3.h>
#include "funct.h"
#include <ctime>
#include <cstdlib>
#define check(x) ( cout<<endl<<#x )

class parameters{
public:
double dx;//krok iteracyjny
double L;//szerokość sandboxu (2xL)
int n; //numer stanu wzbudzonego
int N; //ilość cząstek
double g; //coupling constant
double tau; //parametr czasu urojonego
int I; //ilość iteracji
double t; //czas wykonywania SGPE
double dt; //krok czasowy
double gamma; //parametr gamma
double T; //temperatura
double mu; //potencjał chemiczny
double t1; //czas wykonywania równań hybrydowych
parameters(int n=0, double dx=1e-1, double L=40,int I=10000,double tau=1e-6,double g=0.1,double N=1,double t=60,double dt=0.0001,double gamma=0.01,double T=13.89,double mu=22.41,double t1=60):
n(n),L(L),dx(dx),N(N),g(g),tau(tau),I(I),t(t),dt(dt),gamma(gamma),T(T),mu(mu),t1(t1){}

void show(void);
bool variable(string in, double &var, string str);
bool variable(string in, int &var, string str);
bool fill(int argc, char** argv);
};



class psi{
map<double,complex> value;//wartości funkcji psi
map<double,complex> value_c; //sprzęgnięta funkcja psi
public:
void add_conj(complex*);
psi create_conjugate(void);//funkcja tworząca sprzężoną funkcję psi
parameters p;
int show_NL(void){return 2*p.L/p.dx;}
void add(double x, complex f);//podstawowa funkcja dodawania nowych wartości. muszę stworzyć op.
double show_L(void){return p.L;}
double show_dx(void){return p.dx;}
int show_n(void){return p.n;}
psi(double (*f)(double,double),parameters p);//tworzenie nowego elementu klasy
psi(parameters p):p(p){}//tworzenie pustego elementu klasy
psi(string name, parameters p);
psi(void){};
psi(parameters p, double); //tworzenie bazowej ffal, domyślnie stanu próżni
complex operator()(const double x);//operator zwracania wartości funkcji
complex val(const double x){return operator()(x);}
complex scalar_product(psi &ket);//iloczyn skalarny
complex operator*(psi &ket){return scalar_product(ket);}//op. iloczynu skalarnego
void norm(void);//normalizacja.
double energy_product(psi &ket2);//zwraca <psi_1|H|psi_2>
double operator%(psi &ket2){return energy_product(ket2);}//operator średniej energii

void prnt(void);//drukuje wszystkie wartości funkcji w formie (x,f(x)).
void prnt(string name);

void prnt(double);

void shift(complex*,int=-1);//zmienia tablicę na pierwszym argumencie w funkcję w przestrzeni pędów
void fourier_iterate_stochastic();
psi(complex*,parameters p);//tworzenie ffal z tablicy liczb
double hybrid();
double hybrid(int);
bool hybrid(psi&);//wariacja na temat wcześniejszych hybryd
complex eta(void);
friend ostream& operator<<(ostream& os, psi& ket);

bool ifex();
string name();
double entropy(void);
double iterate(void);
};



#include "psi_body.h"
#endif
