#include <iostream>
#include <cstdlib>
#include <stdio.h>
#define READ(...) # __VA_ARGS__
using namespace std;


int main(){

char a[500];
sprintf(a,"%d%c",500,'a');
cout<<a<<endl;

}
