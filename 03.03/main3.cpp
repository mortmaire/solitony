#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;
int main(int argc,char** argv){
fstream file(argv[1],ios::in);
double x;
while(file>>x)cout<<x<<endl;
}
