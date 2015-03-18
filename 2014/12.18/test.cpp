#include <iostream>
#include <sstream>
#include "activator.cpp"
#include <ctime>
#include <chrono>
int main(){
stringstream bufor;
//activator akt(&bufor,1,0);
//for(int i=0;i<5;i++)akt.run(0);
//cout<<bufor.str().c_str()<<endl;
time_t result = time(NULL);
cout << ctime(&result);

}
