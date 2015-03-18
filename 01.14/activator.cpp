using namespace std;
class activator{
int z,n;
double thr;
stringstream *bufor;
int total;
public:
activator(stringstream *bufor,double thr=2.5,int n=2):bufor(bufor),thr(thr),n(n),z(0),total(0){}

void set_threshold(double thr2){this->thr=thr2;}
void run(double val){
if(val<thr&&z<n){z++;(*bufor)<<" "<<1;}
else if(z==n&&val<thr){total++;(*bufor)<<" "<<0;z++;}
else if(z>n&&val<thr){(*bufor)<<" "<<0;}//nie musisz zwiększać z
else if(val>=thr){z=0;(*bufor)<<" "<<1;}
}
void how_many(void){
(*bufor)<<" "<<total;
}
};
