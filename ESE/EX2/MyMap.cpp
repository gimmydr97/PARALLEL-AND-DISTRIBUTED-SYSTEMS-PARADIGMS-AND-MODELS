//Assigment 2 (March 16): provide a 
//vector<float> map(vector<float>v, int mode, function<float(float)> f, int pardegree) 
//computing in parallel (using C++ threads or asyncs) the map(f) over v.
// Mode is 0 (block) or 1 (cyclic) scheduling. 

#include<iostream>
#include<thread>
#include<vector>
#include<functional>
#include <chrono>
#include <ctime>


using namespace std;

#define INITTIME \
  auto start = std::chrono::high_resolution_clock::now();\
  auto elapsed =std::chrono::high_resolution_clock::now()-start;\
  auto usec = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

#define BEGINTIME start = std::chrono::high_resolution_clock::now();

#define ENDTIME(s, nw) elapsed =std::chrono::high_resolution_clock::now()-start; \
  usec = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();\
  cout<<s<<"\t"<<usec<<" usec with "<<nw<<" threads"<<endl;
 
typedef enum{ BLOCK = 0, CYCLIC = 1} scatter_t;

vector<float> map(vector<float>v, scatter_t mode, function<float(float)> f, int pardegree){
    
    int dim = v.size();
    vector<float> result(dim);
    vector<thread> threads;

    if( mode == BLOCK){

        int blockDim = dim/pardegree;

        auto compute_block = [&](pair<int,int> range){
            for(int i = range.first; i < range.second; i++ ){
                result[i] = f(v[i]);
            }
            return;    
        };

        vector<pair<int,int>> pairs(pardegree);
        

        for(int i = 0; i < pardegree; i++){
            pairs[i] = make_pair(i*blockDim, ( i!= (pardegree-1) ? (i+1)*blockDim : dim));
        }

        for(int i = 0; i < pardegree; i++ ){
            threads.push_back(thread(compute_block, pairs[i]));
        }

        for(int i = 0; i < pardegree; i++ ){
            threads[i].join();
        }

        return result;
    
    }else{
        auto compute_cyclic = [&](int tno){
            for(int i = tno; i < dim; i+=pardegree ){
                result[i] = f(v[i]);
            }
            return;    
        };

        for(int i = 0; i < pardegree; i++ ){
            threads.push_back(thread(compute_cyclic, i));
        }

        for(int i = 0; i < pardegree; i++ ){
            threads[i].join();
        }

        return result;
    }


}

void activewait(std::chrono::milliseconds ms) {
  long msecs = ms.count();
  auto start = std::chrono::high_resolution_clock::now();
  auto end   = false;
  while(!end) {
    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
    if(msec>msecs)
      end = true;
  }
  return;
}

int f(int x) {
  activewait(1ms);
  return(x+x);
}
//g++ -pthread -O3 -std=c++17 MyMap.cpp -o es2
//./es2 24 1000  8 BLOCK

//se la grandezza dei numeri nel vettore è distribuita in modo bilanciato [1,50,2,70,25,89,17,63]
//allora il CHUNCK ha una miglior resa altrimenti se avremo un array sbilanciato del tipo [1,2,3,4,70,69,80,100]
//CYCLIC ha una miglior performance (naturalemnte perchè nel secondo caso gli ultimi thread 
//si troveranno a sostenere un carico di lavoro molto maggiore rispetto ai primi)
int main(int argc, char * argv[]) {

  if(argc == 1) {
    cout << "Usage is:" << argv[0] << " seed n pardegree [BLOCK|CYCLIC]" << endl;
    return(0);
  }  

  int seed = atoi(argv[1]);
  int n = atoi(argv[2]);
  int nw = atoi(argv[3]);
  scatter_t scatterMode = ((argv[4][0]=='B') ? BLOCK : CYCLIC); 

  const int max = 1024;
  srand(seed);
  vector<float> v(n),r1(n),r2(n);
  for(int i=0; i<n; i++)
    v[i] = rand()%max;
  

  // sequential execution
  {
    
    INITTIME
    BEGINTIME
    for(int i = 0; i < n; i++){
        r1[i] = f(v[i]);
    }
    ENDTIME("seq", 1)
  }

  // pthread parallel execution
  {                                       
    INITTIME
    BEGINTIME
    r2 = map(v, scatterMode, f, nw);
    ENDTIME("par", nw)
  }

  // merge is empty
  
  for(int i=0; i<n; i++)
    if(r1[i] != r2[i]) 
      cout << r1[i] << " != " << r2[i] << " at " << i << endl;

  return(0);
}
