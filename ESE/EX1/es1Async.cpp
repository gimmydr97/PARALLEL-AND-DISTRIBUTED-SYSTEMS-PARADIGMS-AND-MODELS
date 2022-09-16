/*Assignment 1 (March 9): evaluate overhead needed to start/join a thread and 
  to start an async (launch::async policy) and get the result of the computed future.
  Measures should be taken as averages of a large number of events. 
  Possibly, measure times on your machine AND on the remote virtual machine.*/

#include<thread>
#include<atomic>
#include<chrono>
#include<vector>
#include<future>
#include<iostream>

using namespace std;

#define INITTIME \
  auto start = std::chrono::high_resolution_clock::now();\
  auto elapsed =std::chrono::high_resolution_clock::now()-start;\
  auto usec = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

#define BEGINTIME start = std::chrono::high_resolution_clock::now();

#define ENDTIME(s, nw) elapsed =std::chrono::high_resolution_clock::now()-start; \
  usec = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();\
  cout<<s<<"\t"<<usec<<" usec with "<<nw<<" threads"<<endl;

//g++ -pthread -O3 es1Async.cpp -o es1
// ./es1 (n =)1000 (nw =)16


int main(int argc, char** argv){

  int n = atoi(argv[1]);
  int nw = atoi(argv[2]);

  INITTIME

  vector<future<int>> x(nw);

  BEGINTIME

  for(int i = 0; i<n; i++){
    for(int j = 0; j<nw; j++){
      x[j]= async(launch::async,
                  [] (int i){return i;},
                  j);
    }
    int sum = 0;
    //faccio qualcosa per aspettare la fine delle async 
    //cercando di fare qualcosa (sum) che impedisca che il codice sia eliso durante la fase di miglioramento
    for(int j = 0; j<nw; j++){
      sum += x[j].get();
    }
  }
  ENDTIME("raw time", nw)

   cout << "Average per async (future + get) " 
        << ((((float) usec) / ((float) n)) / ((float) nw)) <<endl;
  return 0;
}