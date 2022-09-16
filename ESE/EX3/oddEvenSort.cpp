/*Assigment 3: implement in parallel the odd even sort.
  Use as input a random generated integer vector. 
  The odd even sort algorithm looks like the following (sequential version) 
  (n is the length of the vector, 
   v is the vector to sort, 
   swap is a procedure that exchanges the two positions in the vector):

while(not(sorted)) {
  sorted = true;
  for(i=1; i<n-1; i+=2)
    if(v[i] > v[i+1]) { swap(v[i], v[i+1]); sorted = false;}
  for(i = 0; i<n-1; i+=2)
    if(v[i] > v[i+1]) { swap(v[i], v[i+1]); sorted = false;}
}

with the intent to discuss scalability of the code.*/
#include<iostream>
#include<vector>
#include<string>
#include<thread>
#include<atomic>
#include<barrier>

#include "utimer.cpp"

using namespace std;

//initialize a vector with random numbers in
void initV(vector<int>& v, int max, int seed){
  srand(seed);
  for(auto &e: v)
    e = rand() % max;
  return;
}

//show vector
void printV(string msg, vector<int> v){
  cout << msg << ": ";
  for(auto &e : v)
    cout << e << " ";
  cout << endl;
  return;
}

//sequential oesort
int oesortSeq(vector<int>& v){
  bool sorted = false;
  int n = v.size();
  int iter = 0;
  while(!sorted){
    iter++;
    sorted = true;
    //odd phase
    for(int i=1; i <n-2; i+=2){
      if(v[i]>v[i+1]){
          auto temp = v[i];
          v[i] = v[i+1];
          v[i+1] = temp;
          sorted = false;
      }
    }
    //even phase
    for(int i=0; i <n-1; i+=2){
      if(v[i]>v[i+1]){
          auto temp = v[i];
          v[i] = v[i+1];
          v[i+1] = temp;
          sorted = false;
      }
    }
  }
  return iter;
}

//THR oesort
int oesortParallel(vector<int>& v, int nw){
  atomic<bool> sorted=false;
  //parallelism degree, callback
  barrier ba(nw, [&](){sorted = true; return;});
  barrier bb(nw, [&](){return;});

  vector<pair<int,int>> chunks(nw);
  auto n = v.size();
  auto d = n/nw; // assume is multiple (e.g. n= 2^k)
  for(int i = 0; i<nw; i++){
    auto start = i*d;
    auto stop =(i==(n-1) ? n : d*(i+1));
    chunks[i] = make_pair(start, stop); 
  }
  atomic<int> globaliters;
  globaliters = 0;

  auto body = [&] (int ti){
                auto start = chunks[ti].first;
                auto stop = chunks[ti].second; 
                auto last = (ti == (nw-1));
                if(last) stop = chunks[ti].second-1;
                auto iters = 0;

                while(!sorted){
                  iters++;
                  //odd step
                  auto localsorted = true;
                  for(int i= start+1; i < stop; i+=2)
                    if(v[i]>v[i+1]){
                      auto temp = v[i];
                      v[i] = v[i+1];
                      v[i+1] = temp;
                      localsorted = false;
                    }
                  
                  //wait all
                  ba.arrive_and_wait();
                  //even step, start with global sorted true
                  for(int i= start; i < stop; i+=2)
                    if(v[i]>v[i+1]){
                      auto temp = v[i];
                      v[i] = v[i+1];
                      v[i+1] = temp;
                      localsorted = false;
                    }
                 
                  //wait all
                  if(localsorted == false) // first update global sorted
                      if(sorted){
                        sorted = false;
                      }
                  bb.arrive_and_wait();
                  //go to the next iteration
                }
                globaliters = iters;
                return;             
  };

  //now create threads
  vector<thread*> tids(nw);
  for(int i= 0; i<nw; i++){
    tids[i] = new thread(body,i);
  }
  //and await their termination
  for(int i=0; i<nw; i++){
    tids[i]->join();
  }

  return globaliters;
}

//g++-11 -std=c++20 -pthread -O3  oddEvenSort.cpp -o es3
//./es3 123 1024 65536(deve essere una potenza del 2 in questo caso 2^16) 0 2
int main(int argc, char* argv[]){
  //dbg = 1 : print vectors
  if(argc == 1){
    cout << "Usage is: " <<argv[0]<<" seed max n dgb [nw]" <<endl;
    return 0;
  }

  int seed = atoi(argv[1]); //seed for random number generator
  int max = atoi(argv[2]);  //max value in the vector
  int n = atoi(argv[3]);   //length of the vector
  int dbg = (atoi(argv[4]) == 0 ? false : true); //dbg flag
  int nw = (argc == 6 ? atoi(argv[5]) : 1); //par degree (default 1 )

  //initialize the vector

  vector<int> v(n);
  initV(v,max,seed);
  if(dbg) printV("Int", v);

  long t1, t2, iters;
  {
    utimer tseq("Seq", &t1);
    iters = oesortSeq(v);
    //cout << "Sorted in " << iters << " iterations " <<endl;

  }
  //cout << "Iteration took average  " << t1 / iters << " usecs " <<endl;

  if(dbg) printV("Seq sort", v);
  cout << "Seq oesort " << (is_sorted(v.begin(), v.end()) ? "works" : "error! ") <<endl;
  
  // re-initialize, same numers
  initV(v,max,seed);
  {
    utimer tthread("THR", &t2);
    iters = oesortParallel(v,nw);
    // cout << "Sorted in " << iters << "iterations " <<endl,
  }

  if(dbg) printV("THR sort", v);
  cout << "THR oesort " << (is_sorted(v.begin(), v.end()) ? "works" : "error! ") <<endl;

  cout << "THR speedup " << nw << " " << ((float) t1)/((float) t2) << endl;
  return 0;
}