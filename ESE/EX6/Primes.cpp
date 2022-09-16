/*Assignment 6: write a parallel  program that finds all prime numbers 
  in a given range of values, passed as parameters of the command line. 
  In order to check if a number is prime, please use the following code:

// see http://en.wikipedia.org/wiki/Primality_test

static bool is_prime(int n) {
   if (n <= 3) return n > 1; // 1 is not prime !

   if (n % 2 == 0 || n % 3 == 0) return false;
   for (ull i = 5; i * i <= n; i += 6) {
        if (n % i == 0 || n % (i + 2) == 0)
            return false;
    }
   return true;
}

Consider your favorite parallel programming framework among the ones we've seen so far 
(C++ native threads, OpenMP, GrPPI) and possibly implement more that one version 
(e.g. a native C++ and an OpenMP version) comparing the relative performances.*/

#include <iostream>
#include "utimer.cpp"
#include <vector>
#include "ff/ff.hpp"
#include <ff/parallel_for.hpp>

using namespace std;
//using namespace ff;

static bool is_prime(int n) {
   if (n <= 3) return n > 1; // 1 is not prime !

   if (n % 2 == 0 || n % 3 == 0) return false;
   for (int i = 5; i * i <= n; i += 6) {
        if (n % i == 0 || n % (i + 2) == 0)
            return false;
    }
   return true;
}

//g++-11 -O4 -pthread -fopenmp Primes.cpp -o es6
//./es6 1 1000000 2(= pardeg)
//for((i=0;i<8;i++)); do ./es6 1 1000000 4 ; done | grep OMP
//OMP is a little bit faster than ff
int main(int argc, char * argv[]){

    int from = atoi(argv[1]);
    int to = atoi(argv[2]);
    int nw = atoi(argv[3]);
    
    //seq versin
    int ns = 0;
    {
        utimer tseq("Seq");
        for(int i = from; i<= to; i++)
            if(is_prime(i))
                ns++;
    }
    cout << ns << endl;
    
    //OMP version
    int no = 0;
    {
        utimer tomp("OMP");
        /*parallel => we want to use nw threads
          for => we want to parallelize the iteration of this for
          firstprivate => deals with the visibility of from and to in the parallel section
                          with fp every thread has his own copy of the value
          reduction => for accumulate the resalt in no */
        #pragma omp parallel for num_threads(nw) firstprivate(from, to) reduction(+:no)
            for(int i = from; i<= to; i++)
            if(is_prime(i))
                no++;
    }
    cout << no << endl;
    /*
    // ff version
    int nff = 0;
    {
        utimer tff("FF");
        //true => for use spean lock instead of the normal locks => tell the sistem to be fastest possible
        //the setup of the parallelForReduce spent a lot of time so 
        //is better declare it befor respect when you have to use and reuse it
        ParallelForReduce<int> pfr( nw, true);
        //we don't use the chunkSize(=0 Staticaly divided) 
        pfr.parallel_reduce(nff, 0,
                            from, to, 1, 0,
                            [&](const int i, int& nff) {if(is_prime(i)) nff++;},
                            [] (int& s, const int e) {s+=e;},
                            nw);
    }
    cout << nff << endl;

    //static ff version
    //we use chunk size
    int snff = 0;
    {   //chunkSize = 1000 => the iter are divided in 1000 by 1000
        utimer stff("StFF");
        ff::parallel_reduce(nff, 0,
                            from, to, 1, 1000,
                            [&](const int i, int& nff) {if(is_prime(i)) snff++;},
                            [] (int& s, const int e) {s+=e;},
                            nw);
    }
    cout << snff << endl;
    
    //ff Farm verison
    int nfarm = 0;
    {
        utimer tfarm("FFarm");
    
        //emitter (source)
        class src : public ff_node_t<int> {
            private:
                int f, t;
            public: 
                src(int f, int t): f(f),t(t){}

                
                int * svc(int *){
                    for(int i=f; i<=t; i++)
                        ff_send_out(new int(i));
                    return EOS;
                }
        };

        //worker
        class filter : public ff_node_t<int> {
            int * svc(int * t){
                if(is_prime(*t))
                    *t = 1;
                else
                    *t = 0;
                return t;
            }
        };

        //collector
        class drain: public ff_node_t<int>{
            private:
                int sum;
                int * res;
            public:
                drain(int* res):res(res), sum(0) {}
                
                int * svc(int * t){
                    sum+= *t;
                    return GO_ON;
                }

                void svc_end(){
                    *res = sum;
                }
        };
    
    // to compute single item we have
    // one thread that send to the worker
    // the worke that compute and send to the collector 
    // and the collector that make a sum so
    // we have 2 comunication inspite of 1/n° of chunks that we use in the parFor
    // this is the reasone why this solution is very bed for this problem also with scheduling ondemand

        vector <ff_node*> w;
        for(int i=0; i<nw; i++)
            w.push_back(new filter());
        ff_farm fa;
        src e(from, to);
        fa.add_emitter(e);
        fa.add_workers(w);
        drain d(&nfarm);
        fa.add_collector(&d);
        // fa.set_scheduling_ondemand();

        fa.run_and_wait_end();
    }

    cout << nfarm << endl;
*/
    return 0;
}