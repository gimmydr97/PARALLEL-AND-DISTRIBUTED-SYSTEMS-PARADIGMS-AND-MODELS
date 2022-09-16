/*Assignment 5: implementation of a farm using OpenMP. 
  Tasks to be computed have to be provided through a for loop with iterations 
  producing one of the input tasks and then awaiting for Ta time (parameter) before executing next iteration. 
  The parallelism degree of the farm should be a parameter. 
  Each task should spent some parametric amount of time (possibly using the active_wait functioN) 
  to produce the final result.
*/

#include<iostream>
#include<chrono>
#include<thread>
#include<omp.h>

using namespace std;
using namespace std::literals::chrono_literals;

void offload(int task, int tw){
    #pragma omp task
    {
        this_thread::sleep_for(tw * 1ms);
        #pragma omp critical //non cambia il risultato
            cout<< "task " << task << " = " << task * task << endl;
    }
    return;
}

//g++-11  -fopenmp  OpenMPFarm.cpp -o es5
//time ./es5 10 10 1000 4(par degree)

int main(int argc , char* argv[]){

    if(argc == 1) cout << "Usage: " << argv[0] << "nw n ta tw" << endl; 

    int  n = atoi(argv[1]); //number of task to be computed
    int ta = atoi(argv[2]); //interarrival time
    int tw = atoi(argv[3]); //time used to compute the single items
    int nw = atoi(argv[4]); //par degree

    #pragma omp parallel num_threads(nw)
    {
        #pragma omp master //can be also omp single
        {
            for(int i=0; i < n; i++){
                int task = i;
                this_thread::sleep_for(ta * 1ms);
                offload(task,tw);
            }
        }
        #pragma omp taskwait
    }
    return 0;
}
