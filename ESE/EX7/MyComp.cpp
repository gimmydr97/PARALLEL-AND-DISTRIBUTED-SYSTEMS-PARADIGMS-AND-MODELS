//Assignment 7: write a parallel program getting a vector of n floating point numbers 
//with first position 0.0, last position 100.0 and intermediate positions 25.0 
//and computing (in parallel) k iterations each re-computing position i in the vector (i in [1,n-2]) 
//as the average of previous values of positions i, i-1 and i+1 (positions 0 and n never change). 
//Two versions are required, one using FastFlow and the other one using any other framework (C++ threads, OpenMP, GrPPi).

#include <iostream>
#include <vector>
#include "utimer.cpp"
//#include <ff/ff.hpp>
//#include <ff/parallel_for.hpp>

using namespace std;


//g++ -fopenmp  -pthread -O3 MyComp.cpp -o ex7
//./ex7 100000000 10 2(= parDeg)
int main(int argc, char * argv[]){

    int n = atoi(argv[1]);
    int k = atoi(argv[2]);
    int w = atoi(argv[3]);

    //25.0 is the costant used to initialize every position
    vector<double> vseq1(n,25.0);
    vseq1[0] = 0.0;
    vseq1[n-1] = 100.0;

    //vector uset to write the new value
    vector<double> vseq2(n);
    vseq2[0] = 0.0;
    vseq2[n-1] = 100.0;

    auto pv = [&] (vector<double> v){
        for(int i=0; i<n; i++)
            cout << v[i] << " ";
        cout << endl;
        return;
    };

    //pv(v1);

    long seq;
    {
        utimer tseq("Seq", &seq);
        for(int kk=0; kk<k; kk++){
            for(int i=1; i<n-1; i++)
                vseq2[i] = (vseq1[i-1]+vseq1[i]+vseq1[i+1])/3.0;
            kk++;
            //pv(v1);
            for(int i=1; i<n-1; i++)
                vseq1[i] = (vseq2[i-1]+vseq2[i]+vseq2[i+1])/3.0;
            //pv(v2);
        }
    }

    cout << "usec / iter = " << seq/k << endl;
    
    //25.0 is the costant used to initialize every position
    vector<double> vomp1(n,25.0);
    vomp1[0] = 0.0;
    vomp1[n-1] = 100.0;

    //vector uset to write the new value
    vector<double> vomp2(n);
    vomp2[0] = 0.0;
    vomp2[n-1] = 100.0;

    //pv(v2);

     long omp;
    {
        utimer tomp("OMP", &omp);
        for(int kk=0; kk<k; kk++){
            #pragma omp parallel for num_threads(w)
            for(int i=1; i<n-1; i++)
                vomp2[i] = (vomp1[i-1]+vomp1[i]+vomp1[i+1])/3.0;
            kk++;
            //pv(v1);
            #pragma omp parallel for num_threads(w)
            for(int i=1; i<n-1; i++)
                vomp1[i] = (vomp2[i-1]+vomp2[i]+vomp2[i+1])/3.0;
            //pv(v2);
        }
    }

    cout << "usec / iter = " << omp/k << endl;
    
    /*
    //25.0 is the costant used to initialize every position
    vector<double> vff1(n,25.0);
    vff1[0] = 0.0;
    vff1[n-1] = 100.0;

    //vector uset to write the new value
    vector<double> vff2(n);
    vff2[0] = 0.0;
    vff2[n-1] = 100.0; 

    long ff;
    {
        utimer tff("FF", &ff);
            //ff::ParallelFor pf(nw,true);
        for(int kk=0; kk<k; kk++){
            //to be computed in parallel
            //pf.parallel_for(...)
            ff::parallel_for(1,n-1,
                             1, //increment
                             0, //chunkSize
                            [&](const int i) {
                                //for(int i =1; i<n-1; i++)
                                  vff2[i] = (vff1[i-1]+vff1[i]+vff1[i+1])/3.0;
                            },nw);
            k++;
            //pf.parallel_for(...)
            ff::parallel_for(1,n-1,
                             1, //increment
                             0, //chunkSize
                            [&](const int i) {
                                //for(int i =1; i<n-1; i++)
                                  vff1[i] = (vff2[i-1]+vff2[i]+vff2[i+1])/3.0;
                            }, nw);
        }
    }

    cout << "usec / iter = " << ff/k << endl;
    */
    return(0);


}