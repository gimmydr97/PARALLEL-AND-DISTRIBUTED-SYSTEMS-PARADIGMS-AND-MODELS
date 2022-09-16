/*
Jacobi (AJ):

The Jacobi iterative method computes the result of a system of equations Ax = B 
    with x vector of variable of length n, 
    A matrix of coefficients of dimension n by n
    B vector of known terms of length n
iteratively computing a new approximation of the values of the different variables 
according to the formula:
    
starting from some initial assignment of each x_i (e.g. 0).
We require to implement the Jacobi method with both native C++ threads and FastFlow.
*/
#include <iostream>
#include<functional>
#include <vector>
#include <math.h>
#include "utimer.cpp" 
#include<optional>
#include<future>
#include<queue>
#include <mutex>
#include<barrier>
#include<thread>
#include <ff/ff.hpp>
#include <ff/parallel_for.hpp>
#include <grppi/common/patterns.h>
#include <grppi/mapreduce.h>
#include <grppi/dyn/dynamic_execution.h>


using namespace std;
using namespace ff;
using namespace grppi;

#define ERRTOL pow(10,-6)

//norm for stop criteria
float xnorm (vector<float> xatt, vector<float> xprev, int n){
        float xn = 0;
        for(int j = 0; j < n; j++)
            xn = xn + pow(xatt[j]-xprev[j],2);
        return  sqrt(xn);
}

//create A matrix with random values but like a Dominant diagonal matrix 
vector<vector<float>> randAmatrix(int seed, int max, int n){

    srand(seed);
    vector<vector<float>> a(n, vector<float>(n));
    //vector for conteining the sum of the element in a row without the diagonal element
    vector<float> sumOnRow(n);

    for(int i = 0; i < n; ++i){
        for(int j = 0;  j < n; ++j){
            a[i][j] = rand() % max;
            if( i != j)
                sumOnRow[i] += a[i][j];
        }
        //if the diagonal elem in a row is smmaler than the sumOnRow 
        //adjust the diagonal value to be greater than the sumOnRow
        if(a[i][i] < sumOnRow[i])
            a[i][i] = sumOnRow[i] + (rand() % max);
    }
    /*
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j)
            cout<<a[i][j]<<'\t';
        cout<<'\n';
    }*/

    return a;
}

//create b random vector 
vector<float> randBvect(int seed, int max, int n){

    srand(seed);
    vector<float> b(n);

    for(int i = 0; i < n; ++i)
        b[i] = rand() % max;
    /*
    for(int i = 0; i < n; ++i)
            cout<<b[i]<<'\t';
    */
    return b;
}



//g++-11 -O3 -pthread -std=c++20 jacobi.cpp -o jacobi 
int main(int argc, char * argv[]){

    int n = atoi(argv[1]); //size of the vectors and matrix
    int nw = atoi(argv[2]); //par degree
    int k = atoi(argv[3]); //number of iteration
    int seed = atoi(argv[4]); //seed for the creation of matrix
    int max = atoi(argv[5]); //max value for the random matrix

    int pardeg = nw;
    int pardegRed = nw;

    //print vector 
    auto pv = [&] (vector<float> v){
        for(int i=0; i<n; i++)
            cout << v[i] << " ";
        cout << endl;
        return;
    };

    vector<vector<float>> A(n, vector<float>(n)); //matrix of coefficient of dimension n by n
    A = randAmatrix(seed, max, n);

    vector<float> b(n);   // vector of known terms of lenght n
    b = randBvect(seed, max, n);

    //function that compute the position i of the x vector in a certain iteration
    auto computeXi = [&] (const int i, vector<float> xprev, float* xi, int policy){

        float sum = 0;
        
        
        switch ( policy )
        {
         case 0: //seq solution

            for(int j = 0; j < n; j++)
                sum += (A[i][j] * xprev[j]);
            break;

         case 1: //parallel reduce with grppi solution
         {  //execution policy
            dynamic_execution par = grppi::parallel_execution_native{pardegRed};
            //the input data sets xprev and A[i] are specified by ranges packaged into a grppi::zip_view
            /*alternative deprecated solution is grppi::map_reduce(par,
                                                                   begin(xprev),end(xprev),
                                                                   0.0F,
                                                                   [](float a, float x) { return a*x; },
                                                                   [](float sum, float elem) { return sum+elem; },
                                                                   begin(A[i])
                                                                   );
            */
            sum = grppi::map_reduce(par,
                                    grppi::zip_view(xprev,A[i]),
                                    0.0F,
                                    [](float a, float x) { return a*x; },
                                    [](float sum, float elem) { return sum+elem; }
                                    );
            break;
         }
         case 2: //paralle reduce with fastflow solution
         {  
        
            ParallelForReduce<float> pfr(pardegRed,true,true);
            pfr.parallel_reduce(sum,0,0,n,1,0,
                                [&](const float j, float& mysum){ mysum += A[i][j] * xprev[j];},
                                [&](float& s, const float ps){ s += ps;},
                                pardegRed);
            break;
         }
         
        }
         
        sum = sum - (A[i][i] * xprev[i]);
        *xi = (float)(1 / A[i][i]) * (b[i] - sum);
        
        return i;
    };
    
    //sequential solution
    vector<float> x1(n,0);
    vector<float> x2(n);
    int kk = 1;

    long seq;
    {
        utimer tseq("Seq", &seq);
        do{ 
            for(int i=0; i < n; i++)
                if( kk%2 == 1 ){
                    computeXi(i,x1, &x2[i],0);
                    
                }else{
                    computeXi(i,x2, &x1[i],0);
                }
            kk++;
        } while ((xnorm(x1, x2, n) > ERRTOL) && ( kk <= k));  
    }

    /*
    cout <<" iter-seq " << kk-2 << " = " ;
    pv(x2);
    cout <<" iter-seq " << kk-1 << " = " ;
    pv(x1);
    */
    cout << "sequsec / iter = " << (float)seq/(kk -1) << endl;
    cout << "number of iteration : " << (kk -1) << endl;
    
    //parallel solution with static scheduling
    vector<float> x1parS(n,0);
    vector<float> x2parS(n);
    vector<thread> tidsS(nw);

    barrier bar(pardeg, [&](){return;}); 

    //task for the threads of the static scheduling solution
    auto compute = [&] (int i){
        
        for(int kpar = 1; kpar < kk; kpar++){
            
            for(int j=i; j < n ; j+=pardeg){
                if( kpar%2 == 1 ){
                    computeXi(j,x1parS, &x2parS[j],0);
                }else{
                    computeXi(j,x2parS, &x1parS[j],0);
                }
            }
            bar.arrive_and_wait();
            
        }
        return;
    };
    
    long parStatic;
    {
        utimer tparStat("ParStatic", &parStatic);
        
        //spone the threads
        for(int i = 0; i<pardeg; i++){
            tidsS[i]= thread(compute,i);
        }

        //join for threads
        for(int i = 0; i<pardeg; i++)
            tidsS[i].join();

    }
    /*
     cout <<" iter-parS " << kk-2 << " = " ;
     pv(x2parS);
     cout <<" iter-parS " << kk-1 << " = " ;
     pv(x1parS);
    */
    cout << "parSusec / iter = " << (float)parStatic/(kk -1) << endl;

    //solution with FastFlow with obecj ParallelFor 
    vector<float> x1ffo(n,0);
    vector<float> x2ffo(n);

    long parffo;
    {
        utimer tparffO("ParFFObject", &parffo);
        ParallelFor pf(pardeg,false);

        for(int kpar=1; kpar < kk; kpar++){

            if( kpar%2 == 1 ){
                pf.parallel_for(0,n,
                             1, //increment
                             0, //chunkSize
                             [&] (const int i){
                                computeXi(i,x1ffo, &x2ffo[i],0);
                             }, pardeg);
            }else{
                pf.parallel_for(0,n,
                             1, //increment
                             0, //chunkSize
                             [&] (const int i){ 
                                computeXi(i,x2ffo, &x1ffo[i],0);
                              }, pardeg);
            }

        }
    }
    /*
    cout <<" iter-ffO" << kk-2 << " = " ;
    pv(x2ffo);
    cout <<" iter-ffO" << kk-1 << " = " ;
    pv(x1ffo);
    */
    cout << "parFFObject / iter = " << (float)parffo/(kk -1) << endl;

    /*
    mutex gl; //global lock
    //ptask return a future. We us it because we want to be sure to close the threadpoll only when all the result has been computed
    deque<packaged_task<void()>> tasks; // queue of task
    condition_variable cond; //condition used for notify and wait of threads
    bool stop = false; //flag used to avoid infinite wait
    
    //function used by the emitter for push tasks in the dequeue
    auto submit = [&] (packaged_task<void()> f){
        {
            unique_lock<mutex> lock(gl);
            tasks.push_back(move(f));
        }
        //if anyone is waiting for a task one of this is awaken
        cond.notify_one();
    };
    
    //function used by the consumer for pop ane task in mutual exclusion
    auto consume = [&] (int i){
        while(true){
            packaged_task<void()> t;
            {
                unique_lock<mutex> lock(gl);
                //cond on witch the ts are wating. If the queue is an empty queue
                //the wait function relise the lock and the thread is descheduled waiting for a notify_one
                // if the stop == true the thread are awaken from his wait and go to finish his life
                cond.wait( lock,
                           [&](){return(!tasks.empty()||stop);}); // || stop agginto per evitare l'attesa infinita
                //if the queue is not empty the thread take one task and delete the task from the dequeue
                if(!tasks.empty()){
                    t = move(tasks.front());
                    tasks.pop_front();
                }
                //condizione aggiunta per evitarel'attesa infinita; se i tasks da svolgere sono finiti termino 
                if(stop)
                    return;
            }  
            //t is function that compute the new element of the x vector 
            t();
        }
    };

    //function that in mutual exclusion put stop to true for stop the ThreadPool
    auto stopTp= [&](){
            {
                unique_lock<mutex> lock(gl);
                stop = true;
            }
            //tells oll the threads blocked on the wait cond that there was a change
            cond.notify_all();
    };
    
    //parallel solution with dynamic scheduling 
    vector<float> x1par(n,0);
    vector<float> x2par(n);
    vector<thread> tids(nw);
   
    long parDinamic;
    {
        utimer tparDin("ParDinamic", &parDinamic);

        //initializing the threads and all waiting on the cond because in the deq there is't task
        for(int i = 0; i<pardeg; i++){
            tids[i]= thread(consume,i);
        }

        //vector for store the futures
        vector<future<void>> vf(n);
        
        for(int kpar = 1; kpar < kk; kpar++){
            
            // fx is an object of function<float()> type so a function without arguments
            function<float()> fx; 
            
            //adding the task into the dequeue as an emitter
            for(int i = 0; i<n; i++){
                
                if( kpar%2 == 1 )
                    //bind the function in fx 
                    fx = (bind(computeXi,i,x1par,&x2par[i],0));
                else
                    fx = (bind(computeXi,i,x2par,&x1par[i],0));

                //create the packaged task with fx as task
                packaged_task<void()> pt(fx);
                //store the promise of the resutl of pt in the futures vector
                vf[i] = pt.get_future();
                //submit the packaged task in the dequeue
                submit(move(pt));
            
            }
            //we are waiting to complete the vector x of this iteration to move on to the next one
            //is the syncronization phase
            for(int i = 0; i<n; i++)
                //the get method wait for the result of this vf[i] future
                vf[i].get();

        }
        
        //stop the ThreadPool
        stopTp();

        //join the threads
        for(int i = 0; i<pardeg; i++)
            tids[i].join();
        
    }
    
    cout <<" iter-par " << kk-2 << " = " ;
    pv(x2par);
    cout <<" iter-par " << kk-1 << " = " ;
    pv(x1par);
    
    cout << "parusec / iter = " << (float)parDinamic/(kk -1) << endl;

    //FastFlow solution with static call of the ParallelFor 
    vector<float> x1ffs(n,0);
    vector<float> x2ffs(n);

    long parffs;
    {
        utimer tparffS("ParFFStatic", &parffs);

        for(int kpar=1; kpar < kk; kpar++){

            if( kpar%2 == 1 ){
                parallel_for(0,n,
                             1, //increment
                             10, //chunkSize
                             [&] (const int i){
                                computeXi(i,x1ffs, &x2ffs[i],0);
                             }, pardeg);
            }else{
                parallel_for(0,n,
                             1, //increment
                             10, //chunkSize
                             [&] (const int i){ 
                                computeXi(i,x2ffs, &x1ffs[i],0);
                              }, pardeg);
            }

        }
    }
    
    cout <<" iter-ffS" << kk-2 << " = " ;
    pv(x2ffs);
    cout <<" iter-ffS" << kk-1 << " = " ;
    pv(x1ffs);
    
    cout << "parFFStatic / iter = " << (float)parffs/(kk -1) << endl;
    */
                                        //comparation with staitc scheduling solution results
    cout<< "the results are equal? " << ((std::equal(x1.begin(),x1.end(),x1parS.begin()) && 
                                           std::equal(x2.begin(),x2.end(),x2parS.begin()))
                                        &&
                                        //comparation with FastFlow with Object ParalleFor solution results
                                          (std::equal(x1.begin(),x1.end(),x1ffo.begin()) && 
                                           std::equal(x2.begin(),x2.end(),x2ffo.begin()))
                                        /*
                                        &&  
                                          //comparation with dynamic scheduling solution results
                                          (std::equal(x1.begin(),x1.end(),x1par.begin()) && 
                                           std::equal(x2.begin(),x2.end(),x2par.begin()))
                                        &&
                                        //comparation with FastFlow with static call of ParalleFor solution results
                                          (std::equal(x1.begin(),x1.end(),x1ffs.begin()) && 
                                           std::equal(x2.begin(),x2.end(),x2ffs.begin())) 
                                        */
                                        ? "yes" : "no") << endl;
                                        
                                   
    return 0;
}
