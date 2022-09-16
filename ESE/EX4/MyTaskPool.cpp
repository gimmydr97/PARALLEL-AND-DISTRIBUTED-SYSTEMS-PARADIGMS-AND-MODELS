/*Assignement 4: implement a task pool using C++ threads and STDLIB. 
  The task pool should support parallelism degree as creation parameter and
  then the submission of tasks that are computed and produce results through side effects. 
  Alternatively they have to produce a result that can be retrieved by the task submitter.*/

#include<iostream>
#include<functional>
#include<vector>
#include<chrono>
#include<optional>
#include<ctime>
#include<future>
#include<queue>
#include <mutex>
#include <cmath>


using namespace std;
using namespace std::literals::chrono_literals;

//g++-11 -std=c++20 -pthread -O3  MyTaskPool.cpp -o es4
//./es4 4 10 10 1000 
int main(int argc , char* argv[]){

    if(argc == 1) cout << "Usage: " << argv[0] << "nw m ta tf" << endl; 

    int nw = atoi(argv[1]); //par degree
    int  n = atoi(argv[2]); //number of task to be computed
    int ta = atoi(argv[3]); //interrarrival time
    int tf = atoi(argv[4]); //time used to compute the single items

    mutex ll; //global lock

    //ptask return a future. We us it because we want to be sure to close the threadpoll only when all the result has been computed
    deque<packaged_task<float()>> tasks; // queue of task
    condition_variable cond; //cond used to notify and wait the different ts
    bool stop = false; //aggiunta per evitare l'attesa infinita

    //function used for put task in dequeue
    auto submit = [&] (packaged_task<float()>& f){
        {
            unique_lock<mutex> lock(ll);
            tasks.push_back(move(f));
        }
        //if anyone is waiting for a task one of this is awaken
        cond.notify_one();
    };
    
    //function used by the consumer
    auto body = [&] (int i){
        while(true){
            packaged_task<float()> t;
            {
                unique_lock<mutex> lock(ll);
                //cond on witch the ts are wating. If the queue is an empty queue
                //the wait function relise the lock and the thread is descheduled waiting for a notify_one
                // if the stop == true the thread are awaken from his wait and go to finish his life
                cond.wait( lock,
                           [&](){return(!tasks.empty()||stop);}); // || stop agginto per evitare l'attesa infinita
                //if the queue is not empty the thread take one task and delete she from the dequeue
                if(!tasks.empty()){
                    t = move(tasks.front());
                    tasks.pop_front();
                }
                //condizione aggiunta per evitarel'attesa infinita; se i tasks da svolgere sono finiti termino 
                if(stop)
                    return;
            }  
            //t is function that return a float without parameters
            t();
        }
    };


    auto f = [tf](float x){
        this_thread::sleep_for(tf * 10ms);//wait did by sleap
        auto res = sqrt(x);
        //side effect did by printing on the standard output
        cout << res << endl;
        return(res);
    };

    vector<thread> tids(nw);

    //initializing the threads and all waiting on the cond because in the deq there is't task
    for(int i = 0; i<nw; i++){
        tids[i]= thread(body,i);
    }
    
    //vector for store the futures
    vector<future<float>> vf(n);

    //adding the task into the dequeue
    for(int i = 0; i<n; i++){
        this_thread::sleep_for(ta* 1ms); //simulation of interarrival time
        auto x = (double) i;
        cout <<"Task " << x <<endl;
        // fx è un oggetto di tipo function<float()> quindi una funzine senza argomenti
        auto fx = (bind(f,x));
        //encapsulate the binded function in a ptask
        packaged_task<float()> pt(fx);
        //and we store her future in vf
        vf[i] = pt.get_future();
        submit(pt);
    }

    //function that in mutual exclusion put stop to true
    auto stopTp= [&](){
        {
            unique_lock<mutex> lock(ll);
            stop = true;
        }
        //tells oll the threads blocked on the wait cond that somethis appened
        cond.notify_all();
    };

    //we are waiting for the results
    for(int i = 0; i<n; i++)
        cout << "Task computed " << vf[i].get() <<endl;

    //when we have finish to submit the task we declare stop to true
    stopTp();
    cout << "Awaiting threads ..." << endl;
    for(int i = 0; i<nw; i++)
        tids[i].join();
    
    return 0;

} 

//se non aggiungiamo un qualcosa per far uscire i thread dal while perenne
//quando le n task finiscono i threads rimarranno in attesa all'infinito sulla cond 
//non potendo mai fare la join()

