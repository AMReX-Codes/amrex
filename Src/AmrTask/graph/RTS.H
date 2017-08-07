#include "AMReX_AbstractTask.H"
#include "AMReX_TaskGraph.H"
#include <iostream>
#include <queue>
using namespace std;

namespace amrex{

    class RTS{
	private:
	    typedef std::queue<Task*> TaskQueue;
	    typedef std::queue<Data*> MessageQueue;
	    //Task Queues
	    TaskQueue WaitingQueue; 
	    TaskQueue DataFetchingQueue;  //used in Pull model
	    TaskQueue ReadyQueue; 
	    TaskQueue RunningQueue; 
	    MessageQueue MsgQueue; 

	public:
	    int RTS_ProcCount();
	    int RTS_MyProc();
	    int RTS_WorkerThreadCount();
	    int RTS_MyWorkerThread();
	    void RTS_Init();
	    void RTS_Init(int *rank, int *nProcs);
	    void RTS_Run(void *graph, bool keepInitialGraph);
	    void RTS_Finalize();
    };
}
