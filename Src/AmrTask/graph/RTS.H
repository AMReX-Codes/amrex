#ifndef _RTS
#define _RTS
#include "AMReX_AbstractTask.H"
#include "AMReX_TaskGraph.H"
#include <iostream>
#include <queue>
#include "rts_graphimpl.H"
#include <sys/time.h>

using namespace std;

namespace amrex{
    struct _workerThreadInfo{
	int _tid; //thread id in local group
        int _size; //number of threads in the group
    };

    struct _threadInfo{
        bool _isComm; //whether this thread handles communication
        int _wtid; //worker thread id (-1 if this thread is decicated to communication)
        int _nWts; //number of thread groups
    };

    class RTS{
	private:
            int _nWrks;
            void RTS_Init();
	    int _rank, _nProcs;

	public:
            RTS(){
		_nWrks=1; 
		char* nWrks= getenv("NWORKERS");
		if(nWrks) _nWrks= atoi(nWrks);
	    }
            RTS(int nWrks):_nWrks(nWrks){}
	    int ProcCount();
	    int MyProc();
	    int WorkerThreadCount();
	    int MyWorkerThread();
	    void Init(); //Build the runtime system from scratch
	    void Init(int rank, int nProcs);//Build the runtime system on pre-existing MPI processes
	    void Iterate(void *graph);
	    void Finalize();
            double Time();
            void Barrier();
            template<typename T> void ReductionSum(T *local, T *global, int length, int root){
                ReductionSum_impl(local, global, length, root);
            }
    };
}
#endif
