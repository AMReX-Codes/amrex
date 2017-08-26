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

	public:
            RTS():_nWrks(1){}
            RTS(int nWrks):_nWrks(nWrks){}
	    int RTS_ProcCount();
	    int RTS_MyProc();
	    int RTS_WorkerThreadCount();
	    int RTS_MyWorkerThread();
	    void RTS_Init();
	    void RTS_Init(int *rank, int *nProcs);
	    void RTS_Run(void *graph);
	    void RTS_Finalize();
            double Time();
            void Barrier();
            template<typename T> void ReductionSum(T *local, T *global, int length, int root){
                ReductionSum_impl(local, global, length, root);
            }
    };
}
