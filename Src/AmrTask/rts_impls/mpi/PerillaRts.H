#ifndef _PERILLA_RTS_
#define _PERILLA_RTS_

#include <AMReX_Amr.H>
#include <PerillaConfig.H>
#include <RGIter.H>

namespace perilla{

struct _workerThreadInfo{
    int _tid; //thread id in local group
    int _size; //number of threads in the group
};

struct _threadInfo{
    bool _isComm; //whether this thread handles communication
    int _wtid; //worker thread id (-1 if this thread is decicated to communication)
    int _nWts; //number of thread groups
};

class RTS
{
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
	void Iterate(void *graph, int max_step, Real stop_time);
	void Finalize();
//	double Time();
	void Barrier();
	void runAMR(Amr* amrptr, int max_step, Real stop_time);
        void invokeOnDemand(std::vector<RegionGraph *> rg, RGIter *rgi);
};

}
#endif
