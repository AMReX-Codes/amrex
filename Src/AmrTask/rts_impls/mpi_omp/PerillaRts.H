#include <AMReX_Amr.H>
#include <PerillaConfig.H>


namespace perilla{

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
};

}
