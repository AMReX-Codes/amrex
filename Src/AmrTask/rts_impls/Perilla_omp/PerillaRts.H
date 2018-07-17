#include <AMReX_Amr.H>
#include <PerillaConfig.H>


namespace perilla{

class PerillaRTS
{
    private:
	int _nWrks;
	void RTS_Init();
	int _rank, _nProcs;

    public:
	PerillaRTS(){
	    _nWrks=1;
	    char* nWrks= getenv("NWORKERS");
	    if(nWrks) _nWrks= atoi(nWrks);
	}
	PerillaRTS(int nWrks):_nWrks(nWrks){}
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
