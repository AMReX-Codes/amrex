//Question? email tannguyen@lbl.gov
//Created 07-19-2017
//ompodification 08-14-2017
#include <mpi.h>
#include <sched.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <omp.h>
#include "PerillaRts.H"

using namespace perilla;
#ifdef PERILLA_DEBUG
#include <PerillaMemCheck.H>
PerillaMemCheck memcheck;
#endif

#include <iostream>
#include <queue>
using namespace std;
#include <cassert>

namespace perilla{
    Amr* amrptr;

    int RTS::ProcCount(){
	return _nProcs;
    }

    int RTS::MyProc(){
	return _rank;
    }

    int RTS::WorkerThreadCount(){
	return _nWrks;
    }

    int RTS::MyWorkerThread(){
	return 0;
    }

    void RTS::runAMR(Amr* amr, int max_step, Real stop_time){
        while ( amr->okToContinue() &&
              (amr->levelSteps(0) < max_step || max_step < 0) &&
              (amr->cumTime() < stop_time || stop_time < 0.0) )
            
        {
            // Do a coarse timestep, which calls one or multiple timestep updates (i.e. timeStep()) at each AMR level
            amr->coarseTimeStep(stop_time);
        }
    }

    void InitializeMPI(){
	int provided;
	MPI_Init_thread(0, 0, MPI_THREAD_FUNNELED, &provided);
	if(provided == MPI_THREAD_SINGLE){//with this MPI, process can't spawn threads
	    cerr << "Spawning threads is not allowed by the MPI implementation" << std::endl;;
	}
    }

    void RTS::RTS_Init(){
	amrptr= NULL;
    }

    void RTS::Init(){
        InitializeMPI();
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &_nProcs);
        RTS_Init();
    }

    void RTS::Init(int rank, int nProcs){
        _rank= rank;
	_nProcs= nProcs;
	RTS_Init();
    }

    void RTS::Finalize(){
#ifdef PERILLA_DEBUG
        memcheck.report();
#endif
    }

    void RTS::Iterate(void* amrGraph, int max_step, Real stop_time){
	    Perilla::max_step=max_step;
	    assert(amrGraph);
	    amrptr= (Amr*)amrGraph;
            runAMR(amrptr, max_step, stop_time);
    }

#if 0
    const double kMicro = 1.0e-6;
    double RTS::Time()
    {
	struct timeval TV;

	const int RC = gettimeofday(&TV, NULL);
	if(RC == -1)
	{
	    printf("ERROR: Bad call to gettimeofday\n");
	    return(-1);
	}
	return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );
    } 
#endif

    void RTS::Barrier(){
	//nothing
    }

}//end namespace

