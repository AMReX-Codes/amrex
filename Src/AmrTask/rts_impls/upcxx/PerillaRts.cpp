//Question? email tannguyen@lbl.gov
//Created 07-19-2017
//Last modification 08-14-2017
#include <mpi.h>
#include <sched.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <mylock.h>
#include <pthread.h>
#include "PerillaRts.H"

#include <iostream>
#include <queue>
using namespace std;
#include <cassert>

using namespace perilla;
#ifdef PERILLA_DEBUG
#include <PerillaMemCheck.H>
PerillaMemCheck memcheck;
#endif

namespace perilla{
    Amr* amrptr;
    struct RtsDomain{
	pthread_t *_threads;
	int _size;
	MyLock _lock;
	RtsDomain():_threads(NULL), _size(0){};
	~RtsDomain(){
	    free(_threads);
	}
    };
    int numa_nodes;
    RtsDomain *dom;
    MyLock _l;
    volatile char startSignal=0;
    pthread_mutex_t startLock= PTHREAD_MUTEX_INITIALIZER;

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

    struct argT {
	int numaID;
	int tid;
	int g_tid;
	int nThreads;
	int nTotalThreads;
	int max_step;
	Real stop_time;
	RTS* thisRTS;
    };

    void RTS::runAMR(Amr* amr, int max_step, Real stop_time){
        while (amr->okToContinue() &&
              (amr->levelSteps(0) < max_step || max_step < 0) &&
              (amr->cumTime() < stop_time || stop_time < 0.0) )
            
        {
            // Do a coarse timestep, which calls one or multiple timestep updates (i.e. timeStep()) at each AMR level
            amr->coarseTimeStep(stop_time);
        }
    }

#ifdef USE_PERILLA_PTHREADS
    void run(void* threadInfo){
	argT *args= (argT*)threadInfo;
	int numaID= args->numaID;
	int tid= args->tid;
	int g_tid= args->g_tid;
	int nThreads= args->nThreads;
	int nTotalThreads= args->nTotalThreads;
	int max_step= args->max_step;
	Real stop_time= args->stop_time;
	RTS* rts= args->thisRTS;
	Perilla::registerId(g_tid);
	//done with thread id setup, now wait for the start signal from master
        pthread_mutex_lock(&startLock);
	startSignal++;
        pthread_mutex_unlock(&startLock);
	while(startSignal!= nTotalThreads){}
        rts->runAMR(amrptr, max_step, stop_time);
    }
#endif

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
	    assert(amrGraph);
	    Perilla::max_step= max_step;
	    amrptr= (Amr*)amrGraph;
            WorkerThread::init();
#ifndef USE_PERILLA_PTHREADS
            runAMR(amrptr, max_step, stop_time);
#else
	    int numa_nodes= perilla::NUM_THREAD_TEAMS;
	    int worker_per_numa = perilla::NUM_THREADS_PER_TEAM;
            int _nWrks= numa_nodes*worker_per_numa;
	    int base=0; 
	    int localID=-1;
	    //create a list of persistent threads for each NUMA node
	    cpu_set_t cpuset;
	    pthread_attr_t attr;
	    pthread_attr_init(&attr);
	    dom= new RtsDomain[numa_nodes];
	    for(int i=0; i<numa_nodes; i++){
		dom[i]._threads= new pthread_t[worker_per_numa];
	    }
	    for(int i=0, domNo=-1; i<_nWrks; i++){
		localID++;
		if(localID==0){
		    domNo++;
		}
		CPU_ZERO(&cpuset);
		CPU_SET(base+localID, &cpuset);
		if(! (localID==0 && domNo==0)){
		    pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpuset);
		    argT* arg= new argT;
		    arg->numaID= domNo;
		    arg->tid= localID;
		    arg->g_tid= domNo*worker_per_numa+localID;
		    arg->nThreads= worker_per_numa;
		    arg->nTotalThreads= _nWrks;
		    arg->thisRTS= this;
		    arg->max_step= max_step;
		    arg->stop_time= stop_time;
		    int err = pthread_create(&(dom[domNo]._threads[localID]), &attr, (void*(*)(void*))run, arg);
		}else{ //master thread
		    dom[domNo]._threads[localID]= pthread_self();
		    Perilla::registerId(0);
		    //enable worker threads to start computing
        	    pthread_mutex_lock(&startLock);
	   	    startSignal++;
	            pthread_mutex_unlock(&startLock);
                }
		dom[domNo]._size++;
		if(localID == (worker_per_numa-1)){
		    localID=-1;
		    base+= worker_per_numa;
		}
	    }
	    while(startSignal!= _nWrks){}//wait until all threads have done the setup phase
            runAMR(amrptr, max_step, stop_time);
	    for(int i=1; i<_nWrks; i++) pthread_join(dom[i/worker_per_numa]._threads[i%worker_per_numa], NULL);
#endif
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

