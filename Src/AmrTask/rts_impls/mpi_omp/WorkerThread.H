#ifndef P_WORKERTHREAD_H
#define P_WORKERTHREAD_H

#include <PerillaConfig.H>
#include <Barrier.H>
#include <omp.h>

namespace perilla{

    class WorkerThread
    {
	static void* team_shared_memory[perilla::NUM_THREAD_TEAMS];  
        int tid;
	public:
        static void init();
        static Barrier *globalBarrier;
        static Barrier localBarriers[perilla::NUM_THREAD_TEAMS];
	static int perilla_tid();
	static int perilla_wtid();
	static int perilla_wid();
	static int perilla_nWorkerThreads();
	static int perilla_nWorkers();
	static int perilla_nTeamThreads();
	static bool perilla_isMasterWorkerThread();  
	static bool perilla_isMasterThread();
	static bool perilla_isCommunicationThread();
	static bool isMyRegion(int workerID, int regionID);
	static void setTeamSharedMemory(void* dummy, int tid, int tg);  
	static void* getTeamSharedMemory(int tg);
        static void syncWorkers();
        static void syncThreads();
        static void syncComputeThreads();
        static void syncTeamThreads();
        static void syncAllThreads();
        static void syncAllComputeThreads();
        static void syncWorkerThreads();
        static void syncWorkerThreads(int numthreads);
        static void syncComputeWorkerThreads();
        static void syncComputeWorkerThreads(int numthreads);
    };

    static int tid(){return WorkerThread::perilla_tid();}
    static int wtid(){return WorkerThread::perilla_wtid();}
    static int wid(){return WorkerThread::perilla_wid();}
    static int nWorkerThreads(){return WorkerThread::perilla_nWorkerThreads();}
    static int nWorkers(){return WorkerThread::perilla_nWorkers();}
    static int nThreads(){return perilla::NUM_THREAD_TEAMS*perilla::NUM_THREADS_PER_TEAM;}
    static bool isMasterWorkerThread(){return WorkerThread::perilla_isMasterWorkerThread();}
    static bool isMasterThread(){return WorkerThread::perilla_isMasterThread();}
    static bool isCommunicationThread(){return WorkerThread::perilla_isCommunicationThread();}
    static bool isMyRegion(int workerID, int regionID){return WorkerThread::isMyRegion(workerID, regionID);}
    static void setTeamSharedMemory(void* dummy, int tid, int tg){WorkerThread::setTeamSharedMemory(dummy, tid, tg);}
    static void* getTeamSharedMemory(int tg){WorkerThread::getTeamSharedMemory(tg);}
    static void syncWorkers(){WorkerThread::syncWorkers();}
    static void syncThreads(){WorkerThread::syncThreads();}
    static void syncComputeThreads(){WorkerThread::syncComputeThreads();}
    static void syncComputeWorkerThreads(){WorkerThread::syncComputeWorkerThreads();}
    static void syncComputeWorkerThreads(int numthreads){WorkerThread::syncComputeWorkerThreads(numthreads);}
    static void syncWorkerThreads(){WorkerThread::syncWorkerThreads();}
    static void syncWorkerThreads(int numthreads){WorkerThread::syncWorkerThreads(numthreads);}
    static void syncAllWorkerThreads(){
        perilla::syncWorkerThreads();
        perilla::syncWorkers();
    }
    static void syncAllComputeThreads(){
        perilla::syncAllComputeThreads();	
    }
    static void syncAllThreads(){
        WorkerThread::syncAllThreads();
    }
}//end namespace

#endif
