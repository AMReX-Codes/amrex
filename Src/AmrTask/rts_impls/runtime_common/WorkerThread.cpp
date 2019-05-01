#include <WorkerThread.H>
#include <PerillaConfig.H>
#include <Perilla.H>
#include<stdio.h>

using namespace perilla;
//namespace perilla
//{
    //void* WorkerThread::team_shared_memory[perilla::NUM_THREAD_TEAMS];
    Barrier* WorkerThread::globalBarrier;
    Barrier* WorkerThread::localBarriers[perilla::NUM_THREAD_TEAMS];
    Barrier* WorkerThread::localBarriers1[perilla::NUM_THREAD_TEAMS];

    void WorkerThread::init(){
	WorkerThread::globalBarrier= new Barrier(perilla::NUM_THREAD_TEAMS);
	for(int i=0; i<perilla::NUM_THREAD_TEAMS; i++) WorkerThread::localBarriers[i]= new Barrier(perilla::NUM_THREADS_PER_TEAM);
	for(int i=0; i<perilla::NUM_THREAD_TEAMS; i++) WorkerThread::localBarriers1[i]= new Barrier(perilla::NUM_THREADS_PER_TEAM);
    }

    void WorkerThread::syncWorkers(){
	if(isMasterWorkerThread()) WorkerThread::globalBarrier->sync(perilla::NUM_THREAD_TEAMS);
    }

    void WorkerThread::syncTeamThreads(){
        WorkerThread::localBarriers[perilla_wid()]->sync(perilla::NUM_THREADS_PER_TEAM);
    }

    void WorkerThread::syncWorkerThreads(){
        WorkerThread::localBarriers1[perilla_wid()]->sync(perilla::NUM_THREADS_PER_TEAM-1);
    }
    void WorkerThread::syncWorkerThreads(int numthreads){
	assert(numthreads== perilla::NUM_THREADS_PER_TEAM-1);
        WorkerThread::localBarriers1[perilla_wid()]->sync(numthreads);
    }


#ifdef USE_PERILLA_OMP
    void WorkerThread::syncAllThreads(){
        #pragma omp barrier
    }
#else
    void WorkerThread::syncAllThreads(){
        syncTeamThreads();
        syncWorkers();
    }
#endif

    void WorkerThread::syncAllComputeThreads(){
        syncComputeWorkerThreads();
        syncWorkers();
    }

    void WorkerThread::syncThreads(){
        syncWorkerThreads();
        syncWorkers;
    }

    void WorkerThread::syncComputeWorkerThreads(){
        WorkerThread::localBarriers1[perilla_wid()]->sync(perilla::NUM_THREADS_PER_TEAM-1);
    }

    void WorkerThread::syncComputeWorkerThreads(int numthreads){
	assert(numthreads== perilla::NUM_THREADS_PER_TEAM-1);
        WorkerThread::localBarriers1[perilla_wid()]->sync(numthreads);
    }

#ifdef USE_PERILLA_OMP 
    int WorkerThread::perilla_tid(){
        return omp_get_thread_num();
    }
#else
    int WorkerThread::perilla_tid(){
	return Perilla::tid();
    }
#endif

    int WorkerThread::perilla_nTeamThreads(){
	return perilla::NUM_THREADS_PER_TEAM;
    }

    int WorkerThread::perilla_nWorkerThreads(){
	return perilla::NUM_THREADS_PER_TEAM-1;
    }

    int WorkerThread::perilla_nWorkers(){
	return perilla::NUM_THREAD_TEAMS;
    }

    int WorkerThread::perilla_wtid()
    {
	int tid= perilla_tid();
	return (tid % perilla::NUM_THREADS_PER_TEAM) -1;    
    }

    int WorkerThread::perilla_wid()
    {
	int tid= perilla_tid();
	return tid / perilla::NUM_THREADS_PER_TEAM;    
    }

    bool WorkerThread::perilla_isMasterWorkerThread()
    {
	int tid= perilla_tid();
	if((tid % perilla::NUM_THREADS_PER_TEAM)==1)
	    return true;
	else
	    return false;
    }

    bool WorkerThread::perilla_isMasterThread(){ //pick the first one among master worker threads
	return perilla_tid()==1;
    }

    bool WorkerThread::perilla_isCommunicationThread()
    {
	int tid= perilla_tid();
	return (tid % perilla::NUM_THREADS_PER_TEAM)==0 ;
    }

    bool WorkerThread::isMyRegion(int workerID, int regionID)
    {
	return ((regionID) % perilla::NUM_THREAD_TEAMS)==workerID;
    }

#if 0
    void WorkerThread::setTeamSharedMemory(void* dummy, int tid, int tg)
    {
	if((tid % perilla::NUM_THREADS_PER_TEAM)==1)
	    team_shared_memory[tg] = dummy;    
    }

    void* WorkerThread::getTeamSharedMemory(int tg)
    {    
	return team_shared_memory[tg];
    }
#endif
//}//end namepsace
