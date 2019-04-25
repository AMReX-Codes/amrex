#include <WorkerThread.H>
#include <PerillaConfig.H>
#include <Perilla.H>
#include<stdio.h>

namespace perilla
{
    void* WorkerThread::team_shared_memory[perilla::NUM_THREAD_TEAMS];
    Barrier* WorkerThread::globalBarrier;
    Barrier WorkerThread::localBarriers[perilla::NUM_THREAD_TEAMS];

    void WorkerThread::init(){
        WorkerThread::globalBarrier= new Barrier(perilla::NUM_THREAD_TEAMS);
    }

    void WorkerThread::syncWorkers(){
	if(isMasterWorkerThread()) WorkerThread::globalBarrier->sync(perilla::NUM_THREAD_TEAMS);
    }

    void WorkerThread::syncTeamThreads(){
        WorkerThread::localBarriers[perilla_wid()].sync(perilla::NUM_THREADS_PER_TEAM);
    }

    void WorkerThread::syncWorkerThreads(){
        WorkerThread::localBarriers[perilla_wid()].sync(perilla::NUM_THREADS_PER_TEAM-1);
    }
    void WorkerThread::syncWorkerThreads(int numthreads){
        WorkerThread::localBarriers[perilla_wid()].sync(numthreads);
    }

    void WorkerThread::syncAllComputeThreads(){
        syncWorkerThreads();
        syncWorkers();
    }

    void WorkerThread::syncAllThreads(){
	#pragma omp barrier
    }

    void WorkerThread::syncThreads(){
        syncWorkerThreads();
        syncWorkers;
    }

    void WorkerThread::syncComputeWorkerThreads(){
        WorkerThread::localBarriers[perilla_wid()].sync(perilla::NUM_THREADS_PER_TEAM-1);
    }

    void WorkerThread::syncComputeWorkerThreads(int numthreads){
        WorkerThread::localBarriers[perilla_wid()].sync(numthreads);
    }
 
    int WorkerThread::perilla_tid(){
	return omp_get_thread_num();
    }

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
	int tid= omp_get_thread_num();
	return (tid % perilla::NUM_THREADS_PER_TEAM) -1;    
    }

    int WorkerThread::perilla_wid()
    {
	int tid= omp_get_thread_num();
	return tid / perilla::NUM_THREADS_PER_TEAM;    
    }

    bool WorkerThread::perilla_isMasterWorkerThread()
    {
	int tid= omp_get_thread_num();
	if((tid % perilla::NUM_THREADS_PER_TEAM)==1)
	    return true;
	else
	    return false;
    }

    bool WorkerThread::perilla_isMasterThread(){
	return perilla_tid()==1;
    }

    bool WorkerThread::perilla_isCommunicationThread()
    {
	int tid= omp_get_thread_num();
	return (tid % perilla::NUM_THREADS_PER_TEAM)==0 ;
    }

    bool WorkerThread::isMyRegion(int workerID, int regionID)
    {
	return ((regionID) % perilla::NUM_THREAD_TEAMS)==workerID;
    }

    void WorkerThread::setTeamSharedMemory(void* dummy, int tid, int tg)
    {
	if((tid % perilla::NUM_THREADS_PER_TEAM)==1)
	    team_shared_memory[tg] = dummy;    
    }

    void* WorkerThread::getTeamSharedMemory(int tg)
    {    
	return team_shared_memory[tg];
    }
}//end namepsace
