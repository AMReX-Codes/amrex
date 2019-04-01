#include <RegionGraph.H>
#include <WorkerThread.H>
#include <stdlib.h>
#include <omp.h>

using namespace std;
using namespace amrex;
using namespace perilla;

int RegionGraph::graphCnt = 0;

RegionGraph::RegionGraph(int numtasks)
{
    sCopyMapHead = 0;
    rCopyMapHead = 0;
    srcLinkGraph = 0;
    isDepGraph = false;
    numFabs = numtasks;
    numTasks = numtasks;
    graphID = ++graphCnt;
    worker.resize(perilla::NUM_THREAD_TEAMS);
    task.resize(numTasks);
    totalFinishes=0;
    okToReset = new bool[perilla::NUM_THREAD_TEAMS];
    omp_init_lock(&finishLock);
    Initialize();
#ifdef PERILLA_DEBUG
    memcheck.add(memcheck.genKey(this), (void*)this, "Package");
#endif
}

void RegionGraph::Initialize()
{
    int numfabs = numTasks;        
    int numthreads = omp_get_num_threads();

    if(numthreads==1)
    {
#pragma omp parallel shared(numfabs)
	{
	    int tg = WorkerThread::perilla_wid();

	    if(WorkerThread::perilla_isMasterWorkerThread())
	    {
		worker[tg] = new Worker();
		worker[tg]->barr = new Barrier(perilla::NUM_THREADS_PER_TEAM-1);
		worker[tg]->l_barr = new Barrier(perilla::NUM_THREADS_PER_TEAM-2);
		if(numfabs <= perilla::TASK_QUEUE_DEFAULT_MAXSIZE)
		{
		    worker[tg]->fireableRegionQueue = new RegionQueue();
		    worker[tg]->unfireableRegionQueue = new RegionQueue();
		    worker[tg]->computedRegionQueue = new RegionQueue();
		    worker[tg]->completedRegionQueue = new RegionQueue();
		}
		else
		{
		    worker[tg]->fireableRegionQueue = new RegionQueue(numfabs);
		    worker[tg]->unfireableRegionQueue = new RegionQueue(numfabs);
		    worker[tg]->computedRegionQueue = new RegionQueue(numfabs);
		    worker[tg]->completedRegionQueue = new RegionQueue(numfabs);
		}
		worker[tg]->totalTasks = 0;
		worker[tg]->computedTasks = 0;
		for(int f=0; f < numfabs; f++)
		    if(WorkerThread::isMyRegion(tg,f))
		    {
			task[f] = new Task();
			worker[tg]->unfireableRegionQueue->addRegion(f);
			worker[tg]->totalTasks++;
			for(int i=0; i<16; i++)
			    task[f]->state[i] = 0;
			task[f]->init = true;
		    }
		worker[tg]->init = true;
		okToReset[tg] = false;	      
	    }
	}// omp parallel end
    }
    else // numthread are > 1, so already in parallel region
    {
	int tg = WorkerThread::perilla_wid();
	if(WorkerThread::perilla_isMasterWorkerThread() && worker[tg]->init == false )
	{
	    worker[tg]->barr = new Barrier(perilla::NUM_THREADS_PER_TEAM-1);
	    worker[tg]->l_barr = new Barrier(perilla::NUM_THREADS_PER_TEAM-2);
	    worker[tg]->fireableRegionQueue = new RegionQueue();
	    worker[tg]->unfireableRegionQueue = new RegionQueue();
	    worker[tg]->completedRegionQueue = new RegionQueue();
	    worker[tg]->totalTasks = 0;
	    worker[tg]->computedTasks = 0;
	    for(int f=0; f < numfabs; f++)
		if(WorkerThread::isMyRegion(tg,f))
		{
		    worker[tg]->unfireableRegionQueue->addRegion(f);
		    worker[tg]->totalTasks++;
		    for(int i=0; i<16; i++)
			task[f]->state[i] = 0;
		    task[f]->init = true;
		}
	    worker[tg]->init = true;
	}	
    }
}

void RegionGraph::Reset()
{
    int tg= perilla::wid();
    omp_set_lock(&finishLock);
    if(okToReset[tg])
	totalFinishes--;
    omp_unset_lock(&finishLock);

    if(okToReset[tg])
    {
	worker[tg]->totalTasks = 0;
	worker[tg]->computedTasks = 0;
	while(worker[tg]->completedRegionQueue->queueSize(true) > 0)
	{
	    int r = worker[tg]->completedRegionQueue->removeRegion(true);
	    if(WorkerThread::isMyRegion(tg, r))
	    {
		worker[tg]->unfireableRegionQueue->addRegion(r,true);
		worker[tg]->totalTasks++;
		for(int i=0; i<16; i++)
		    task[r]->state[i] = 0;
		task[r]->init = true;
		if(task[r]->depTaskIDs.size() > 0)
		    task[r]->depTasksCompleted = false; 
	    }
	    else
		break;
	}
    }
}

bool RegionGraph::isGraphEmpty()
{
    int tg= perilla::wid();
    //worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); 
    perilla::syncWorkerThreads();
    if(worker[tg]->completedRegionQueue->queueSize(true)== worker[tg]->totalTasks)
	return true;
    return false;	       
}

bool RegionGraph::isGraphEmptyV2()
{
    int tg=perilla::wid();

    if(worker[tg]->completedRegionQueue->queueSize(true) == worker[tg]->totalTasks || worker[tg]->computedTasks == worker[tg]->totalTasks)
	return true;
	return false;	       
}

void RegionGraph::finalizeGraph()
{
    omp_set_lock(&finishLock);
    totalFinishes++;
    int tg=perilla::wid();
    okToReset[tg]=true;
    omp_unset_lock(&finishLock);
}

void RegionGraph::regionGraphReset(int numfabs)
{
    int nt;
    int tg;
    int r;
    //#pragma omp parallel private(r,tg,nt,tid) shared(numfabs)
    {
        tg = perilla::wid();
        nt = perilla::wtid();
        if(perilla::isMasterThread())
	    totalFinishes=0;	
	//#pragma omp barrier
        if(perilla::isMasterWorkerThread())
	{
	    worker[tg]->totalTasks = 0;
	    worker[tg]->computedTasks = 0;
	    while(worker[tg]->completedRegionQueue->queueSize(true) > 0)
	    {
		r = worker[tg]->completedRegionQueue->removeRegion(true);
		if(WorkerThread::isMyRegion(tg, r))
		{
		    worker[tg]->unfireableRegionQueue->addRegion(r,true);
		    worker[tg]->totalTasks++;
		    for(int i=0; i<16; i++)
			task[r]->state[i] = 0;
		    task[r]->init = true;
		}
		else
		    break;
	    }
	    okToReset[tg] = false;
	}
    }// omp parallel end
}


void RegionGraph::regionGraphMinReset(void)
{
    int nt;
    int tg;
    int r;
    {
        tg = perilla::wid();
        nt = perilla::wtid();
	if(perilla::isMasterThread())	
	    totalFinishes=0;	
	if(perilla::isMasterWorkerThread())
	{
	    while(worker[tg]->completedRegionQueue->queueSize(true) > 0)
	    {
		r = worker[tg]->completedRegionQueue->removeRegion(true);
		if(WorkerThread::isMyRegion(tg, r))
		{
		    worker[tg]->unfireableRegionQueue->addRegion(r,true);
		}
		else
		    break;
	    }
	    okToReset[tg] = false;
	}
    }
}


void RegionGraph::enableAllRegions()
{
    int numfabs = numTasks;
    int r;
    int tg = WorkerThread::perilla_wid();
    perilla::syncWorkerThreads();
    if(perilla::isMasterWorkerThread())
	for(int f=0; f<numfabs; f++)
	    if(WorkerThread::isMyRegion(tg, f))
	    {
		r = worker[tg]->unfireableRegionQueue->removeRegion(true);
		worker[tg]->fireableRegionQueue->addRegion(r,true);
	    }    
    //worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads        
    perilla::syncWorkerThreads();
}

void RegionGraph::disableRegion(int r, int tg)
{
    //int tg = perilla::wid();
    if(perilla::isMasterWorkerThread())
	if(WorkerThread::isMyRegion(tg, r))
	{
	    int rID = worker[tg]->fireableRegionQueue->removeRegion(true);
	    worker[tg]->unfireableRegionQueue->addRegion(rID,true);
	}
}

void RegionGraph::regionComputed(int r)
{
    int tg= perilla::wid();
    worker[tg]->l_barr->sync(perilla::NUM_THREADS_PER_TEAM-2);
    if(perilla::isMasterWorkerThread())
	if(WorkerThread::isMyRegion(tg, r))
	{
	    int rr = worker[tg]->fireableRegionQueue->removeRegion();
	    if(r != rr)
	    {
		std::cout << "ERROR: In computedeRegion" << std::endl;
		exit(EXIT_FAILURE);
	    }
	    worker[tg]->computedRegionQueue->addRegion(rr);
	    worker[tg]->computedTasks++;
	}
    worker[tg]->l_barr->sync(perilla::NUM_THREADS_PER_TEAM-2);
}

void RegionGraph::finalizeRegion(int r)
{
    int tg= perilla::wid();
    int ntid=perilla::wtid();
    perilla::syncWorkerThreads();
    if(perilla::isMasterWorkerThread())
	if(WorkerThread::isMyRegion(tg, r))
	{
	    int rr = worker[tg]->fireableRegionQueue->removeRegion(true);
	    if(r != rr)
	    {
		std::cout << "ERROR: In completeRegion" << std::endl;
		exit(EXIT_FAILURE);
	    }
	    worker[tg]->completedRegionQueue->addRegion(rr,true);
	}
    perilla::syncWorkerThreads();
}

void RegionGraph::finalizeRegionGraph()
{
    int tg= perilla::wid();
    omp_set_lock(&finishLock);
    totalFinishes++;
    okToReset[tg]=true;
    omp_unset_lock(&finishLock);
}

bool RegionGraph::isFireableRegion(int r)
{
    int myProc = ParallelDescriptor::MyProc();
    FabCopyAssoc *cpDst = task[r]->cpAsc_dstHead;
    if(lMap.size() > 0)
	if(lMap[r]->l_con.firingRuleCnt != lMap[r]->l_con.ndcpy)
	{
	    return false;
	}
	while(cpDst != 0)
	{
	    if(cpDst->l_con.firingRuleCnt != cpDst->l_con.ndcpy)
	    {
		return false;
	    }
	    cpDst = cpDst->next;
	}

	if(srcLinkGraph != 0)
	{
	    if(!task[r]->depTasksCompleted)
	    {
		for(int i=0; i<task[r]->depTaskIDs.size(); i++)
		    if(!srcLinkGraph->isFireableRegion(task[r]->depTaskIDs[i]))
			return false;
		task[r]->depTasksCompleted = true;
	    }
	}

    if(ParallelDescriptor::NProcs() == 1) return true;

    if(lMap.size() > 0)
	if(lMap[r]->r_con.firingRuleCnt != lMap[r]->r_con.nrcv)
	{
	    return false;
	}

	cpDst = task[r]->cpAsc_dstHead;
	while(cpDst != 0)
	{
	    if(cpDst->r_con.firingRuleCnt != cpDst->r_con.nrcv)
	    {
		return false;
	    }
	    cpDst = cpDst->next;
	}
    return true;
}

int RegionGraph::getFireableRegion(bool isSingleThread)
{
    int r = -1;
    bool fireable;
    int tg= perilla::wid();

    if(worker[tg]->unfireableRegionQueue->queueSize(true)!=0 && worker[tg]->fireableRegionQueue->queueSize() == 0)
    {
	fireable = false;
	r = worker[tg]->unfireableRegionQueue->removeRegion(true);
	while(!fireable)
	{
	    fireable = isFireableRegion(r);
	    if(!fireable)
	    {
		worker[tg]->unfireableRegionQueue->addRegion(r,true);
		r = worker[tg]->unfireableRegionQueue->removeRegion(true);
	    }
	}
    }
    else if(worker[tg]->unfireableRegionQueue->queueSize(true)!=0)
    {
	int unfQsize = worker[tg]->unfireableRegionQueue->queueSize(true);
	for(int i = 0; i < unfQsize; i++)
	{
	    int tr = worker[tg]->unfireableRegionQueue->removeRegion(true);
	    if(isFireableRegion(tr))
	    {
		r = tr;
		break;
	    }
	    else
		worker[tg]->unfireableRegionQueue->addRegion(tr,true);
	}
    }

    return r;
}

#if 0
int RegionGraph::getFireableRegion(bool patchFilled, bool isSingleThread)
{   
    int r = -1;
    bool fireable;
    int tg= perilla::wid();
    int nt= perilla::wtid();

    //if(worker[tg]->unfireableRegionQueue->queueSize(true)!=0 && worker[tg]->fireableRegionQueue->queueSize() == 0)      
    //{
    if(!isSingleThread)worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads
    if(perilla::isMasterWorkerThread())
    { 
	if(worker[tg]->fireableRegionQueue->queueSize()==0){
	    fireable = false;
	    assert(worker[tg]->unfireableRegionQueue->queueSize()>0);
	    r = worker[tg]->unfireableRegionQueue->removeRegion(true);
	    while(!fireable)
	    {   
		fireable = isFireableRegion(r, patchFilled);
		//fireable = true;
		if(!fireable)
		{   
		    worker[tg]->unfireableRegionQueue->addRegion(r,true);
		    r = worker[tg]->unfireableRegionQueue->removeRegion(true);
		}
		else worker[tg]->fireableRegionQueue->addRegion(r,true);
	    }
	}
    }
#if 0
    else if(worker[tg]->unfireableRegionQueue->queueSize(true)!=0)
    {
	int unfQsize = worker[tg]->unfireableRegionQueue->queueSize(true);
	for(int i = 0; i < unfQsize; i++)
	{
	    int tr = worker[tg]->unfireableRegionQueue->removeRegion(true);
	    if(isFireableRegion(tr))
	    {
		r = tr;
		break;
	    }
	    else
		worker[tg]->unfireableRegionQueue->addRegion(tr,true);
	}
    }
#endif
    if(!isSingleThread)worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads
    std::cout<<"FOUND A REGION"<<r<<std::endl;
    r = worker[tg]->fireableRegionQueue->getFrontRegion(true);
    return r;
}
#endif



#if 0
int RegionGraph::getFireableRegion(bool patchFilled)
{
    int r = -1;
    bool fireable;
    int tg= perilla::wid();
    int nt= perilla::wtid();

    //if(worker[tg]->unfireableRegionQueue->queueSize(true)!=0 && worker[tg]->fireableRegionQueue->queueSize() == 0)      
    //{
    worker[tg]->barr->sync(); // Barrier to synchronize team threads
    if(nt == 0 && worker[tg]->fireableRegionQueue->queueSize()==0){      
	fireable = false;
	assert(worker[tg]->unfireableRegionQueue->queueSize()>0);
	r = worker[tg]->unfireableRegionQueue->removeRegion(true);
	while(!fireable)
	{
	    fireable = isFireableRegion(r, patchFilled);
	    //fireable = true;
	    if(!fireable)
	    {
		worker[tg]->unfireableRegionQueue->addRegion(r,true);
		r = worker[tg]->unfireableRegionQueue->removeRegion(true);
	    }
	    else worker[tg]->fireableRegionQueue->addRegion(r,true);
	}	  
    }
#if 0
    else if(worker[tg]->unfireableRegionQueue->queueSize(true)!=0)
    {
	int unfQsize = worker[tg]->unfireableRegionQueue->queueSize(true);
	for(int i = 0; i < unfQsize; i++)
	{
	    int tr = worker[tg]->unfireableRegionQueue->removeRegion(true);
	    if(isFireableRegion(tr))
	    {
		r = tr;
		break;
	    }
	    else
		worker[tg]->unfireableRegionQueue->addRegion(tr,true);
	}
    }
#endif
    worker[tg]->barr->sync(); // Barrier to synchronize team threads
    r = worker[tg]->fireableRegionQueue->getFrontRegion(true);
    return r;
}
#endif

void RegionGraph::setFireableRegion(int r)
{
    worker[perilla::wid()]->fireableRegionQueue->addRegion(r);
}


int RegionGraph::getAnyFireableRegion()
{
    int myProc = ParallelDescriptor::MyProc();
    int tg = perilla::wid();
    int nt = perilla::wtid();
    int r;
    perilla::syncWorkerThreads();
    if(nt ==0)
    if(worker[tg]->fireableRegionQueue->queueSize()==0)      
    {
	bool fireable = false;
	r = worker[tg]->unfireableRegionQueue->removeRegion(true);
	while(!fireable)
	{
	    fireable = isFireableRegion(r);
	    if(!fireable)
	    {
		worker[tg]->unfireableRegionQueue->addRegion(r,true);
		r = worker[tg]->unfireableRegionQueue->removeRegion(true);
	    }
	    else
		worker[tg]->fireableRegionQueue->addRegion(r,true);
	}
    }
    perilla::syncWorkerThreads();
    return worker[tg]->fireableRegionQueue->getFrontRegion(true);
}

int RegionGraph::getAnyFireableRegion(RegionGraph& depGraph)
{
    int nt;
    int tg;
    int r;
    bool fireable;

    int myProc = amrex::ParallelDescriptor::MyProc();

    tg = perilla::wid();
    nt = perilla::wtid();
    if(nt == perilla::NUM_COMM_THREADS && worker[tg]->fireableRegionQueue->queueSize()==0)
    {
        fireable = false;
        r = worker[tg]->unfireableRegionQueue->removeRegion(true);
        while(!fireable)
        {
            fireable = isFireableRegion(r);
            fireable &= depGraph.isFireableRegion(r);
            if(!fireable)
            {
                worker[tg]->unfireableRegionQueue->addRegion(r,true);
                r = worker[tg]->unfireableRegionQueue->removeRegion(true);
            }
            else
                worker[tg]->fireableRegionQueue->addRegion(r,true);
        }
    }
    worker[tg]->barr->sync();
    //worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS); // Barrier to synchronize team threads
    r = worker[tg]->fireableRegionQueue->getFrontRegion(true);
    return r;
}


int RegionGraph::getPulledFireableRegion()
{
    int myProc = ParallelDescriptor::MyProc();
    int tg = WorkerThread::perilla_wid();
    int nt = WorkerThread::perilla_wtid();
    if(nt == 0 && worker[tg]->fireableRegionQueue->queueSize()==0)      
    {
	while(worker[tg]->fireableRegionQueue->queueSize()==0);
    }
    worker[tg]->l_barr->sync(perilla::NUM_THREADS_PER_TEAM-2);
    return worker[tg]->fireableRegionQueue->getFrontRegion(true);
}

void RegionGraph::graphTeardown()
{
    MPI_Status status;
    Package* package;
    int numfabs = numTasks;
    int tg = WorkerThread::perilla_wid();

#if 0
    for(int f=0; f<numfabs; f++)
    {
	if(WorkerThread::isMyRegion(tg,f))
	{
	    FabCopyAssoc *cpDst = task[f]->cpAsc_dstHead;
	    while(cpDst != 0)
	    {
		cpDst->l_con.firingRuleCnt = 0;

		for(int i=0; i<cpDst->l_con.ndcpy; i++)
		{
		    while(cpDst->l_con.dcpy[i].pQueue.queueSize() >= 1)
		    {
			package = cpDst->l_con.dcpy[i].pQueue.dequeue();
			//package->completed = false;
			//package->served = false;
			//package->notified = false;
			//package->request = MPI_REQUEST_NULL;		    
			cpDst->l_con.dcpy[i].recycleQueue.enqueue(package);
		    }
		}

		cpDst = cpDst->next;
	    }
	}
    }


    for(int f=0; f<numfabs; f++)
    {
	if(WorkerThread::isMyRegion(tg,f))
	{
	    FabCopyAssoc *cpSrc = task[f]->cpAsc_srcHead;
	    while(cpSrc != 0)
	    {
		//cpSrc->l_con.firingRuleCnt = 0;

		for(int i=0; i<cpSrc->l_con.nscpy; i++)
		{
		    while(cpSrc->l_con.scpy[i].pQueue.queueSize() >= 1)
		    {
			package = cpSrc->l_con.scpy[i].pQueue.dequeue();

			FabCopyAssoc* cpDst = cpSrc->graphPartner->task[cpSrc->l_con.scpy[i].nd]->cpAsc_dstHead;
			while(cpDst != 0)
			{
			    if(cpDst->graphPartner == this) //graphArray[g])
				break;
			    cpDst = cpDst->next;
			}			    
			//Package* sPackage = cpSrc->l_con.scpy[i].pQueue.dequeue(true);
			omp_set_lock(&(cpDst->l_con.dLock));
			int dPartner = cpSrc->l_con.scpy[i].dPartner;
			Package* dPackage = cpDst->l_con.dcpy[dPartner].recycleQueue.dequeue(true);
			/*
			   for(int j=0; j<dPackage->bufSize; j++)
			   {
			   dPackage->databuf[j] = sPackage->databuf[j];
			   }
			 */
			std::memcpy(dPackage->databuf, package->databuf, dPackage->bufSize * sizeof(double));
			//std::swap(dPackage->databuf, sPackage->databuf);

			cpDst->l_con.dcpy[dPartner].pQueue.enqueue(dPackage,true);
			if(cpDst->l_con.dcpy[dPartner].pQueue.queueSize(true) == 1)
			    cpDst->l_con.firingRuleCnt++;
			omp_unset_lock(&(cpDst->l_con.dLock));
			//cpSrc->l_con.scpy[i].recycleQueue.enqueue(sPackage,true);

			//package->completed = false;
			//package->served = false;
			//package->notified = false;
			//package->request = MPI_REQUEST_NULL;		    
			cpSrc->l_con.scpy[i].recycleQueue.enqueue(package);
		    }
		}

		cpSrc = cpSrc->next;
	    }
	}
    }



    for(int f=0; f<numfabs; f++)
    {
	if(WorkerThread::isMyRegion(tg,f))
	{
	    if(lMap.size() > 0)
	    {
		lMap[f]->l_con.firingRuleCnt = 0;
	    }
	}
    }
#endif

    if(ParallelDescriptor::NProcs() == 1) return;



    for(int f=0; f<numfabs; f++)
    {
	if(WorkerThread::isMyRegion(tg,f))
	{
	    FabCopyAssoc *cpDst = task[f]->cpAsc_dstHead;
	    while(cpDst != 0)
	    {
		cpDst->r_con.firingRuleCnt = 0;
		for(int i=0; i<cpDst->r_con.nrcv; i++)
		{
		    while(cpDst->r_con.rcv[i].pQueue.queueSize() >= 1)
		    {
			package = cpDst->r_con.rcv[i].pQueue.dequeue();
			package->completed = false;
			package->served = false;
			package->notified = false;
			package->request = MPI_REQUEST_NULL;		    
			cpDst->r_con.rcv[i].recycleQueue.enqueue(package);
		    }
		}

		cpDst = cpDst->next;
	    }
	}
    }


    for(int f=0; f<numfabs; f++)
    {
	if(WorkerThread::isMyRegion(tg,f))
	{
	    FabCopyAssoc *cpSrc = task[f]->cpAsc_srcHead;
	    while(cpSrc != 0)
	    {
		//cpSrc->r_con.firingRuleCnt = 0;
		for(int i=0; i<cpSrc->r_con.nsnd; i++)
		{
		    while(cpSrc->r_con.snd[i].pQueue.queueSize() >= 1)
		    {
			package = cpSrc->r_con.snd[i].pQueue.dequeue();
			package->completed = false;
			package->served = false;
			package->notified = false;
			package->request = MPI_REQUEST_NULL;		    
			cpSrc->r_con.snd[i].recycleQueue.enqueue(package);
		    }
		}

		cpSrc = cpSrc->next;
	    }
	}
    }


#if 0
    if(tg == 0)
    {
	CopyMap* cpDst = rCopyMapHead;
	while(cpDst != 0)
	{
	    for(int f=0; f<cpDst->map.size(); f++)
	    {
		cpDst->map[f]->r_con.firingRuleCnt = 0;
		for(int i=0; i<cpDst->map[f]->r_con.nrcv; i++)
		{
		    while(cpDst->map[f]->r_con.rcv[i].pQueue.queueSize() >= 1)
		    {
			package = cpDst->map[f]->r_con.rcv[i].pQueue.dequeue();
			if(package->request != MPI_REQUEST_NULL)
			    MPI_Cancel( &(package->request) );
			package->completed = false;
			package->served = false;
			package->notified = false;
			package->request = MPI_REQUEST_NULL;		    
			cpDst->map[f]->r_con.rcv[i].recycleQueue.enqueue(package);
		    }
		}

	    }

	    cpDst = cpDst->next;
	}


	CopyMap* cpSrc = sCopyMapHead;
	while(cpSrc != 0)
	{
	    for(int f=0; f<cpSrc->map.size(); f++)
	    {
		for(int i=0; i<cpSrc->map[f]->r_con.nsnd; i++)
		{
		    while(cpSrc->map[f]->r_con.snd[i].pQueue.queueSize() >= 1)
		    {

			package = cpSrc->map[f]->r_con.snd[i].pQueue.dequeue();
			/*		
					int ns = cpSrc->map[f]->r_con.snd[i].ns;
					int nd = cpSrc->map[f]->r_con.snd[i].nd;
					int r_gid = cpSrc->map[f]->r_con.snd[i].r_gid;
					int r_grids = cpSrc->map[f]->r_con.snd[i].r_grids;
			//int tag = tagGen(ns, nd, r_gid-1, np*r_grids, nGraphs);
			int tag = Perilla::myTagMap[r_gid][nd][ns][cpSrc->map[f]->r_con.snd[i].sz];

			Package* sPackage = lMap[f]->r_con.snd[i].pQueue.getFront(true);
			package->request = ParallelDescriptor::Asend(sPackage->databuf,
			cpSrc->map[f]->r_con.snd[i].sz,
			cpSrc->map[f]->r_con.snd[i].pr, tag).req();  // tag == SeqNum in c++ ver

			 */
			MPI_Wait( &(package->request), &status );
			package->completed = false;
			package->served = false;
			package->notified = false;
			package->request = MPI_REQUEST_NULL;		    
			cpSrc->map[f]->r_con.snd[i].recycleQueue.enqueue(package);
		    }
		}

	    }

	    cpSrc = cpSrc->next;
	}
    }

    //if(WorkerThread::isTeamMasterThread(tid)) commented out b/c its already call by single thread in a team
    //Perilla::globalBarrier->sync(perilla::NUM_THREAD_TEAMS);

    // Parallel Copy Reset on Local tg
    for(int f=0; f<numfabs; f++)
    {
	//if(WorkerThread::isMyRegion(tg,f))
	{
	    if(lMap.size() > 0)
	    {
		lMap[f]->r_con.firingRuleCnt = 0;

		for(int i=0; i<lMap[f]->r_con.nsnd; i++)
		    while(lMap[f]->r_con.snd[i].pQueue.queueSize() >= 1)
		    {
			package = lMap[f]->r_con.snd[i].pQueue.dequeue();
			package->completed = false;
			package->served = false;
			package->notified = false;
			package->request = MPI_REQUEST_NULL;
			lMap[f]->r_con.snd[i].recycleQueue.enqueue(package);
		    }

		for(int i=0; i<lMap[f]->r_con.nrcv; i++)
		    while(lMap[f]->r_con.rcv[i].pQueue.queueSize() >= 1)
		    {
			package = lMap[f]->r_con.rcv[i].pQueue.dequeue();
			package->completed = false;
			package->served = false;
			package->notified = false;
			package->request = MPI_REQUEST_NULL;
			lMap[f]->r_con.rcv[i].recycleQueue.enqueue(package);
		    }
	    }
	}
    }              																             

    // Fill boundary reset on local tg
    if(tg == 0)
    {
	for(int f=0; f<numfabs; f++)
	{
	    if(rMap.size() > 0)
	    {
		// if(WorkerThread::isMyRegion(tg,f))
		{
		    for(int i=0; i< rMap[f]->r_con.nrcv; i++)
			while( rMap[f]->r_con.rcv[i].pQueue.queueSize() >= 1)
			{
			    package =  rMap[f]->r_con.rcv[i].pQueue.dequeue();
			    if(package->request != MPI_REQUEST_NULL)
				MPI_Cancel( &(package->request) );
			    package->completed = false;
			    package->served = false;
			    package->notified = false;
			    package->request = MPI_REQUEST_NULL;
			    rMap[f]->r_con.rcv[i].recycleQueue.enqueue(package);
			}
		    for(int i=0; i< sMap[f]->r_con.nsnd; i++)
			while( sMap[f]->r_con.snd[i].pQueue.queueSize() >= 1)
			{
			    package =  sMap[f]->r_con.snd[i].pQueue.dequeue();
			    MPI_Wait( &(package->request), &status );
			    package->completed = false;
			    package->served = false;
			    package->notified = false;
			    package->request = MPI_REQUEST_NULL;
			    sMap[f]->r_con.snd[i].recycleQueue.enqueue(package);
			}
		}
	    }
	}
    }
#endif

}

void RegionGraph::workerTeardown()
{
    int numfabs = numTasks;
    Package* package;

    regionGraphMinReset();   
}

RegionGraph::~RegionGraph()
{
    delete[] okToReset;
    for(int tg=0; tg<perilla::NUM_THREAD_TEAMS; tg++)delete worker[tg];
    worker.clear();
    for(int i=0; i<task.size(); i++) delete task[i];
    task.clear();

    if(sCopyMapHead != 0)
      delete sCopyMapHead;
    if(rCopyMapHead != 0)
      delete rCopyMapHead;

    for(int i=0; i<lMap.size(); i++) delete lMap[i];
    for(int i=0; i<sMap.size(); i++) delete sMap[i];
    for(int i=0; i<rMap.size(); i++) delete rMap[i];

    lMap.clear();
    sMap.clear();
    rMap.clear();

    for(int i=0; i<fabTiles.size(); i++) delete fabTiles[i];
    for(int i=0; i<fabTiles_gtbx.size(); i++) delete fabTiles_gtbx[i];

    fabTiles.clear();
    fabTiles_gtbx.clear();
#ifdef PERILLA_DEBUG
    memcheck.remove(memcheck.genKey(this));
#endif
}

#if 0

RegionGraph::~RegionGraph()
{
    lMap.clear();
    sMap.clear();
    rMap.clear();
    //fabTiles.clear();
    //if(sCopyMapHead != 0)
    //  delete sCopyMapHead;
    //if(rCopyMapHead != 0)
    //  delete rCopyMapHead;
    //delete[] worker;
    //delete[] task;
    worker.clear();
    task.clear();
    delete[] okToReset;
}
#endif
