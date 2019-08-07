#include <AMReX_MultiFab.H>
#include <AMReX_FabArray.H>
#include <AMReX_Periodicity.H>
#include <WorkerThread.H>
#include <PerillaConfig.H>
#include <RegionGraph.H>
#include <Barrier.H>
#include <vector>
#include <iostream>
#include <limits>
#include <exception>
#include <mpi.h>
#include <Perilla.H>
using namespace std;
using namespace amrex;
using namespace perilla;

#ifdef PERILLA_DEBUG
#include<stdlib.h>
#include <sys/time.h>

const double kMicro = 1.0e-6;
double getTime()
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

double isendDelay=0.0;
double irecvDelay=0.0;
double isendTestDelay=0.0;
double irecvTestDelay=0.0;
double localScheDelay=0.0;
#endif

void Perilla::syncProcesses(){
    MPI_Barrier(MPI_COMM_WORLD);
}

void Perilla::serviceLocalRequests(RegionGraph* rg, int tg)
{
    int numfabs = rg->lMap.size();
    for(int f=0; f<numfabs; f++)
    {
	if(WorkerThread::isMyRegion(tg,f))
	{
	    bool anyReq=false;
	    for(int i=0; i<rg->lMap[f]->l_con.nscpy; i++)
		if(rg->lMap[f]->l_con.scpy[i].pQueue.queueSize(true)>0){
		    anyReq=true;
		    break;
		}
	    if(anyReq){
		pthread_mutex_lock(&(rg->lMap[f]->l_con.sLock));
		for(int i=0; i<rg->lMap[f]->l_con.nscpy; i++){
		    if(rg->lMap[f]->l_con.scpy[i].pQueue.queueSize(true)>0)
		    {
			Package *sPackage = rg->lMap[f]->l_con.scpy[i].pQueue.dequeue(true);
			pthread_mutex_lock(&(rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dLock));
			int dPartner = rg->lMap[f]->l_con.scpy[i].dPartner;
			Package *dPackage = rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dcpy[dPartner].recycleQueue.dequeue(true);
			std::memcpy(dPackage->databuf, sPackage->databuf, dPackage->bufSize * sizeof(double));
			rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dcpy[dPartner].pQueue.enqueue(dPackage,true);
			if(rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dcpy[dPartner].pQueue.queueSize(true)==1)
			    rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.firingRuleCnt++;
			pthread_mutex_unlock(&(rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dLock));
			rg->lMap[f]->l_con.scpy[i].recycleQueue.enqueue(sPackage,true);
		    }
		}
		pthread_mutex_unlock(&(rg->lMap[f]->l_con.sLock));
	    }//if there is any local send request
	}// if my region
    }// for(f<numfabs)    
}//serviceLocalRequests

void Perilla::serviceRemoteRequests(RegionGraph* rg, int graphID, int nGraphs)
{
    bool nextsReq, nextrReq;
    int np = ParallelDescriptor::NProcs();
    int myProc = ParallelDescriptor::MyProc();
    int numfabs = rg->rMap.size();
    int tg = WorkerThread::perilla_wid();

    // !we first pre-post receive  
    for(int f=0; f<numfabs; f++)
    {
	for(int i=0; i<rg->lMap[f]->r_con.nrcv; i++)
	{
	    if(rg->rMap[f]->r_con.rcv[i].pQueue.queueSize(true) == 0) //!no message has been received or all received messages have been claimed
		nextrReq = true;
	    else
	    {
		//we buffer at most 2 packages per send task - recv task pair, but 1 must be completed before we buffer the next to allow for tag reuse
		Package *rearPackage = rg->rMap[f]->r_con.rcv[i].pQueue.getRear(true);
		if(rearPackage) 
		    if(rearPackage->completed && rg->rMap[f]->r_con.rcv[i].pQueue.queueSize(true) == 1) //!latest receive request has been completed
			nextrReq = true;
		    else //!expected message is still on the way
			nextrReq = false;
	    }
	    if(nextrReq) //!take a message from recycle pool and post a receive
	    {
		pthread_mutex_lock(&(rg->rMap[f]->r_con.rcvLock));
		pthread_mutex_lock(&(rg->lMap[f]->r_con.rcvLock));
		int ns = rg->rMap[f]->r_con.rcv[i].ns;
		int nd = rg->rMap[f]->r_con.rcv[i].nd;
		int lnd = rg->rMap[f]->r_con.rcv[i].lnd;
		int r_grids = rg->rMap[f]->r_con.rcv[i].r_grids;
		//!create a package to keep track of receive requests
		Package *rMetaPackage = rg->rMap[f]->r_con.rcv[i].recycleQueue.dequeue(true);
		//!extract a package from the recycle pool at the destination NUMA node to buffer incoming data
		Package *rPackage = rg->lMap[f]->r_con.rcv[i].recycleQueue.dequeue(true);
		int tag = tagMap[rg->rMap[f]->r_con.rcv[i].pr][graphID][nd][ns][rg->rMap[f]->r_con.rcv[i].sz];

		rMetaPackage->request = MPI_REQUEST_NULL;
		rg->lMap[f]->r_con.rcv[i].pQueue.enqueue(rPackage,true);   //!this is not done yet
		rg->rMap[f]->r_con.rcv[i].pQueue.enqueue(rMetaPackage,true);   //!this is not done yet
		rMetaPackage->request = ParallelDescriptor::Arecv(rPackage->databuf,
			rg->rMap[f]->r_con.rcv[i].sz,
			rg->rMap[f]->r_con.rcv[i].pr, tag).req(); // tag == SeqNum in c++ ver
		pthread_mutex_unlock(&(rg->lMap[f]->r_con.rcvLock));
		pthread_mutex_unlock(&(rg->rMap[f]->r_con.rcvLock));
	    }
	}//for num messages in each Fab
    }// for(f<numfabs)

    //handle send requests. We don't need to lock the queues because no one changes the front queue data.
    for(int f=0; f<numfabs; f++)
    {
	for(int i=0; i<rg->sMap[f]->r_con.nsnd; i++)
	{
	    if(rg->sMap[f]->r_con.snd[i].pQueue.queueSize(true) == 0) //then !no message has been issued or all send requests have been fulfilled
		nextsReq = false;
	    else
		nextsReq = true;

	    if(nextsReq)
	    {
		Package *sMetaPackage = rg->sMap[f]->r_con.snd[i].pQueue.getFront(true);
		if(!sMetaPackage->served)
		{
		    Package *sPackage = rg->lMap[f]->r_con.snd[i].pQueue.getFront(true);
		    sMetaPackage->completed = false;
		    sMetaPackage->served = true;
		    sMetaPackage->request = MPI_REQUEST_NULL;
		    int ns = rg->sMap[f]->r_con.snd[i].ns;
		    int nd = rg->sMap[f]->r_con.snd[i].nd;
		    int r_gid = rg->sMap[f]->r_con.snd[i].r_gid;
		    int r_grids = rg->sMap[f]->r_con.snd[i].r_grids;
		    int tag = Perilla::myTagMap[r_gid][nd][ns][rg->sMap[f]->r_con.snd[i].sz];
		    sMetaPackage->request = ParallelDescriptor::Asend(sPackage->databuf,
			    rg->sMap[f]->r_con.snd[i].sz,
			    rg->sMap[f]->r_con.snd[i].pr, tag).req();  
		}
	    }
	} // for(i<nsnd)
    } // for(f<numfabs)

    //!now we test if send and receive requests have been serviced
    for(int f=0; f<numfabs; f++)
    {
	//!receive requests
	for(int i=0; i<rg->rMap[f]->r_con.nrcv; i++)
	{
	    if(rg->rMap[f]->r_con.rcv[i].pQueue.queueSize(true) > 0) //by our policy, we can assume that all messages before rear have completed
	    {
		//we dont need to lock the queue, because other consumers just take front messages. A circular buffer guarantees that the rear of the queue can be safely accessed when other queue data is modified
		Package *rearPackage =  rg->rMap[f]->r_con.rcv[i].pQueue.getRear(true);
		if(rearPackage) 
		    if(!(rearPackage->completed))
		    {
			bool flag = false;
			int ret_flag;
			MPI_Status status;
			ParallelDescriptor::Test(rearPackage->request, ret_flag, status);
			flag = (ret_flag == 0) ? false : true;
			if(flag)
			{
			    pthread_mutex_lock(&(rg->lMap[f]->r_con.rcvLock));
			    rearPackage->completeRequest(true);
			    rg->lMap[f]->r_con.rcv[i].pQueue.getRear()->completeRequest(true);
			    if(rg->rMap[f]->r_con.rcv[i].pQueue.queueSize(true) == 1)
				rg->lMap[f]->r_con.firingRuleCnt++;
			    pthread_mutex_unlock(&(rg->lMap[f]->r_con.rcvLock));
			}
		    }
	    } // if(queueSize > 0)
	} // for(i<nrcv)
    } // for(f<numfabs)

    for(int f=0; f<numfabs; f++)
    {
	//!send requests
	for(int i=0; i<rg->lMap[f]->r_con.nsnd; i++)
	{
	    if(rg->sMap[f]->r_con.snd[i].pQueue.queueSize(true) > 0)
	    {
		Package *frontPackage = rg->sMap[f]->r_con.snd[i].pQueue.getFront(true);
		if(frontPackage->served && !frontPackage->completed) //!latest receive request has NOT been completed
		{
		    bool flag = false;
		    int ret_flag;
		    MPI_Status status;
		    ParallelDescriptor::Test(frontPackage->request, ret_flag, status);
		    flag = (ret_flag == 0) ? false : true;
		    if(flag)
		    {
			pthread_mutex_lock(&(rg->sMap[f]->r_con.sndLock));
			frontPackage = rg->sMap[f]->r_con.snd[i].pQueue.dequeue(true);
			frontPackage->completed = false;
			frontPackage->served = false;
			frontPackage->request = MPI_REQUEST_NULL;
			rg->sMap[f]->r_con.snd[i].recycleQueue.enqueue(frontPackage,true);
			pthread_mutex_unlock(&(rg->sMap[f]->r_con.sndLock));
			pthread_mutex_lock(&(rg->lMap[f]->r_con.sndLock));
			frontPackage = rg->lMap[f]->r_con.snd[i].pQueue.dequeue(true);
			frontPackage->completed = false;
			frontPackage->served = false;
			frontPackage->request = MPI_REQUEST_NULL;
			rg->lMap[f]->r_con.snd[i].recycleQueue.enqueue(frontPackage,true);
			pthread_mutex_unlock(&(rg->lMap[f]->r_con.sndLock));			
		    }
		}
	    } // if(queueSize > 0)
	} // for(i<nsnd)
    }// for(f<numfabs)

}// serviceRemoteRequests

void Perilla::serviceMultipleGraphCommDynamic(std::vector<RegionGraph*> graphArray, bool cpyAcross, int tid)
{
    int tg = WorkerThread::perilla_wid();
    int np = ParallelDescriptor::NProcs();
    int nGraphs = graphArray.size();

    for(int g=0; g<nGraphs; g++)
    {
#ifdef PERILLA_DEBUG
	double time= -getTime();
#endif
	serviceLocalRequests(graphArray[g], tg);
	if(cpyAcross)
	    serviceLocalGridCopyRequests(graphArray,g,tg);
#ifdef PERILLA_DEBUG
	time+= getTime();
	if(tg==0)localScheDelay+= time;
#endif
	if(np > 1)
	{
	    if(tg==0)
	    {
		serviceRemoteRequests(graphArray[g],g,nGraphs);
		if(cpyAcross)
		    serviceRemoteGridCopyRequests(graphArray,g,nGraphs,tg);
	    }
	}
    }
}//serviceMultipleGraphCommDynamic

void Perilla::multifabCopyPull(RegionGraph* destGraph, RegionGraph* srcGraph, MultiFab* mfDst, MultiFab* mfSrc, int f, int dstcomp, int srccomp, int nc, int ng, int ngsrc, bool singleT)
{
    int myProc = ParallelDescriptor::MyProc();

    int ntid = WorkerThread::perilla_wtid();
    int tg = WorkerThread::perilla_wid();
    if(nc<1) cout <<"MULTIFAB_COPY_C: nc must be >= 1"<< endl;
    if(mfDst->nComp() < (dstcomp-1)) cout <<"MULTIFAB_COPY_C: nc too large for dst multifab"<< endl;

    if(true)//if(!(*mfDst == *mfSrc))
    {
	if(ng > mfDst->nGrow()) cout <<"MULTIFAB_COPY_C: ng > 0 not supported in parallel copy"<< endl;
	FabCopyAssoc* cpDst = destGraph->task[f]->cpAsc_dstHead;
	while(cpDst != 0)
	{   
	    if(cpDst->graphPartner == srcGraph)
		break;
	    cpDst = cpDst->next;
	}
	if(cpDst == 0) cout <<"Metadata for across grid copy not found"<< endl;
	if(singleT)
	{
	    pthread_mutex_lock(&(cpDst->l_con.dLock));
	    for(int i=0; i<cpDst->l_con.ndcpy; i++)
	    {
		Package* rcvPackage = cpDst->l_con.dcpy[i].pQueue.getFront(true); // corrected from recycleQ to pQ
		mfDst->m_fabs_v[f]->copyFromMem(cpDst->l_con.dcpy[i].dbx,dstcomp,nc,rcvPackage->databuf);
	    }
	    for(int i=0; i<cpDst->l_con.ndcpy; i++)
		cpDst->l_con.dcpy[i].recycleQueue.enqueue(cpDst->l_con.dcpy[i].pQueue.dequeue()); // corrected from pQ to recycleQ and from recycleQ to pQ
	    cpDst->l_con.firingRuleCnt = cpDst->l_con.firingRuleCnt - cpDst->l_con.ndcpy;
	    pthread_mutex_unlock(&(cpDst->l_con.dLock));
	}
	else
	{
	    if(ntid==0){
		pthread_mutex_lock(&(cpDst->l_con.dLock));
		for(int i=0; i<cpDst->l_con.ndcpy; i++)
		{
		    Package* rcvPackage = cpDst->l_con.dcpy[i].pQueue.getFront(true); // corrected from recycleQ to pQ
		    mfDst->m_fabs_v[f]->copyFromMem(cpDst->l_con.dcpy[i].dbx,dstcomp,nc,rcvPackage->databuf);
		}
		for(int i=0; i<cpDst->l_con.ndcpy; i++)
		    cpDst->l_con.dcpy[i].recycleQueue.enqueue(cpDst->l_con.dcpy[i].pQueue.dequeue()); // corrected from pQ to recycleQ and from recycleQ to pQ
		cpDst->l_con.firingRuleCnt = cpDst->l_con.firingRuleCnt - cpDst->l_con.ndcpy;
		pthread_mutex_unlock(&(cpDst->l_con.dLock));
	    }
	}

	int np = ParallelDescriptor::NProcs();
	if(np == 1)
	    return;

	if(singleT)
	{
	    pthread_mutex_lock(&(cpDst->r_con.rcvLock));
	    for(int i=0; i<cpDst->r_con.nrcv; i++)
	    {
		///*
		Package* rcvPackage = cpDst->r_con.rcv[i].pQueue.dequeue(true);                               // corrected from recycleQ to pQ
		mfDst->m_fabs_v[f]->copyFromMem(cpDst->r_con.rcv[i].dbx,dstcomp,nc,rcvPackage->databuf); 
		rcvPackage->completed = false;
		rcvPackage->served = false;
		rcvPackage->request = MPI_REQUEST_NULL;
		cpDst->r_con.rcv[i].recycleQueue.enqueue(rcvPackage, true);                         // corrected from pQ to recycleQ           
	    }
	    cpDst->r_con.firingRuleCnt = cpDst->r_con.firingRuleCnt - cpDst->r_con.nrcv;

	    cpDst->r_con.remotePullDone = true;
	    ///*
	    for(int i=0; i<cpDst->r_con.nrcv; i++)
		if(cpDst->r_con.rcv[i].pQueue.queueSize(true) >= 1)
		    if(cpDst->r_con.rcv[i].pQueue.getFront(true)->checkRequest())
			cpDst->r_con.firingRuleCnt++;
	    //*/
	    pthread_mutex_unlock(&(cpDst->r_con.rcvLock));
	}
	else
	{
	    if(ntid==0)
	    {
		pthread_mutex_lock(&(cpDst->r_con.rcvLock));
		for(int i=0; i<cpDst->r_con.nrcv; i++)
		{
		    Package* rcvPackage = cpDst->r_con.rcv[i].pQueue.dequeue(true);                               // corrected from recycleQ to pQ
		    mfDst->m_fabs_v[f]->copyFromMem(cpDst->r_con.rcv[i].dbx,dstcomp,nc,rcvPackage->databuf);
		    rcvPackage->completed = false;
		    rcvPackage->served = false;
		    rcvPackage->request = MPI_REQUEST_NULL;
		    cpDst->r_con.rcv[i].recycleQueue.enqueue(rcvPackage, true);                         // corrected from pQ to recycleQ       
		}
		cpDst->r_con.firingRuleCnt = cpDst->r_con.firingRuleCnt - cpDst->r_con.nrcv;

		cpDst->r_con.remotePullDone = true;
		for(int i=0; i<cpDst->r_con.nrcv; i++)
		    if(cpDst->r_con.rcv[i].pQueue.queueSize() >= 1)
			if(cpDst->r_con.rcv[i].pQueue.getFront()->checkRequest())
			    cpDst->r_con.firingRuleCnt++;
		pthread_mutex_unlock(&(cpDst->r_con.rcvLock));
	    }
	}
    } // if(!(*mfDst == *mfSrc))
} // multifabCopyPull


void Perilla::fillBoundaryPull(RegionGraph* graph, MultiFab* mf, int f, bool singleT)
{

    int nComp = mf->nComp();
    int tg= WorkerThread::perilla_wid();
    int ntid = WorkerThread::perilla_wtid();

    if(ntid==0)
	pthread_mutex_lock(&(graph->lMap[f]->l_con.dLock));
    if(!singleT)
	graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads    

    for(int i=0; i<graph->lMap[f]->l_con.ndcpy; i++)
	if( (i%(perilla::NUM_THREADS_PER_TEAM-1)) == ntid)
	{
	    Package *dPackage = graph->lMap[f]->l_con.dcpy[i].pQueue.getFront(true);
	    mf->m_fabs_v[f]->copyFromMem(graph->lMap[f]->l_con.dcpy[i].dbx,0,nComp,dPackage->databuf);
	}

    if(!singleT)
	graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads
    if(ntid==0)
    {
	for(int i=0; i<graph->lMap[f]->l_con.ndcpy; i++)
	    graph->lMap[f]->l_con.dcpy[i].recycleQueue.enqueue( graph->lMap[f]->l_con.dcpy[i].pQueue.dequeue(true),true );

	graph->lMap[f]->l_con.firingRuleCnt = graph->lMap[f]->l_con.firingRuleCnt - graph->lMap[f]->l_con.ndcpy;

	graph->lMap[f]->l_con.scpyCnt = 0;
	for(int i=0; i<graph->lMap[f]->l_con.ndcpy; i++)
	    if(graph->lMap[f]->l_con.dcpy[i].pQueue.queueSize(true) >= 1)
		graph->lMap[f]->l_con.firingRuleCnt++;
	pthread_mutex_unlock(&(graph->lMap[f]->l_con.dLock));
    }
    if(!singleT)
	graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads

    int np = ParallelDescriptor::NProcs();
    if (np==1) return;

    if(ntid==0)
    {
	pthread_mutex_lock(&(graph->rMap[f]->r_con.rcvLock));
	pthread_mutex_lock(&(graph->lMap[f]->r_con.rcvLock));
    }
    if(!singleT)
	graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads

    for(int i=0; i<graph->lMap[f]->r_con.nrcv; i++)
	if( (i%(perilla::NUM_THREADS_PER_TEAM-1)) == ntid)
	{
	    Package *rcvMetaPackage = graph->rMap[f]->r_con.rcv[i].pQueue.dequeue(true);
	    rcvMetaPackage->completed = false;
	    rcvMetaPackage->served = false;
	    rcvMetaPackage->request = 0;
	    graph->rMap[f]->r_con.rcv[i].recycleQueue.enqueue(rcvMetaPackage,true);
	    Package *rcvPackage = graph->lMap[f]->r_con.rcv[i].pQueue.dequeue(true);
	    mf->m_fabs_v[f]->copyFromMem(graph->lMap[f]->r_con.rcv[i].dbx,0,nComp,rcvPackage->databuf);
	    rcvPackage->completed = false;
	    graph->lMap[f]->r_con.rcv[i].recycleQueue.enqueue(rcvPackage,true);
	}
    if(!singleT)
	graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads
    if(ntid==0)
    {
	graph->lMap[f]->r_con.firingRuleCnt = graph->lMap[f]->r_con.firingRuleCnt - graph->lMap[f]->r_con.nrcv;
	for(int i=0; i<graph->lMap[f]->r_con.nrcv; i++)
	    if(graph->lMap[f]->r_con.rcv[i].pQueue.queueSize(true) >= 1)
		if(graph->lMap[f]->r_con.rcv[i].pQueue.getFront(true)->checkRequest())
		    graph->lMap[f]->r_con.firingRuleCnt++;
	pthread_mutex_unlock(&(graph->lMap[f]->r_con.rcvLock));
	pthread_mutex_unlock(&(graph->rMap[f]->r_con.rcvLock));
    }

} // fillBoundaryPull



void Perilla::fillBoundaryPull_1Team(RegionGraph* graph, amrex::MultiFab& mf, int f)
{
    exit(0);
} // fillBoundaryPull


void Perilla::multifabCopyPushAsync(RegionGraph* destGraph, RegionGraph* srcGraph, MultiFab* mfDst, MultiFab* mfSrc, int f, int dstcomp, int srccomp, int nc, int ng, int ngsrc, bool singleT)
{
    int ntid = WorkerThread::perilla_wtid();
    int tg = WorkerThread::perilla_wid();
    int myProc = ParallelDescriptor::MyProc();
    // MultiFab* mfDst = destGraph->assocMF;
    // MultiFab* mfSrc = srcGraph->assocMF;
    if(nc<1) cout <<"MULTIFAB_COPY_C: nc must be >= 1"<< endl;
    if(mfDst->nComp() < (dstcomp-1)) cout <<"MULTIFAB_COPY_C: nc too large for dst multifab"<< endl;
    if(mfSrc->nComp() < (srccomp-1)) cout <<"MULTIFAB_COPY_C: nc too large for src multifab"<< endl;

    if(true)//if(!(*mfDst == *mfSrc))
    {
	if(ng > mfDst->nGrow()) cout <<"MULTIFAB_COPY_C: ng > 0 not supported in parallel copy"<< endl;
	if(ngsrc > mfSrc->nGrow()) cout <<"MULTIFAB_COPY_C: ngsrc > msrc%ng"<< endl;
	FabCopyAssoc* cpSrc = srcGraph->task[f]->cpAsc_srcHead;

	while(cpSrc != 0)
	{
	    if(cpSrc->graphPartner == destGraph)
		break;
	    cpSrc = cpSrc->next;
	}
	if(cpSrc == 0) cout <<"Metadata for across grid copy not found"<< endl;

	if(singleT)
	{
	    pthread_mutex_lock(&(cpSrc->l_con.sLock));
	    for(int i=0; i<cpSrc->l_con.nscpy; i++)
	    {
		Package* sndPackage = cpSrc->l_con.scpy[i].recycleQueue.getFront(true);
		mfSrc->m_fabs_v[f]->copyToMem(cpSrc->l_con.scpy[i].sbx,srccomp,nc,sndPackage->databuf);
	    }
	    for(int i=0;i<cpSrc->l_con.nscpy; i++)
		cpSrc->l_con.scpy[i].pQueue.enqueue(cpSrc->l_con.scpy[i].recycleQueue.dequeue(true),true);
	    pthread_mutex_unlock(&(cpSrc->l_con.sLock));
	}
	else
	{
	    if(ntid == 0)
		pthread_mutex_lock(&(cpSrc->l_con.sLock));
	    srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
	    for(int i=0; i<cpSrc->l_con.nscpy; i++)
		if((i%(perilla::NUM_THREADS_PER_TEAM-1)) == ntid)
		{
		    Package* sndPackage = cpSrc->l_con.scpy[i].recycleQueue.getFront(true);
		    mfSrc->m_fabs_v[f]->copyToMem(cpSrc->l_con.scpy[i].sbx,srccomp,nc,sndPackage->databuf);
		}
	    srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
	    if(ntid==0)
	    {
		for(int i=0;i<cpSrc->l_con.nscpy; i++)
		    cpSrc->l_con.scpy[i].pQueue.enqueue(cpSrc->l_con.scpy[i].recycleQueue.dequeue(true),true);
		pthread_mutex_unlock(&(cpSrc->l_con.sLock));
	    }
	    srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
	}

	int np = ParallelDescriptor::NProcs();
	if(np == 1)
	    return;

	if(singleT)
	{
	    pthread_mutex_lock(&(cpSrc->r_con.sndLock));
	    for(int i=0; i<cpSrc->r_con.nsnd; i++)
	    {

		Package* sndPackage = cpSrc->r_con.snd[i].recycleQueue.dequeue(true);
		mfSrc->m_fabs_v[f]->copyToMem(cpSrc->r_con.snd[i].sbx,srccomp,nc,sndPackage->databuf);
		cpSrc->r_con.snd[i].pQueue.enqueue(sndPackage,true);
	    }

	    pthread_mutex_unlock(&(cpSrc->r_con.sndLock));

	    cpSrc->r_con.remotePushReady = true;
	    ///*
	    pthread_mutex_lock(&(srcGraph->sCopyMapHead->map[f]->r_con.sndLock));
	    for(int i=0; i<cpSrc->r_con.nsnd; i++)
		srcGraph->sCopyMapHead->map[f]->r_con.snd[i].pQueue.enqueue(srcGraph->sCopyMapHead->map[f]->r_con.snd[i].recycleQueue.dequeue(true),true);
	    pthread_mutex_unlock(&(srcGraph->sCopyMapHead->map[f]->r_con.sndLock));
	}
	else
	{
	    if(ntid == 0)
		pthread_mutex_lock(&(cpSrc->r_con.sndLock));
	    srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);

	    for(int i=0; i<cpSrc->r_con.nsnd; i++)
		if((i%(perilla::NUM_THREADS_PER_TEAM-1)) == ntid)
		{
		    Package* sndPackage = cpSrc->r_con.snd[i].recycleQueue.dequeue(true);
		    mfSrc->m_fabs_v[f]->copyToMem(cpSrc->r_con.snd[i].sbx,srccomp,nc,sndPackage->databuf);
		    cpSrc->r_con.snd[i].pQueue.enqueue(sndPackage,true);
		}

	    srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
	    if(ntid==0)
	    {
		pthread_mutex_unlock(&(cpSrc->r_con.sndLock));

		cpSrc->r_con.remotePushReady = true;
		///*
		pthread_mutex_lock(&(srcGraph->sCopyMapHead->map[f]->r_con.sndLock));
		for(int i=0; i<cpSrc->r_con.nsnd; i++)
		    srcGraph->sCopyMapHead->map[f]->r_con.snd[i].pQueue.enqueue(srcGraph->sCopyMapHead->map[f]->r_con.snd[i].recycleQueue.dequeue(true),true);
		pthread_mutex_unlock(&(srcGraph->sCopyMapHead->map[f]->r_con.sndLock));
		//*/
	    }
	    srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
	}
    } // if(!(*mfDst == *mfSrc))                                                                                                              
} // multifabCopyPushAsync


void Perilla::fillBoundaryPush(RegionGraph* graph, MultiFab* mf, int f)
{

    int nComp = mf->nComp();
    int tg= WorkerThread::perilla_wid();
    int ntid = WorkerThread::perilla_wtid();

    if(ntid == 0)
	pthread_mutex_lock(&(graph->lMap[f]->l_con.sLock));
    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads

    for(int i=0; i<graph->lMap[f]->l_con.nscpy; i++)
	if( (i%(perilla::NUM_THREADS_PER_TEAM-1)) == ntid)
	{
	    Package *sPackage = graph->lMap[f]->l_con.scpy[i].recycleQueue.getFront(true);
	    mf->m_fabs_v[f]->copyToMem(graph->lMap[f]->l_con.scpy[i].sbx,0,nComp,sPackage->databuf);
	}
    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads

    if(ntid==0)
    {
	for(int i=0; i<graph->lMap[f]->l_con.nscpy; i++)
	{
	    graph->lMap[f]->l_con.scpy[i].pQueue.enqueue( graph->lMap[f]->l_con.scpy[i].recycleQueue.dequeue(true),true );
	}
	pthread_mutex_unlock(&(graph->lMap[f]->l_con.sLock));
    }
    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads

    int np = ParallelDescriptor::NProcs();
    if (np==1) return;

    if(ntid==0)
	pthread_mutex_lock(&(graph->lMap[f]->r_con.sndLock));
    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads

    for(int i=0; i<graph->lMap[f]->r_con.nsnd; i++)
	if((i%(perilla::NUM_THREADS_PER_TEAM-1))==ntid)
	{
	    Package *sndPackage = graph->lMap[f]->r_con.snd[i].recycleQueue.dequeue(true);
	    mf->m_fabs_v[f]->copyToMem(graph->lMap[f]->r_con.snd[i].sbx,0,nComp,sndPackage->databuf);
	    graph->lMap[f]->r_con.snd[i].pQueue.enqueue( sndPackage,true );
	}
    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads
    if(ntid==0)
    {
	pthread_mutex_unlock(&(graph->lMap[f]->r_con.sndLock));
	pthread_mutex_lock(&(graph->sMap[f]->r_con.sndLock));
	for(int i=0; i<graph->lMap[f]->r_con.nsnd; i++)
	    graph->sMap[f]->r_con.snd[i].pQueue.enqueue( graph->sMap[f]->r_con.snd[i].recycleQueue.dequeue(true),true );
	pthread_mutex_unlock(&(graph->sMap[f]->r_con.sndLock));
    }
} // fillBoundaryPush



void Perilla::fillBoundaryPush(amrex::RGIter& rgi, amrex::MultiFab& mf)
{
    if(rgi.currentItr != rgi.totalItr)
	return;

    int f = rgi.currentRegion;
    fillBoundaryPush(rgi.itrGraph, &mf, f);
}

void Perilla::fillBoundaryPush(amrex::RGIter& rgi, RegionGraph* rg, amrex::MultiFab& mf)
{
    if(rgi.currentItr != rgi.totalItr)
	return;

    int f = rgi.currentRegion;
    fillBoundaryPush(rg, &mf, f);
}


void Perilla::multifabCopyPush(RegionGraph* destGraph, RegionGraph* srcGraph, amrex::MultiFab* mfDst, amrex::MultiFab* mfSrc, int f, int dstcomp, int srccomp, int nc, int ng, int ngsrc, bool singleT)
{
    if(nc<1) cout <<"MULTIFAB_COPY_C: nc must be >= 1"<< endl;
    if(mfDst->nComp() < (dstcomp-1)) cout <<"MULTIFAB_COPY_C: nc too large for dst multifab"<< endl;
    if(mfSrc->nComp() < (srccomp-1)) cout <<"MULTIFAB_COPY_C: nc too large for src multifab"<< endl;

    multifabCopyPush_1Team(destGraph,srcGraph,mfDst,mfSrc,f,dstcomp,srccomp,nc,ng,ngsrc,singleT);
}


void Perilla::multifabCopyPush_1Team(RegionGraph* destGraph, RegionGraph* srcGraph, amrex::MultiFab* mfDst, amrex::MultiFab* mfSrc, int f, int dstcomp, int srccomp, int nc, int ng, int ngsrc, bool singleT)
{
    int ntid = perilla::wtid();// - perilla::NUM_COMM_THREADS;
    int tg = perilla::wid();
    int myProc = amrex::ParallelDescriptor::MyProc();

    if(true)//if(!(*mfDst == *mfSrc))
    {
	if(ng > mfDst->nGrow()) cout <<"MULTIFAB_COPY_C: ng > 0 not supported in parallel copy"<< endl;
	if(ngsrc > mfSrc->nGrow()) cout <<"MULTIFAB_COPY_C: ngsrc > msrc%ng"<< endl;
	FabCopyAssoc* cpSrc = srcGraph->task[f]->cpAsc_srcHead;

	while(cpSrc != 0)
	{
	    if(cpSrc->graphPartner == destGraph)
		break;
	    cpSrc = cpSrc->next;
	}
	if(cpSrc == 0) cout <<"Metadata for across grid copy not found"<< endl;

	if(singleT)
	{
	    pthread_mutex_lock(&(cpSrc->l_con.sLock));
	    for(int i=0; i<cpSrc->l_con.nscpy; i++)
	    {
		Package* sndPackage = cpSrc->l_con.scpy[i].recycleQueue.getFront(true);
		mfSrc->m_fabs_v[f]->copyToMem(cpSrc->l_con.scpy[i].sbx,srccomp,nc,sndPackage->databuf);
	    }
	    for(int i=0;i<cpSrc->l_con.nscpy; i++)
		cpSrc->l_con.scpy[i].pQueue.enqueue(cpSrc->l_con.scpy[i].recycleQueue.dequeue(true));
	    pthread_mutex_unlock(&(cpSrc->l_con.sLock));
	}
	else
	{
	    if(ntid == 0)
	    {
		pthread_mutex_lock(&(cpSrc->l_con.sLock));
		for(int i=0; i<cpSrc->l_con.nscpy; i++)
		{
		    Package* sndPackage = cpSrc->l_con.scpy[i].recycleQueue.getFront(true);
		    mfSrc->m_fabs_v[f]->copyToMem(cpSrc->l_con.scpy[i].sbx,srccomp,nc,sndPackage->databuf);
		}

		for(int i=0;i<cpSrc->l_con.nscpy; i++)
		    cpSrc->l_con.scpy[i].pQueue.enqueue(cpSrc->l_con.scpy[i].recycleQueue.dequeue(true));
		pthread_mutex_unlock(&(cpSrc->l_con.sLock));
	    }
	}

	int np = amrex::ParallelDescriptor::NProcs();
	if(np == 1)
	    return;
	if(singleT)
	{
	    pthread_mutex_lock(&(cpSrc->r_con.sndLock));
	    for(int i=0; i<cpSrc->r_con.nsnd; i++)
	    {
		Package* sndPackage = cpSrc->r_con.snd[i].recycleQueue.dequeue(true);
		mfSrc->m_fabs_v[f]->copyToMem(cpSrc->r_con.snd[i].sbx,srccomp,nc,sndPackage->databuf);
		sndPackage->served = false;
		sndPackage->completed = false;
		cpSrc->r_con.snd[i].pQueue.enqueue(sndPackage, true);
	    }
	    cpSrc->r_con.remotePushReady = true;
	    pthread_mutex_unlock(&(cpSrc->r_con.sndLock));
	}
	else
	{
	    if(ntid == 0)
	    {
		pthread_mutex_lock(&(cpSrc->r_con.sndLock));
		for(int i=0; i<cpSrc->r_con.nsnd; i++)
		{
		    Package* sndPackage = cpSrc->r_con.snd[i].recycleQueue.dequeue(true);
		    mfSrc->m_fabs_v[f]->copyToMem(cpSrc->r_con.snd[i].sbx,srccomp,nc,sndPackage->databuf);
		    sndPackage->served = false;
		    sndPackage->completed = false;
		    cpSrc->r_con.snd[i].pQueue.enqueue(sndPackage, true);
		}
		cpSrc->r_con.remotePushReady = true;
		pthread_mutex_unlock(&(cpSrc->r_con.sndLock));
	    }
	}
    } // if(!(*mfDst == *mfSrc))                                                                                                                    
} // multifabCopyPush


void Perilla::serviceLocalGridCopyRequests(std::vector<RegionGraph*> graphArray, int g, int tg)
{
    int nfabs = graphArray[g]->numTasks;

    for(int f=0; f<nfabs; f++)
    {
	if(WorkerThread::isMyRegion(tg,f)) //tg == fg
	{
	    FabCopyAssoc* cpSrc = graphArray[g]->task[f]->cpAsc_srcHead;
	    while(cpSrc != 0)
	    {
		bool anyReq=false;
		for(int i=0; i<cpSrc->l_con.nscpy; i++)
		    if(cpSrc->l_con.scpy[i].pQueue.queueSize(true)>0){
			anyReq=true;
			break;
		    }
		if(anyReq)
		{
		    pthread_mutex_lock(&(cpSrc->l_con.sLock));
		    for(int i=0; i<cpSrc->l_con.nscpy; i++)
		    {
			if(cpSrc->l_con.scpy[i].pQueue.queueSize(true)>0)
			{
			    FabCopyAssoc* cpDst = cpSrc->graphPartner->task[cpSrc->l_con.scpy[i].nd]->cpAsc_dstHead;
			    while(cpDst != 0)
			    {
				if(cpDst->graphPartner == graphArray[g])
				    break;
				cpDst = cpDst->next;
			    }
			    Package* sPackage = cpSrc->l_con.scpy[i].pQueue.dequeue(true);
			    pthread_mutex_lock(&(cpDst->l_con.dLock));
			    int dPartner = cpSrc->l_con.scpy[i].dPartner;
			    Package* dPackage = cpDst->l_con.dcpy[dPartner].recycleQueue.dequeue(true);
			    std::memcpy(dPackage->databuf, sPackage->databuf, dPackage->bufSize * sizeof(double));
			    cpDst->l_con.dcpy[dPartner].pQueue.enqueue(dPackage,true);
			    if(cpDst->l_con.dcpy[dPartner].pQueue.queueSize(true) == 1)
				cpDst->l_con.firingRuleCnt++;
			    pthread_mutex_unlock(&(cpDst->l_con.dLock));
			    cpSrc->l_con.scpy[i].recycleQueue.enqueue(sPackage,true);
			}
		    } // for
		    pthread_mutex_unlock(&(cpSrc->l_con.sLock));
		}//anyReq
		cpSrc = cpSrc->next;
	    } // while(cpSrc != 0)
	} // if(tg==fg)
    } // for(f<nfabs)
} // serviceLocalGridCopyRequests

void Perilla::serviceRemoteGridCopyRequests(std::vector<RegionGraph*> graphArray, int g, int nGraphs, int tg)
{
    bool nextsReq, nextrReq;
    int np = ParallelDescriptor::NProcs();
    int myProc = ParallelDescriptor::MyProc();
    int numfabs = graphArray[g]->numTasks;
    int graphID = graphArray[g]->graphID;

#ifdef PERILLA_DEBUG
    double time= -getTime();
#endif
    for(int f=0; f<numfabs; f++)
    {	
	FabCopyAssoc* cpDst = graphArray[g]->task[f]->cpAsc_dstHead;
	while(cpDst != 0)
	{
	    for(int i=0; i<cpDst->r_con.nrcv; i++)
	    {
		if(cpDst->r_con.rcv[i].pQueue.queueSize(true)==0)
		{
		    nextrReq = true;
		}
		else
		{
		    Package *rearPackage = cpDst->r_con.rcv[i].pQueue.getRear(true);
		    if(rearPackage) 
			if(rearPackage->completed && cpDst->r_con.rcv[i].pQueue.queueSize(true) == 1) //!latest receive request has been completed
			{
			    nextrReq = true;
			}
			else //!expected message is still on the way
			    nextrReq = false;
		}
		if(nextrReq) //!take a message from recycle pool and post a receive
		{
		    pthread_mutex_lock(&(cpDst->r_con.rcvLock));
		    int ns = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].ns;
		    int nd = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].nd;
		    int lnd = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].lnd;
		    int r_grids = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].r_grids;
		    Package *rPackage = cpDst->r_con.rcv[i].recycleQueue.dequeue(true);
		    int tag = tagMap[graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pr][g][nd][ns][graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].sz];
		    rPackage->request = MPI_REQUEST_NULL;
		    rPackage->completed=false;
		    cpDst->r_con.rcv[i].pQueue.enqueue(rPackage, true);   //!this is not done yet
		    rPackage->request = ParallelDescriptor::Arecv(rPackage->databuf,
			    graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].sz,
			    graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pr, tag).req(); // tag == SeqNum in c++ ver
		    pthread_mutex_unlock(&(cpDst->r_con.rcvLock));
		}						
	    } // for (i<i<cpDst->r_con.nrcv)
	    cpDst = cpDst->next;
	} // while(cpDst != 0)	
    } // for(f<nfabs)
#ifdef PERILLA_DEBUG
    time+= getTime();
    irecvDelay+= time;

    time= -getTime();
#endif
    for(int f=0; f<numfabs; f++)
    {	
	FabCopyAssoc* cpSrc = graphArray[g]->task[f]->cpAsc_srcHead;
	while(cpSrc != 0)
	{
	    for(int i=0; i<cpSrc->r_con.nsnd; i++)
	    {
		if(cpSrc->r_con.snd[i].pQueue.queueSize(true) == 0)
		    nextsReq = false;
		else
		    nextsReq = true;

		if(nextsReq) 
		{
		    //there is no need to lock the queue because we only touch the front to initialize the send
		    //During this time, workers can produce more messages into the queue, but a circular queue ensures that the front of the queue will not be modified
		    Package *sPackage = cpSrc->r_con.snd[i].pQueue.getFront(true);
		    if(!sPackage->served)
		    {  		    
			sPackage->completed = false;
			sPackage->served = true;
			sPackage->request = MPI_REQUEST_NULL;
			int ns = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].ns;
			int nd = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].nd;
			int r_gid = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].r_gid;
			int r_grids = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].r_grids;
			int tag = Perilla::myTagMap[r_gid][nd][ns][graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].sz];
			sPackage->request = ParallelDescriptor::Asend(sPackage->databuf,
				graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].sz,
				graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].pr, tag).req();  // tag == SeqNum in c++ ver
		    }
		}		
	    } // for (i<i<cpSrc->r_con.nsnd)	    
	    cpSrc = cpSrc->next;
	} // while(cpSrc != 0)	
    } // for(f<nfabs)
#ifdef PERILLA_DEBUG
    time+= getTime();
    isendDelay+= time;

    time= -getTime();
#endif
    for(int f=0; f<numfabs; f++)
    {	
	FabCopyAssoc* cpDst = graphArray[g]->task[f]->cpAsc_dstHead;
	while(cpDst != 0)
	{
	    for(int i=0; i<cpDst->r_con.nrcv; i++)
	    {
		if(cpDst->r_con.rcv[i].pQueue.queueSize(true) > 0) 
		{		    
		    Package *rearPackage =  cpDst->r_con.rcv[i].pQueue.getRear(true);
		    //Note: all messages before rear have completed
		    if(rearPackage)
			if(!rearPackage->completed)
			{
			    bool flag = false;
			    int ret_flag=0;
			    MPI_Status status;
			    ParallelDescriptor::Test(rearPackage->request, ret_flag, status);

			    flag = (ret_flag == 0) ? false : true;
			    if(flag)
			    {
				pthread_mutex_lock(&(cpDst->r_con.rcvLock));
				rearPackage->completeRequest(true);				
				if(cpDst->r_con.rcv[i].pQueue.queueSize(true) == 1)
				{
				    cpDst->r_con.firingRuleCnt++;
				}
				pthread_mutex_unlock(&(cpDst->r_con.rcvLock));
			    }
			}		   		    
		} // if(pQueue.queueSize(true) > 0)		    
	    } // for (i<i<cpDst->r_con.nrcv)
	    cpDst = cpDst->next;
	} // while(cpDst != 0)	
    } // for(f<nfabs)
#ifdef PERILLA_DEBUG
    time+= getTime();
    irecvTestDelay+= time;

    time= -getTime();
#endif
    for(int f=0; f<numfabs; f++)
    {	
	FabCopyAssoc* cpSrc = graphArray[g]->task[f]->cpAsc_srcHead;
	while(cpSrc != 0)
	{
	    for(int i=0; i<cpSrc->r_con.nsnd; i++)
	    {		
		if(cpSrc->r_con.snd[i].pQueue.queueSize(true) >0)
		{
		    Package *frontPackage = cpSrc->r_con.snd[i].pQueue.getFront(true);
		    if(frontPackage->served)
		    {
			bool flag = false;
			int ret_flag;
			MPI_Status status;
			ParallelDescriptor::Test(frontPackage->request, ret_flag, status);
			flag = (ret_flag == 0) ? false : true;
			if(flag)
			{
			    //we have to lock the queue before removing the front 
			    pthread_mutex_lock(&(cpSrc->r_con.sndLock));
			    frontPackage = cpSrc->r_con.snd[i].pQueue.dequeue(true);
			    frontPackage->completed = false;
			    frontPackage->served = false;
			    frontPackage->request = MPI_REQUEST_NULL;
			    cpSrc->r_con.snd[i].recycleQueue.enqueue(frontPackage, true);
			    pthread_mutex_unlock(&(cpSrc->r_con.sndLock));			
			}
		    }
		} // if(queueSize > 0)				
	    } // for (i<i<cpSrc->r_con.nsnd)	    
	    cpSrc = cpSrc->next;
	} // while(cpSrc != 0)	
    } // for(f<nfabs)
#ifdef PERILLA_DEBUG
    time+= getTime();
    isendTestDelay+= time;
#endif
} // serviceRemoteGridCopyRequests

