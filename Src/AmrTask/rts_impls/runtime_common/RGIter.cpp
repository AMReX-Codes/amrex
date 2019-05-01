#include <AMReX_Box.H>
#include <RGIter.H>
#include <WorkerThread.H>
#include <Perilla.H>
#include <cmath>

#include <AMReX_AmrLevel.H>
#include <PerillaConfig.H>
using namespace perilla;
#include <PerillaRts.H>

#ifdef USE_PERILLA_ON_DEMAND
    pthread_mutex_t teamFinishLock=PTHREAD_MUTEX_INITIALIZER;
#endif

namespace amrex{

    RGIter::RGIter(RegionGraph* rg
#ifdef USE_PERILLA_ON_DEMAND
	,std::vector<RegionGraph*> graphArray
#endif    
	, bool enableAllTasks
	):
	itrGraph(rg),
	implicit(false),
	ppteams(true),
	//typ(rg->typ),
	haveDepGraph(false),
	depGraph(NULL),
	getFireableTime(0.)
    {
	tid = perilla::tid();
	tg = perilla::wid();
	ntid = perilla::wtid();

#ifdef USE_PERILLA_ON_DEMAND
        if(perilla::isCommunicationThread())
        {
            while(true){
                Perilla::serviceMultipleGraphCommDynamic(graphArray,true,perilla::tid());
                if( Perilla::numTeamsFinished == perilla::NUM_THREAD_TEAMS)
                {
		    graphArray.clear();
		    Perilla::numTeamsFinished=0;
                    break;
                }
            }
        }else{
#endif
	itrGraph->worker[tg]->l_barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
	if(perilla::isMasterWorkerThread())
	    itrGraph->Reset();
	itrGraph->worker[tg]->l_barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
	if(enableAllTasks)
	    itrGraph->enableAllRegions();
	itrGraph->worker[tg]->l_barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
	init();
#ifdef USE_PERILLA_ON_DEMAND
}
#endif
    }

    RGIter::RGIter(RegionGraph* rg
#ifdef USE_PERILLA_ON_DEMAND
        ,std::vector<RegionGraph*> graphArray
#endif
	, RegionGraph* drg, bool isDep
    ):
	itrGraph(rg),
	implicit(false),
	ppteams(true),
	//typ(rg->typ),
	haveDepGraph(isDep),
	depGraph(drg),
	getFireableTime(0.)
    {
	tid = perilla::tid();
	tg = perilla::wid();
	ntid = perilla::wtid();

#ifdef USE_PERILLA_ON_DEMAND
        if(perilla::isCommunicationThread())
        {
            Perilla::flattenGraphHierarchy(m_level_afpi[iteration-1]->m_amrlevel.parent->graphArray, graphArray);
            while(true){
                Perilla::serviceMultipleGraphCommDynamic(graphArray,true,perilla::tid());
                if( Perilla::numTeamsFinished == perilla::NUM_THREAD_TEAMS)
                {
		    graphArray.clear();
		    Perilla::numTeamsFinished=0;
                    break;
                }
            }
        }else{
#endif
	itrGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
	if(perilla::isMasterWorkerThread()) itrGraph->Reset();
	itrGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
	init();
#ifdef USE_PERILLA_ON_DEMAND
}   
#endif
    }

    RGIter::RGIter(amrex::AsyncFillPatchIterator* afpi, bool enableAllTasks):
	itrGraph(afpi->destGraph),
	implicit(false),
	ppteams(true),
	//typ(afpi->destGraph->typ),
	haveDepGraph(false),
	depGraph(NULL),
	getFireableTime(0.)
    {
	tid = perilla::tid();
	tg = perilla::wid();
	ntid = perilla::wtid();
#ifdef USE_PERILLA_ON_DEMAND
        if(perilla::isCommunicationThread())
        {
            std::vector<RegionGraph*> flattenedGraphArray;
            Perilla::flattenGraphHierarchy(m_level_afpi[iteration-1]->m_amrlevel.parent->graphArray, flattenedGraphArray);
            while(true){
                Perilla::serviceMultipleGraphCommDynamic(flattenedGraphArray,true,perilla::tid());
                if( Perilla::numTeamsFinished == perilla::NUM_THREAD_TEAMS)
                {
		    flattenedGraphArray.clear();
		    Perilla::numTeamsFinished=0;
                    break;
                }
            }
        }else{
#endif  
	itrGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
	if(perilla::isMasterWorkerThread())
	    afpi->Reset();
	itrGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
	if(enableAllTasks)
	    itrGraph->enableAllRegions();
	itrGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);

	init();
#ifdef USE_PERILLA_ON_DEMAND
}   
#endif

    }

#ifndef USE_PERILLA_ON_DEMAND
    RGIter::RGIter(Vector<amrex::AsyncFillPatchIterator*> afpi, Vector<amrex::AsyncFillPatchIterator*> upper_afpi, 
	    amrex::MultiFab& dest, int  bG, double tm, int  ind, int  sc, int nc, int itr):
	itrGraph(afpi[itr-1]->destGraph),
	m_level_afpi(afpi),
	m_upper_level_afpi(upper_afpi),
	boxGrow(bG), 
	time(tm), 
	index(ind), 
	scomp(sc), 
	ncomp(nc), 
	iteration(itr),
	implicit(true),
	ppteams(true),
	//typ(afpi[itr-1]->destGraph->typ),
	haveDepGraph(false),
	depGraph(NULL),
	getFireableTime(0.)
    {
	int myProc = amrex::ParallelDescriptor::MyProc();
	bool push = true;

	tid = perilla::tid();
	tg = perilla::wid();
	ntid = perilla::wtid();
        AsyncFillPatchIterator::initialSend(afpi, upper_afpi, bG, tm, ind, 0, nc, itr);

	itrGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
	if(perilla::isMasterWorkerThread())
	    m_level_afpi[iteration-1]->Reset();
	itrGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);

	if(ntid == perilla::NUM_THREADS_PER_TEAM-2)
	{
	    int f;
	    int level = m_level_afpi[iteration-1]->m_amrlevel.level;
	    double dt = m_level_afpi[iteration-1]->m_amrlevel.parent->dtLevel(level);
	    this->currentItr = 1;
	    this->totalItr = 1;

	    //////////////////////////////////////Push Pull Thread Start/////////////////////////
	    while(m_level_afpi[iteration-1]->destGraph->worker[tg]->completedRegionQueue->queueSize(true) != m_level_afpi[iteration-1]->destGraph->worker[tg]->totalTasks ||
		    m_level_afpi[iteration-1]->destGraph->worker[tg]->computedTasks != m_level_afpi[iteration-1]->destGraph->worker[tg]->totalTasks)
	    {
		f = m_level_afpi[iteration-1]->destGraph->getFireableRegion(tg);
		if(f != -1)
		{
		    m_level_afpi[iteration-1]->Receive(this,dest,boxGrow,time,index,scomp,ncomp,f,true);
		    m_level_afpi[iteration-1]->destGraph->setFireableRegion(f);
		    if(m_level_afpi[iteration-1]->destGraph->worker[tg]->unfireableRegionQueue->queueSize(true) !=0 && 
			    m_level_afpi[iteration-1]->destGraph->worker[tg]->fireableRegionQueue->queueSize(true) < 2)
			continue;
		}

		if(m_level_afpi[iteration-1]->destGraph->worker[tg]->computedRegionQueue->queueSize() != 0)
		{
		    f = m_level_afpi[iteration-1]->destGraph->worker[tg]->computedRegionQueue->removeRegion();

		    if(push & level == m_level_afpi[iteration-1]->m_amrlevel.parent->finestLevel() && iteration < m_level_afpi[iteration-1]->m_amrlevel.parent->nCycle(level))
			m_level_afpi[iteration]->SendIntraLevel(*(this),boxGrow,time+dt,index,scomp,ncomp,iteration,f,true);

		    if(push & level < m_level_afpi[iteration-1]->m_amrlevel.parent->finestLevel())
		    {
			for(int i=0; i < m_level_afpi[iteration-1]->m_amrlevel.parent->nCycle(level+1); i++)
			{
			    m_upper_level_afpi[i]->SendInterLevel(this,boxGrow,time+(i*m_level_afpi[iteration-1]->m_amrlevel.parent->dtLevel(level+1)),index,scomp,ncomp,i+1,f,true);
			}
		    }	
		    m_level_afpi[iteration-1]->destGraph->worker[tg]->completedRegionQueue->addRegion(f,true);
		}
	    }
	    //fout.close();
	    ////////////////////////////////////////////////////////Push Pull Thread End////////////////////
	}
	else
	{
	    //fout << "Calling init "<< std::endl;
	    //fout.close();
	    init();
	}
    }

#else

    RGIter::RGIter(Vector<amrex::AsyncFillPatchIterator*> afpi, Vector<amrex::AsyncFillPatchIterator*> upper_afpi,
            amrex::MultiFab& dest, int  bG, double tm, int  ind, int  sc, int nc, int itr)
:
        itrGraph(afpi[itr-1]->destGraph),
        m_level_afpi(afpi),
        m_upper_level_afpi(upper_afpi),
        _dest(&dest),
        boxGrow(bG),
        time(tm),
        index(ind),
        scomp(sc),
        ncomp(nc),
        iteration(itr),
        implicit(true),
        ppteams(true),
        haveDepGraph(false),
        depGraph(NULL),
        getFireableTime(0.)
   {
        int myProc = amrex::ParallelDescriptor::MyProc();
        bool push = true;
        tid = perilla::tid();
        tg = perilla::wid();
        ntid= perilla::wtid();

        if(perilla::isCommunicationThread())
        {
            std::vector<RegionGraph*> flattenedGraphArray;
            Perilla::flattenGraphHierarchy(m_level_afpi[iteration-1]->m_amrlevel.parent->graphArray, flattenedGraphArray);
            while(true){
                Perilla::serviceMultipleGraphCommDynamic(flattenedGraphArray,true,perilla::tid());
                if( Perilla::numTeamsFinished == perilla::NUM_THREAD_TEAMS)
                {
		    flattenedGraphArray.clear();
		    Perilla::numTeamsFinished=0;
                    break;
                }
            }
        }else
{

        AsyncFillPatchIterator::initialSend(m_level_afpi, m_upper_level_afpi, boxGrow, time, index, scomp, ncomp, iteration);

        itrGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
        if(perilla::isMasterWorkerThread())
            m_level_afpi[iteration-1]->Reset();
        itrGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);

        if(ntid == perilla::NUM_THREADS_PER_TEAM-2)
        {
            int f;
            int level = m_level_afpi[iteration-1]->m_amrlevel.level;
            double dt = m_level_afpi[iteration-1]->m_amrlevel.parent->dtLevel(level);
            this->currentItr = 1;
            this->totalItr = 1;
            while(m_level_afpi[iteration-1]->destGraph->worker[tg]->completedRegionQueue->queueSize(true) != m_level_afpi[iteration-1]->destGraph->worker[tg]->totalTasks ||
                    m_level_afpi[iteration-1]->destGraph->worker[tg]->computedTasks != m_level_afpi[iteration-1]->destGraph->worker[tg]->totalTasks)
            {
                f = m_level_afpi[iteration-1]->destGraph->getFireableRegion(tg);
                if(f != -1)
                {
                    m_level_afpi[iteration-1]->Receive(this,*_dest,boxGrow,time,index,scomp,ncomp,f,true);
                    m_level_afpi[iteration-1]->destGraph->setFireableRegion(f);
                    if(m_level_afpi[iteration-1]->destGraph->worker[tg]->unfireableRegionQueue->queueSize(true) !=0 &&
                            m_level_afpi[iteration-1]->destGraph->worker[tg]->fireableRegionQueue->queueSize(true) < 8)
                        continue;
                }

                if(m_level_afpi[iteration-1]->destGraph->worker[tg]->computedRegionQueue->queueSize() != 0)
                {
                    f = m_level_afpi[iteration-1]->destGraph->worker[tg]->computedRegionQueue->removeRegion();

                    if(push & level == m_level_afpi[iteration-1]->m_amrlevel.parent->finestLevel() && iteration < m_level_afpi[iteration-1]->m_amrlevel.parent->nCycle(level))
                        m_level_afpi[iteration]->SendIntraLevel(*(this),boxGrow,time+dt,index,scomp,ncomp,iteration,f,true);

                    if(push & level < m_level_afpi[iteration-1]->m_amrlevel.parent->finestLevel())
                    {
                        for(int i=0; i < m_level_afpi[iteration-1]->m_amrlevel.parent->nCycle(level+1); i++)
                        {
                            m_upper_level_afpi[i]->SendInterLevel(this,boxGrow,time+(i*m_level_afpi[iteration-1]->m_amrlevel.parent->dtLevel(level+1)),index,scomp,ncomp,i+1,f,true);
                        }
                    }
                    m_level_afpi[iteration-1]->destGraph->worker[tg]->completedRegionQueue->addRegion(f,true);
                }
            }
        }
        else
        {
            //fout << "Calling init "<< std::endl;
            //fout.close();
            init();
        }
}

    }

#endif

    using namespace perilla;

    RGIter::~RGIter()
    {
	//fout.close();
    }

    void RGIter::init()
    {
	if(itrGraph->fabTiles.size() == 0)
	    tiling = false;
	else
	    tiling = true;

	int myProc = amrex::ParallelDescriptor::MyProc();
	if(implicit)
	{
	    if(!itrGraph->isGraphEmptyV2())
	    {
		currentRegion = itrGraph->getPulledFireableRegion();
		if(tiling)
		    totalItr = std::ceil( (1.0*itrGraph->fabTiles[currentRegion]->numTiles) / (perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS-1) );
		else
		    totalItr = 1;

		currentItr = 1;

		currentTile = 0;
		if(tiling)
		    for(currentTile = 0; currentTile < itrGraph->fabTiles[currentRegion]->numTiles; currentTile++)
			if(currentTile % (perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS-1) == ntid/*-perilla::NUM_COMM_THREADS*/)
			    break;
	    }
	    else
	    {
	    }
	}
	else
	{
	    if(!itrGraph->isGraphEmpty())
	    {
		if(haveDepGraph)
		    currentRegion = itrGraph->getAnyFireableRegion(*depGraph);
		else
		    currentRegion = itrGraph->getAnyFireableRegion();

		if(tiling)
		    totalItr = std::ceil( (1.0*itrGraph->fabTiles[currentRegion]->numTiles) / (perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS) );
		else
		    totalItr = 1;

		currentItr = 1;

		currentTile = 0;
		if(tiling)
		    for(currentTile = 0; currentTile < itrGraph->fabTiles[currentRegion]->numTiles; currentTile++)
			if(currentTile % (perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS) == ntid/*-perilla::NUM_COMM_THREADS*/)
			    break;
	    }
	    else
	    {
	    }
	}
    }

    //! Increment iterator to the next tile we own.
    void RGIter::operator++ ()
    {

	currentItr++;

	if(tiling)
	    for( (currentTile == itrGraph->fabTiles[currentRegion]->numTiles ? currentTile : ++currentTile); currentTile < itrGraph->fabTiles[currentRegion]->numTiles; currentTile++)
	    {
		if(implicit)
		{
		    if(currentTile % (perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS-1) == ntid/*-perilla::NUM_COMM_THREADS*/)
			break;
		}
		else
		{
		    if(currentTile % (perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS) == ntid/*-perilla::NUM_COMM_THREADS*/)
			break;
		}
	    }

	int myProc = amrex::ParallelDescriptor::MyProc();

	if( currentItr > totalItr )//&& currentTile == itrGraph->fabTiles[currentRegion]->numTiles)
	{
	    //if(WorkerThread::isTeamMasterThread(tid) )
	    //fout << "++B GEmpty " << itrGraph->isGraphEmpty(tg) << std::endl;

	    //fout << "++B CmpReg isGE " << (implicit?  itrGraph->isGraphEmptyV2(tg) : itrGraph->isGraphEmpty(tg)) << " CompleteQ "<< itrGraph->worker[tg]->nompletedRegionQueue->queueSize(true) << " totTasks " << itrGraph->worker[tg]->totalTasks << " FireQ "<< itrGraph->worker[tg]->fireableRegionQueue->queueSize(true) << " UnfireQ "<< itrGraph->worker[tg]->unfireableRegionQueue->queueSize(true) << std::endl;

	    if(implicit)
		itrGraph->regionComputed(currentRegion);
	    else
		itrGraph->finalizeRegion(currentRegion);

	    //if(WorkerThread::isTeamMasterThread(tid) )
	    //fout << "++A GEmpty " << itrGraph->isGraphEmpty(tg) << std::endl;

	    //fout << "++A CmpReg isGE " << (implicit?  itrGraph->isGraphEmptyV2(tg) : itrGraph->isGraphEmpty(tg)) << " CompleteQ "<< itrGraph->worker[tg]->completedRegionQueue->queueSize(true) << " totTasks " << itrGraph->worker[tg]->totalTasks << " FireQ "<< itrGraph->worker[tg]->fireableRegionQueue->queueSize(true) << " UnfireQ "<< itrGraph->worker[tg]->unfireableRegionQueue->queueSize(true) << std::endl;

	    if(implicit)
	    {
		if(!itrGraph->isGraphEmptyV2())
		{
		    currentRegion = itrGraph->getPulledFireableRegion();
		    if(tiling)
			totalItr = std::ceil( (1.0*itrGraph->fabTiles[currentRegion]->numTiles) / (perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS-1) );
		    else
			totalItr = 1;

		    currentItr = 1;

		    currentTile = 0;
		    if(tiling)
			for(currentTile = 0; currentTile < itrGraph->fabTiles[currentRegion]->numTiles; currentTile++)
			    if(currentTile % (perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS-1) == ntid/*-perilla::NUM_COMM_THREADS*/)
				break;
		}
		else
		{
		    //fout << "Graph is Empty" << std::endl;
		    //currentRegion = 0;
		    //currentTile = 0;
		}
	    }
	    else
	    {
		if(!itrGraph->isGraphEmpty())
		{
//		    double start_time_wtime = omp_get_wtime();	  

		    if(haveDepGraph)
			currentRegion = itrGraph->getAnyFireableRegion(*depGraph);
		    else
			currentRegion = itrGraph->getAnyFireableRegion();

//		    double end_time_wtime = omp_get_wtime();
//		    getFireableTime += end_time_wtime - start_time_wtime;

		    if(tiling)
			totalItr = std::ceil( (1.0*itrGraph->fabTiles[currentRegion]->numTiles) / (perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS) );
		    else
			totalItr = 1;

		    currentItr = 1;

		    currentTile = 0;
		    if(tiling)
			for(currentTile = 0; currentTile < itrGraph->fabTiles[currentRegion]->numTiles; currentTile++)
			    if(currentTile % (perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS) == ntid/*-perilla::NUM_COMM_THREADS*/)
				break;	      
		}
	    }
	}

	//fout << "++E Region " << currentRegion << " Tile " << currentTile << " numTile "<< itrGraph->fabTiles[currentRegion]->numTiles <<" tid " << tid << " myP " << myProc <<std::endl;
	//fout.close();
    }

    //! Is the iterator valid, are more regions to iterate over?
    bool RGIter::isValid ()
    {
        if(perilla::isCommunicationThread() ) return false;
	bool valid;
	bool do_remaining = true;

	int myProc = amrex::ParallelDescriptor::MyProc();

	if(implicit)
	{
	    if(ntid != perilla::NUM_THREADS_PER_TEAM-1)
	    {
		valid = !itrGraph->isGraphEmptyV2();
		if(valid)	      
		{
		    do_remaining = false;
		}
	    }

	    if(do_remaining)
	    {
		bool push = false;

		int f;
		int level = m_level_afpi[iteration-1]->m_amrlevel.level;
		double dt = m_level_afpi[iteration-1]->m_amrlevel.parent->dtLevel(level);
		this->currentItr = 1;
		this->totalItr = 1;

#if 0
		while(!itrGraph->isGraphEmpty())
		{
		    f = itrGraph->worker[tg]->computedRegionQueue->getFrontRegion(true);

		    if(push & level == m_level_afpi[iteration-1]->m_amrlevel.parent->finestLevel() && iteration < m_level_afpi[iteration-1]->m_amrlevel.parent->nCycle(level))
			m_level_afpi[iteration]->SendIntraLevel(this,boxGrow,time+dt,index,scomp,ncomp,iteration,f,false);
		    //else if(level == parent->finestLevel() && iteration == ncycle)
		    //SborderAFPI[0]->PushOnly(NUM_GROW, time+dt, State_Type, 0, NUM_STATE, f, tid, 0x02, 1);

		    if(push & level < m_level_afpi[iteration-1]->m_amrlevel.parent->finestLevel())
		    {
			for(int i=0; i < m_level_afpi[iteration-1]->m_amrlevel.parent->nCycle(level+1); i++)
			{
			    m_upper_level_afpi[i]->SendInterLevel(this,boxGrow,time+(i*m_level_afpi[iteration-1]->m_amrlevel.parent->dtLevel(level+1)),index,scomp,ncomp,i+1,f,false);
			    //upperLevel.SborderAFPI[i]->PushOnly(NUM_GROW, time+(i*parent->dtLevel(level+1)), State_Type, 0, NUM_STATE, f, tid, tuc, tempf, false);
			}
		    }		

		    itrGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS-1);
		    if(perilla::isMasterWorkerThread())
		    {
			f = itrGraph->worker[tg]->computedRegionQueue->removeRegion();
			itrGraph->worker[tg]->completedRegionQueue->addRegion(f,true);
		    }	  	
		}
#endif

		//m_level_afpi[iteration-1]->destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
		if(perilla::isMasterWorkerThread())
		{
		    m_level_afpi[iteration-1]->completeRegionGraphs();
#ifdef USE_PERILLA_ON_DEMAND		    
                    pthread_mutex_lock(&teamFinishLock);
                    Perilla::numTeamsFinished++;
                    pthread_mutex_unlock(&teamFinishLock);
#endif
		}
		valid = false;
	    }
	}
	else
	{
	    if(itrGraph->isGraphEmpty())
            {
		if(perilla::isMasterWorkerThread())
		{
		    itrGraph->finalizeRegionGraph();
#ifdef USE_PERILLA_ON_DEMAND		    
                    pthread_mutex_lock(&teamFinishLock);
                    Perilla::numTeamsFinished++;
                    pthread_mutex_unlock(&teamFinishLock);
#endif
                }
            }
	    valid = !(itrGraph->isGraphEmpty());
	}
	/*
	   itrGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
	   if(!isV && tg==0 && myProc==0)
	   if(WorkerThread::isTeamMasterThread(tid))
	   fout << " M " <<std::endl;
	//    else
	//  fout << " W " <<std::endl;
	itrGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
	 */

	/*
	   fout << "isValid Ending " << !(itrGraph->isGraphEmpty(tg)) << " tid " << tid <<std::endl;
	   fout.close();
	 */
	//if(!valid && tid == perilla::NUM_COMM_THREADS)
	{
	    //Perilla::getAnyFRTime += getFireableTime;
	    //if(itrGraph->graphID != -1)
		//Perilla::getAnyFRTimeSplit[itrGraph->graphID-1] += getFireableTime;
	    //if(myProc == 150 && itrGraph->graphID != -1)
	    //{
	    //  std::cout << "gID " << itrGraph->graphID << " getFRTime " << getFireableTime << std::endl;
	    //}
	}

	return  valid;
    }

    amrex::Box RGIter::tileBox()
    {

	int myProc = amrex::ParallelDescriptor::MyProc();
	//fout.open(std::to_string(myProc)+ "_" + std::to_string(tid) + ".txt", std::fstream::app);

	//fout << "nTls " << itrGraph->fabTiles[currentRegion]->numTiles << " cT " << currentTile << std::endl;

	if(currentTile == itrGraph->fabTiles[currentRegion]->numTiles)
	    //if( (currentTile % (perilla::NUM_THREADS_PER_TEAM-1) != ntid-1) )
	{
	    //fout << "invalidBox " << std::endl;
	    //fout.close();
	    return amrex::Box();
	}
	else
	{
	    //fout << "validBox tBxSize " << itrGraph->fabTiles[currentRegion]->tileBx.size() << std::endl;
	    //fout.close();
	    return   *(itrGraph->fabTiles[currentRegion]->tileBx[currentTile]);
	}
    }

    amrex::Box RGIter::validBox() const
    {
	return *(itrGraph->fabTiles[currentRegion]->validBx);
    }

    amrex::Box RGIter::tilebox()
    {
	return this->tileBox();
    }

    amrex::Box RGIter::growntilebox()
    {
	return this->tileBox();	
    }

    amrex::Box RGIter::growntilebox(int ng)
    {

	Box bx = this->tileBox();
	if(currentTile == itrGraph->fabTiles[currentRegion]->numTiles)
	    return bx;

	if (ng < -100) ng = 0;
	const Box& vbx = validBox();
	for (int d=0; d<BL_SPACEDIM; ++d) {
	    if (bx.smallEnd(d) == vbx.smallEnd(d)) {
		bx.growLo(d, ng);
	    }
	    if (bx.bigEnd(d) == vbx.bigEnd(d)) {
		bx.growHi(d, ng);
	    }
	}
	return bx;

    }

    amrex::Box RGIter::nodaltilebox(int dir)
    {
	BL_ASSERT(dir < BL_SPACEDIM);
	//BL_ASSERT(tile_array != 0);

	//Box bx((*tile_array)[currentIndex]);
	Box bx = this->tileBox();
	bx.convert(typ);
	const Box& vbx = this->validBox();
	const IntVect& Big = vbx.bigEnd();
	int d0, d1;
	if (dir < 0) {
	    d0 = 0;
	    d1 = BL_SPACEDIM-1;
	} else {
	    d0 = d1 = dir;
	}
	for (int d=d0; d<=d1; ++d) {
	    if (typ.cellCentered(d)) { // validbox should also be cell-centered in d-direction.
		bx.surroundingNodes(d);
		if (bx.bigEnd(d) <= Big[d]) {
		    bx.growHi(d,-1);
		}
	    }
	}
	return bx;
    }

    void RGIter::sync_workers()
    {

	if(implicit)
	    itrGraph->worker[tg]->l_barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS-1);
	else
	    itrGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);

    }
}
