
#include <AMReX_Box.H>

#include <RGIter.H>
#include <WorkerThread.H>
#include <PerillaConfig.H>

#include <Perilla.H>

#include <cmath>

namespace amrex{

RGIter::RGIter(RegionGraph* rg, bool enableAllTasks)
:
  itrGraph(rg),
  implicit(false),
  ppteams(true),
  typ(rg->typ),
  haveDepGraph(false),
  depGraph(NULL),
  getFireableTime(0.)
{

  tid = omp_get_thread_num();
  tg = WorkerThread::groupID(tid);
  ntid = WorkerThread::numaTID(tid);

  itrGraph->worker[tg]->barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS);
  if(WorkerThread::isTeamMasterThread(tid))
    itrGraph->Reset(tg);
  itrGraph->worker[tg]->barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS);

  if(enableAllTasks)
    itrGraph->enableAllRegions(tid);
  itrGraph->worker[tg]->barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS);

  init();
}

RGIter::RGIter(RegionGraph* rg, RegionGraph* drg, bool isDep)
:
  itrGraph(rg),
  implicit(false),
  ppteams(true),
  typ(rg->typ),
  haveDepGraph(isDep),
  depGraph(drg),
  getFireableTime(0.)
{
  //int myProc = amrex::ParallelDescriptor::MyProc();
  tid = omp_get_thread_num();
  tg = WorkerThread::groupID(tid);
  ntid = WorkerThread::numaTID(tid);

  itrGraph->worker[tg]->barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS);
  if(WorkerThread::isTeamMasterThread(tid))
    itrGraph->Reset(tg);
  itrGraph->worker[tg]->barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS);

  //if(enableAllTasks)
  //itrGraph->enableAllRegions(tid);
  //itrGraph->worker[tg]->barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-1);

  //if(myProc ==0 && WorkerThread::isTeamMasterThread(tid))
  //std::cout << "initializing RGIter hDG " << haveDepGraph << std::endl;

  init();

  //if(myProc ==0 && WorkerThread::isTeamMasterThread(tid))
  //std::cout << "initialized RGIter " << std::endl;

}


RGIter::RGIter(amrex::AsyncFillPatchIterator* afpi, bool enableAllTasks)
  :
  itrGraph(afpi->destGraph),
  implicit(false),
  ppteams(true),
  typ(afpi->destGraph->typ),
  haveDepGraph(false),
  depGraph(NULL),
  getFireableTime(0.)
{
  tid = omp_get_thread_num();
  tg = WorkerThread::groupID(tid);
  ntid = WorkerThread::numaTID(tid);

  itrGraph->worker[tg]->barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS);
  if(WorkerThread::isTeamMasterThread(tid))
    afpi->Reset(tg);
  itrGraph->worker[tg]->barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS);

  if(enableAllTasks)
    itrGraph->enableAllRegions(tid);
  itrGraph->worker[tg]->barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS);

  init();
}

RGIter::RGIter(amrex::Array<amrex::AsyncFillPatchIterator*> afpi, amrex::Array<amrex::AsyncFillPatchIterator*> upper_afpi, 
	       amrex::MultiFab& dest, int  bG, double tm, int  ind, int  sc, int nc, int itr)
   :
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
   typ(afpi[itr-1]->destGraph->typ),
   haveDepGraph(false),
   depGraph(NULL),
   getFireableTime(0.)
 {
    int myProc = amrex::ParallelDescriptor::MyProc();

    bool push = false;

    tid = omp_get_thread_num();
    tg = WorkerThread::groupID(tid);
    ntid = WorkerThread::numaTID(tid);
    
    itrGraph->worker[tg]->barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS);
    if(WorkerThread::isTeamMasterThread(tid))
      m_level_afpi[iteration-1]->Reset(tg);
    itrGraph->worker[tg]->barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS);

    //fout.open(std::to_string(myProc)+ "_" + std::to_string(tid) + ".txt", std::fstream::app);
    //fout << "Spliting up" << std::endl;

    if(ntid == PerillaConfig::NUM_THREADS_PER_TEAM-1)
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

	    //if(myProc == 0 && f != -1)
	    //std::cout << "Got region "<< f << " for Pulling unfireQ Sz " << m_level_afpi[iteration-1]->destGraph->worker[tg]->unfireableRegionQueue->queueSize(true) << " computedQ Sz " << m_level_afpi[iteration-1]->destGraph->worker[tg]->computedRegionQueue->queueSize() << " completedQ Sz " << m_level_afpi[iteration-1]->destGraph->worker[tg]->completedRegionQueue->queueSize(true) << " tot " << m_level_afpi[iteration-1]->destGraph->worker[tg]->totalTasks<< std::endl;


	    if(f != -1)
	      {
		m_level_afpi[iteration-1]->Receive(this,dest,boxGrow,time,index,scomp,ncomp,f,tid,true);
		//m_level_afpi[iteration-1]->PullOnly(NUM_GROW, time, State_Type, 0, NUM_STATE, f, tid, true);
		m_level_afpi[iteration-1]->destGraph->setFireableRegion(tg, f);
		if(m_level_afpi[iteration-1]->destGraph->worker[tg]->unfireableRegionQueue->queueSize(true) !=0 && 
		   m_level_afpi[iteration-1]->destGraph->worker[tg]->fireableRegionQueue->queueSize(true) < 2)
		  continue;
	      }

	    if(m_level_afpi[iteration-1]->destGraph->worker[tg]->computedRegionQueue->queueSize() != 0)
	      {
		f = m_level_afpi[iteration-1]->destGraph->worker[tg]->computedRegionQueue->removeRegion();

		//fout << "Pushing for region "<< f << " lvl " << level << " itr " << iteration << std::endl;

		if(push & level == m_level_afpi[iteration-1]->m_amrlevel.parent->finestLevel() && iteration < m_level_afpi[iteration-1]->m_amrlevel.parent->nCycle(level))
		  m_level_afpi[iteration]->SendIntraLevel(this,boxGrow,time+dt,index,scomp,ncomp,iteration,f,tid,true);
		    //m_level_afpi[iteration]->PushOnly(NUM_GROW, time+dt, State_Type, 0, NUM_STATE, f, tid, 0x02, 1, true);
		//else if(level == parent->finestLevel() && iteration == ncycle)
		//SborderAFPI[0]->PushOnly(NUM_GROW, time+dt, State_Type, 0, NUM_STATE, f, tid, 0x02, 1);
	  
		if(push & level < m_level_afpi[iteration-1]->m_amrlevel.parent->finestLevel())
		  {
		    for(int i=0; i < m_level_afpi[iteration-1]->m_amrlevel.parent->nCycle(level+1); i++)
		      {
			m_upper_level_afpi[i]->SendInterLevel(this,boxGrow,time+(i*m_level_afpi[iteration-1]->m_amrlevel.parent->dtLevel(level+1)),index,scomp,ncomp,i+1,f,tid,true);
			//m_upper_level_afpi[i]->PushOnly(NUM_GROW, time+(i*m_level_afpi[iteration-1]->m_amrlevel.parent->dtLevel(level+1)), State_Type, 0, NUM_STATE, f, tid, tuc, tempf, true);
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
  /*
  fout.open(std::to_string(myProc)+ "_" + std::to_string(tid) + ".txt", std::fstream::app);

  fout << "BR Starting Init: isGE " << itrGraph->isGraphEmpty(tg) << " CompleteQ "<< itrGraph->worker[tg]->completedRegionQueue->queueSize(true) << " totTasks " << itrGraph->worker[tg]->totalTasks << " FireQ "<< itrGraph->worker[tg]->fireableRegionQueue->queueSize(true) << " UnfireQ "<< itrGraph->worker[tg]->unfireableRegionQueue->queueSize(true) << std::endl;

  fout.close();
  */

  //fout.open(std::to_string(myProc)+ "_" + std::to_string(tid) + ".txt", std::fstream::app);
  
  if(implicit)
    {

      //fout << "Starting Init: isGE " << itrGraph->isGraphEmptyV2(tg) << " CompleteQ "<< itrGraph->worker[tg]->completedRegionQueue->queueSize(true) << " totTasks " << itrGraph->worker[tg]->totalTasks << " FireQ "<< itrGraph->worker[tg]->fireableRegionQueue->queueSize(true) << " UnfireQ "<< itrGraph->worker[tg]->unfireableRegionQueue->queueSize(true) << std::endl;
      
      if(!itrGraph->isGraphEmptyV2(tg))
	{
	  currentRegion = itrGraph->getPulledFireableRegion(tid);
	  if(tiling)
	    totalItr = std::ceil( (1.0*itrGraph->fabTiles[currentRegion]->numTiles) / (PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS-1) );
	  else
	    totalItr = 1;
	  
	  currentItr = 1;
	  
	  currentTile = 0;
	  if(tiling)
	    for(currentTile = 0; currentTile < itrGraph->fabTiles[currentRegion]->numTiles; currentTile++)
	      if(currentTile % (PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS-1) == ntid-PerillaConfig::NUM_COMM_THREADS)
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

      //fout << "Starting Init: isGE " << itrGraph->isGraphEmpty(tg) << " CompleteQ "<< itrGraph->worker[tg]->completedRegionQueue->queueSize(true) << " totTasks " << itrGraph->worker[tg]->totalTasks << " FireQ "<< itrGraph->worker[tg]->fireableRegionQueue->queueSize(true) << " UnfireQ "<< itrGraph->worker[tg]->unfireableRegionQueue->queueSize(true) << std::endl;
      
      if(!itrGraph->isGraphEmpty(tg))
	{
	  //if(tid==1)
	  //std::cout << "Getting AFR at " << itrGraph->graphID << std::endl;
	  double start_time_wtime = omp_get_wtime();	  
	  if(haveDepGraph)
	    currentRegion = itrGraph->getAnyFireableRegion(*depGraph, tid);
	  else
	    currentRegion = itrGraph->getAnyFireableRegion(tid);

	  double end_time_wtime = omp_get_wtime();
	  getFireableTime += end_time_wtime - start_time_wtime;

	  //if(tid==1)
	  //std::cout << "Got AFR at " << itrGraph->graphID << " f " << currentRegion<< std::endl;
	  /*
	  if(tid==1 || tid ==13)
	    {
	      std::cout << "init f " << currentRegion << " gID " << itrGraph->graphID << " nFabs " << itrGraph->fabTiles.size() << " tTsk " << itrGraph->worker[tg]->totalTasks << " cmpltRQ " << itrGraph->worker[tg]->completedRegionQueue->queueSize(true) << " fireRQ " << itrGraph->worker[tg]->fireableRegionQueue->queueSize(true) << " unfireRQ " << itrGraph->worker[tg]->unfireableRegionQueue->queueSize(true) << " myP " << myProc << " tg " << tg << std::endl;
	      std::cout << "init f " << currentRegion << " gID " << itrGraph->graphID << " nTiles " << itrGraph->fabTiles[currentRegion]->numTiles << " myP " << myProc << " tg " << tg << std::endl;
	    }
	  */
	  if(tiling)
	    totalItr = std::ceil( (1.0*itrGraph->fabTiles[currentRegion]->numTiles) / (PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS) );
	  else
	    totalItr = 1;
	  
	  currentItr = 1;
	  
	  currentTile = 0;
	  if(tiling)
	    for(currentTile = 0; currentTile < itrGraph->fabTiles[currentRegion]->numTiles; currentTile++)
	      if(currentTile % (PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS) == ntid-PerillaConfig::NUM_COMM_THREADS)
		break;
	}
      else
	{
	  //fout << "Graph is Empty" << std::endl;
	  //currentRegion = 0;
	  //currentTile = 0;
	}
    }
 
  //fout << "Init: Region " << currentRegion << " numTile "<< itrGraph->fabTiles[currentRegion]->numTiles <<" tid " << tid << " myP " << myProc <<std::endl;
  //fout.close();

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
	  if(currentTile % (PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS-1) == ntid-PerillaConfig::NUM_COMM_THREADS)
	    break;
	}
      else
	{
	  if(currentTile % (PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS) == ntid-PerillaConfig::NUM_COMM_THREADS)
	    break;
	}
    }

  int myProc = amrex::ParallelDescriptor::MyProc();
  //fout.open(std::to_string(myProc)+ "_" + std::to_string(tid) + ".txt", std::fstream::app);
  //if(WorkerThread::isTeamMasterThread(tid))
  //fout << "++S Region " << currentRegion << " Tile " << currentTile << " numTile "<< itrGraph->fabTiles[currentRegion]->numTiles <<" tid " << tid << " myP " << myProc <<std::endl;


  if( currentItr > totalItr )//&& currentTile == itrGraph->fabTiles[currentRegion]->numTiles)
    {
      //if(WorkerThread::isTeamMasterThread(tid) )
      //fout << "++B GEmpty " << itrGraph->isGraphEmpty(tg) << std::endl;
      	
      //fout << "++B CmpReg isGE " << (implicit?  itrGraph->isGraphEmptyV2(tg) : itrGraph->isGraphEmpty(tg)) << " CompleteQ "<< itrGraph->worker[tg]->completedRegionQueue->queueSize(true) << " totTasks " << itrGraph->worker[tg]->totalTasks << " FireQ "<< itrGraph->worker[tg]->fireableRegionQueue->queueSize(true) << " UnfireQ "<< itrGraph->worker[tg]->unfireableRegionQueue->queueSize(true) << std::endl;

      if(implicit)
	itrGraph->regionComputed(currentRegion,tg,ntid,tid);
      else
	itrGraph->completeRegion(currentRegion,tg,ntid,tid);

      //if(WorkerThread::isTeamMasterThread(tid) )
      //fout << "++A GEmpty " << itrGraph->isGraphEmpty(tg) << std::endl;

      //fout << "++A CmpReg isGE " << (implicit?  itrGraph->isGraphEmptyV2(tg) : itrGraph->isGraphEmpty(tg)) << " CompleteQ "<< itrGraph->worker[tg]->completedRegionQueue->queueSize(true) << " totTasks " << itrGraph->worker[tg]->totalTasks << " FireQ "<< itrGraph->worker[tg]->fireableRegionQueue->queueSize(true) << " UnfireQ "<< itrGraph->worker[tg]->unfireableRegionQueue->queueSize(true) << std::endl;

  if(implicit)
    {
      if(!itrGraph->isGraphEmptyV2(tg))
	{
	  currentRegion = itrGraph->getPulledFireableRegion(tid);
	  if(tiling)
	    totalItr = std::ceil( (1.0*itrGraph->fabTiles[currentRegion]->numTiles) / (PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS-1) );
	  else
	    totalItr = 1;
	  
	  currentItr = 1;
	  
	  currentTile = 0;
	  if(tiling)
	    for(currentTile = 0; currentTile < itrGraph->fabTiles[currentRegion]->numTiles; currentTile++)
	      if(currentTile % (PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS-1) == ntid-PerillaConfig::NUM_COMM_THREADS)
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
      if(!itrGraph->isGraphEmpty(tg))
	{
	  double start_time_wtime = omp_get_wtime();	  

	  if(haveDepGraph)
	    currentRegion = itrGraph->getAnyFireableRegion(*depGraph, tid);
	  else
	    currentRegion = itrGraph->getAnyFireableRegion(tid);

	  double end_time_wtime = omp_get_wtime();
	  getFireableTime += end_time_wtime - start_time_wtime;
	  
	  if(tiling)
	    totalItr = std::ceil( (1.0*itrGraph->fabTiles[currentRegion]->numTiles) / (PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS) );
	  else
	    totalItr = 1;

	  currentItr = 1;

	  currentTile = 0;
	  if(tiling)
	    for(currentTile = 0; currentTile < itrGraph->fabTiles[currentRegion]->numTiles; currentTile++)
	      if(currentTile % (PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS) == ntid-PerillaConfig::NUM_COMM_THREADS)
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

  bool valid;
  bool do_remaining = true;

  int myProc = amrex::ParallelDescriptor::MyProc();
  /*
  fout.open(std::to_string(myProc)+ "_" + std::to_string(tid) + ".txt", std::fstream::app);
  fout << "isValid Starting " << " CompleteQ "<< itrGraph->worker[tg]->completedRegionQueue->queueSize(true) << " totTasks " << itrGraph->worker[tg]->totalTasks << " FireQ "<< itrGraph->worker[tg]->fireableRegionQueue->queueSize(true) << " UnfireQ "<< itrGraph->worker[tg]->unfireableRegionQueue->queueSize(true) << std::endl;

  bool isV = !(itrGraph->isGraphEmpty(tg));

  fout << "isV "<< isV << " tid " << tid <<" tg " << tg<< " myP " << myProc <<std::endl;
  */

  if(implicit)
    {
      if(ntid != PerillaConfig::NUM_THREADS_PER_TEAM-1)
	{
	  valid = !itrGraph->isGraphEmptyV2(tg);
	  if(valid)	      
	    do_remaining = false;
	}

      if(do_remaining)
	{
	  bool push = false;

	  int f;
	  int level = m_level_afpi[iteration-1]->m_amrlevel.level;
	  double dt = m_level_afpi[iteration-1]->m_amrlevel.parent->dtLevel(level);
	  this->currentItr = 1;
	  this->totalItr = 1;

	  while(!itrGraph->isGraphEmpty(tg))
	    {

	      f = itrGraph->worker[tg]->computedRegionQueue->getFrontRegion(true);
	      
	      if(push & level == m_level_afpi[iteration-1]->m_amrlevel.parent->finestLevel() && iteration < m_level_afpi[iteration-1]->m_amrlevel.parent->nCycle(level))
		m_level_afpi[iteration]->SendIntraLevel(this,boxGrow,time+dt,index,scomp,ncomp,iteration,f,tid,false);
	      //else if(level == parent->finestLevel() && iteration == ncycle)
	      //SborderAFPI[0]->PushOnly(NUM_GROW, time+dt, State_Type, 0, NUM_STATE, f, tid, 0x02, 1);

	      if(push & level < m_level_afpi[iteration-1]->m_amrlevel.parent->finestLevel())
		  {
		    for(int i=0; i < m_level_afpi[iteration-1]->m_amrlevel.parent->nCycle(level+1); i++)
		      {
			m_upper_level_afpi[i]->SendInterLevel(this,boxGrow,time+(i*m_level_afpi[iteration-1]->m_amrlevel.parent->dtLevel(level+1)),index,scomp,ncomp,i+1,f,tid,false);
			//upperLevel.SborderAFPI[i]->PushOnly(NUM_GROW, time+(i*parent->dtLevel(level+1)), State_Type, 0, NUM_STATE, f, tid, tuc, tempf, false);
		      }
		  }		
	
	      itrGraph->worker[tg]->barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS);
	      if(WorkerThread::isTeamMasterThread(tid))
		{
		  f = itrGraph->worker[tg]->computedRegionQueue->removeRegion();
		  itrGraph->worker[tg]->completedRegionQueue->addRegion(f,true);
		}	  	
	    }
    

	  //m_level_afpi[iteration-1]->destGraph->worker[tg]->barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-1);
	  if(WorkerThread::isTeamMasterThread(tid))
	    m_level_afpi[iteration-1]->completeRegionGraphs(tg);

	  
	  valid = false;
	}
    }
  else
    {
      if(itrGraph->isGraphEmpty(tg))
	if(WorkerThread::isTeamMasterThread(tid))
	  {
	    //fout << " MB " << std::endl;
	    itrGraph->completeRegionGraph(tg);
	    //fout << " MA " << std::endl;
	  }
      valid = !(itrGraph->isGraphEmpty(tg));
    }
  /*
  itrGraph->worker[tg]->barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-1);
  if(!isV && tg==0 && myProc==0)
    if(WorkerThread::isTeamMasterThread(tid))
      fout << " M " <<std::endl;
  //    else
  //  fout << " W " <<std::endl;
  itrGraph->worker[tg]->barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-1);
  */

  /*
  fout << "isValid Ending " << !(itrGraph->isGraphEmpty(tg)) << " tid " << tid <<std::endl;
  fout.close();
  */
  if(!valid && tid == PerillaConfig::NUM_COMM_THREADS)
    {
      Perilla::getAnyFRTime += getFireableTime;
      if(itrGraph->graphID != -1)
	Perilla::getAnyFRTimeSplit[itrGraph->graphID-1] += getFireableTime;
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
    //if( (currentTile % (PerillaConfig::NUM_THREADS_PER_TEAM-1) != ntid-1) )
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
    itrGraph->worker[tg]->l_barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS-1);
  else
    itrGraph->worker[tg]->barr->sync(PerillaConfig::NUM_THREADS_PER_TEAM-PerillaConfig::NUM_COMM_THREADS);

  }

}
