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


void Perilla::serviceMultipleGraphCommDynamic(std::vector<std::vector<RegionGraph*> > graphArrayHierarchy, bool cpyAcross, int tid)
{
    int tg = WorkerThread::perilla_wid();
    int np = ParallelDescriptor::NProcs();
    int myProc = ParallelDescriptor::MyProc();
    int graphFinishCnt = 0;
    int nGraphs;
    bool doublechecked = false;

    double maxltime=0;
    double minltime=10;
    double avgltime=0;
    double numloops=0;
    double ltime,lstime,letime;

    int gCnt=0;
    for(int l=0; l<graphArrayHierarchy.size(); l++) gCnt+= graphArrayHierarchy[l].size();
    std::vector<RegionGraph*> graphArray;
    for(int l=0; l<graphArrayHierarchy.size(); l++)
        for(int g=0; g<graphArrayHierarchy[l].size(); g++)
            graphArray.push_back(graphArrayHierarchy[l][g]);

    while(true)
    {   
        //lstime = omp_get_wtime();
        for(int g=0; g<graphArray.size(); g++)
        {   
            nGraphs = graphArray.size();
            //if(graphArray[g]->totalFinishes < perilla::NUM_THREAD_TEAMS)
            {   
                /*try{
                  if(graphArray[g]->assocMF == 0)
                  std::cout<<"Processing Graph with NULL MF "<<g<<" ";
                  }catch (const std::exception& e) {
                  std::cout<<"Processing Graph with NULL MF "<<g<<" ";
                  }*/
                //if(graphArray[g]->graphID==1)
                //std::cout<<"Processing Local Req Graph "<<g+1 << " tg " <<tg <<std::endl;
                serviceLocalRequests(graphArray[g], tg);
                if(cpyAcross)
                {   
                    //if(graphArray[g]->graphID==13)
                    //std::cout<<"Processing Local GridCopy Req Graph "<< g+1 << " tg " << tg <<std::endl;
                    serviceLocalGridCopyRequests(graphArray,g,tg);
                }
                if(np > 1)//if(tg==0)
                {   
                    serviceRemoteRequests(graphArray[g],g,nGraphs);
                    if(cpyAcross)
                    {   
                        //resetRemoteGridCopyRequests(graphArray,g,nGraphs,tg);
                        if(tg==0)
                            serviceRemoteGridCopyRequests(graphArray,g,nGraphs,tg);
                    }
                }
            }
        }
        
        
        if( Perilla::numTeamsFinished == perilla::NUM_THREAD_TEAMS)
        {   
            if(doublechecked) // double check if there are still something to send
                return;
            else
                doublechecked = true;
        }

        //std::cout<<"Teams Completed "<< Perilla::numTeamsFinished << " tid "<< tid << " myProc " << myProc <<std::endl;

        //letime = omp_get_wtime();
        numloops++;
        //ltime = letime - lstime;

        avgltime += ltime;
        //if(ltime < minltime)
        //minltime = ltime;
        //if(ltime > maxltime)
        //maxltime = ltime;
    } // while(true)

    //nGraphs = graphArray.size();
    //if(tg==0)
    //for(int g=0; g<nGraphs; g++)
    //{
        //ParallelDescriptor::Barrier("serviceMultipleGraph-1");
        //graphArray[g]->graphTeardown(tg);
        //graphArray[g]->workerTeardown(tg);
        //ParallelDescriptor::Barrier("serviceMultipleGraph-2");
    //}

} // serviceMultipleGraphCommDynamic


void Perilla::multifabCopyPush(RegionGraph* destGraph, RegionGraph* srcGraph, amrex::MultiFab* mfDst, amrex::MultiFab* mfSrc, int f, int dstcomp, int srccomp, int nc, int ng, int ngsrc, bool singleT)
{

    //double start_time_wtime = omp_get_wtime();

    if(nc<1) cout <<"MULTIFAB_COPY_C: nc must be >= 1"<< endl;
    if(mfDst->nComp() < (dstcomp-1)) cout <<"MULTIFAB_COPY_C: nc too large for dst multifab"<< endl;
    if(mfSrc->nComp() < (srccomp-1)) cout <<"MULTIFAB_COPY_C: nc too large for src multifab"<< endl;

    //mTeams = false; 

//    if(np==1)
      //multifabCopyPush_1Team(destGraph,srcGraph,mfDst,mfSrc,f,dstcomp,srccomp,nc,ng,ngsrc,singleT);
/*    else if(mTeams)
      {
        if(WorkerThread::isLocPPTID(tid))
          multifabCopyLocPush(destGraph,srcGraph,mfDst,mfSrc,f,tid,dstcomp,srccomp,nc,ng,ngsrc);
        else
          multifabCopyRmtPush(destGraph,srcGraph,mfDst,mfSrc,f,tid,dstcomp,srccomp,nc,ng,ngsrc);
      }
    else
      multifabCopyPush_1Team(destGraph,srcGraph,mfDst,mfSrc,f,tid,dstcomp,srccomp,nc,ng,ngsrc,singleT);
*/

    if(!singleT)
      srcGraph->worker[perilla::wid()]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);

    //double end_time_wtime = omp_get_wtime();
    //if(ntid==0)
      //Perilla::getPPPTimeSplit[2] += end_time_wtime - start_time_wtime;
}


void Perilla::multifabCopyPush(RegionGraph* destGraph, RegionGraph* srcGraph, amrex::MultiFab* mfDst, amrex::MultiFab* mfSrc, int f, bool singleT)
  {
    multifabCopyPush(destGraph, srcGraph, mfDst, mfSrc, f, 1, 1, 1, 0, 0, singleT);
  }

  void Perilla::multifabCopyPush_1Team(RegionGraph* destGraph, RegionGraph* srcGraph, amrex::MultiFab* mfDst, amrex::MultiFab* mfSrc, int f, int dstcomp, int srccomp, int nc, int ng, int ngsrc, bool singleT)
  {
    int ntid = perilla::wtid();
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
              cpSrc->l_con.scpy[i].pQueue.enqueue(cpSrc->l_con.scpy[i].recycleQueue.dequeue(true),true);
            pthread_mutex_unlock(&(cpSrc->l_con.sLock));
          }
        else
          { 
            if(ntid == 0)
              pthread_mutex_lock(&(cpSrc->l_con.sLock));      
            srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
            
            //std::ofstream fout;
            //fout.open(std::to_string(myProc)+ "_" + std::to_string(tid) + ".txt", std::fstream::app);
            for(int i=0; i<cpSrc->l_con.nscpy; i++)
              if((i%(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS)) == ntid)
                {
                  Package* sndPackage = cpSrc->l_con.scpy[i].recycleQueue.getFront(true);
                  mfSrc->m_fabs_v[f]->copyToMem(cpSrc->l_con.scpy[i].sbx,srccomp,nc,sndPackage->databuf);
                  /*
                  for(int ii=0; ii < sndPackage->bufSize; ii++)
                    if(sndPackage->databuf[ii] == 0)
                      fout << "MFCPush loc zero at " << f << " i " << i << " ii " << ii << " sbx "<< cpSrc->l_con.scpy[i].sbx << std::endl;
                  */
                }

            //fout.close();

            srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
            if(ntid==0)
              {
                for(int i=0;i<cpSrc->l_con.nscpy; i++)
                  cpSrc->l_con.scpy[i].pQueue.enqueue(cpSrc->l_con.scpy[i].recycleQueue.dequeue(true),true);
                pthread_mutex_unlock(&(cpSrc->l_con.sLock));
              }
            srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
          }

        int np = amrex::ParallelDescriptor::NProcs();
        if(np == 1)
          return;

        //if(myProc==26 && srcGraph->graphID==18  && ntid == 0)
        //std::cout << "Notw its sgID 18,"<< f <<" turn lets see " << cpSrc->r_con.nsnd <<std::endl;

        //if(myProc==28 && srcGraph->graphID==18  && ntid == 0)
        //std::cout << "Notw its sgID 18,"<< f <<" turn lets see " << cpSrc->r_con.nsnd <<std::endl;

        //if(srcGraph->graphID==18 && f ==316)   
        //BL_ASSERT(cpSrc->r_con.nsnd == 177);
        if(singleT)
          {
            pthread_mutex_lock(&(cpSrc->r_con.sndLock));
            for(int i=0; i<cpSrc->r_con.nsnd; i++)
              {
                Package* sndPackage = cpSrc->r_con.snd[i].recycleQueue.dequeue(true);
                mfSrc->m_fabs_v[f]->copyToMem(cpSrc->r_con.snd[i].sbx,srccomp,nc,sndPackage->databuf);
                sndPackage->notified = false;
                sndPackage->notified = false;
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
            srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);

            for(int i=0; i<cpSrc->r_con.nsnd; i++)
              if((i%(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS)) == ntid)
                {
                  Package* sndPackage = cpSrc->r_con.snd[i].recycleQueue.dequeue(true);
                  mfSrc->m_fabs_v[f]->copyToMem(cpSrc->r_con.snd[i].sbx,srccomp,nc,sndPackage->databuf);
                  sndPackage->notified = false;
                  sndPackage->notified = false;
                  cpSrc->r_con.snd[i].pQueue.enqueue(sndPackage,true);
                }

            //fout.close();         
            srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
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
            srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
          }
      } // if(!(*mfDst == *mfSrc))                                                                                                                    
  } // multifabCopyPush

  void Perilla::fillBoundaryPull_1Team(RegionGraph* graph, amrex::MultiFab& mf, int f)
  {
    int myProc = amrex::ParallelDescriptor::MyProc();
    int mfi = mf.IndexArray()[f];

    int nComp = mf.nComp();
    int tg= perilla::wid();
    int ntid = perilla::wtid();

    if(ntid==0)
      pthread_mutex_lock(&(graph->lMap[f]->l_con.dLock));
    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS); // Barrier to synchronize team threads    

    if(perilla::LAZY_PUSH)
      { }
    else
      {
        if(perilla::UNPACKING_FINEGRAIN)
          {}
        else
          {
            for(int i=0; i<graph->lMap[f]->l_con.ndcpy; i++)
              if( (i%(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS)) == ntid)
                {
                  Package *dPackage = graph->lMap[f]->l_con.dcpy[i].pQueue.getFront(true);
                  /*
                  for(int d=0; d<dPackage->bufSize; d++)
                    if(dPackage->databuf[d] == 0)
                      {
                        //std::cout<< "in fbPull Reciving 0 for f "<< f <<std::endl;
                        //BL_ASSERT(dPackage->databuf[d] != 0);
                      }
                  */
                  /*
                  if(f==0)
                  //if(graph->lMap[f]->l_con.dcpy[i].dbx.smallEnd() == graph->lMap[f]->l_con.dcpy[i].dbx.bigEnd())
                  //if(graph->lMap[f]->l_con.dcpy[i].dbx.smallEnd(0)==-1 && graph->lMap[f]->l_con.dcpy[i].dbx.smallEnd(1)==-1 && graph->lMap[f]->l_con.dcpy[i].dbx.smallEnd(2)==4)
                      std::cout<< "Corner Pull for f "<< f << " data0 " <<dPackage->databuf[0]<< " size " <<dPackage->bufSize <<" se " <<graph->lMap[f]->l_con.dcpy[i].dbx.smallEnd()<<std::endl;
                  */
                  /*
                  if(mfi==0)
                    {
                      std::cout<< "LPull " << i <<std::endl;
                      for(int d=0; d<dPackage->bufSize; d++)
                        std::cout << dPackage->databuf[d] << " ";
                      std::cout << std::endl;
                    }
                  */
                  mf.m_fabs_v[f]->copyFromMem(graph->lMap[f]->l_con.dcpy[i].dbx,0,nComp,dPackage->databuf);
                }
          } // if(UNPACKING_FINEGRAIN) - else
      } // if(LAZY_PUSH) - else

    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS); // Barrier to synchronize team threads

    if(ntid==0)
      {
        for(int i=0; i<graph->lMap[f]->l_con.ndcpy; i++)
          {
            graph->lMap[f]->l_con.dcpy[i].recycleQueue.enqueue( graph->lMap[f]->l_con.dcpy[i].pQueue.dequeue(true),true );
          }

        graph->lMap[f]->l_con.firingRuleCnt = graph->lMap[f]->l_con.firingRuleCnt - graph->lMap[f]->l_con.ndcpy;


        graph->lMap[f]->l_con.scpyCnt = 0;
        for(int i=0; i<graph->lMap[f]->l_con.ndcpy; i++)
          if(graph->lMap[f]->l_con.dcpy[i].pQueue.queueSize(true) >= 1)
            {
              graph->lMap[f]->l_con.firingRuleCnt++;
            }

        pthread_mutex_unlock(&(graph->lMap[f]->l_con.dLock));
      }
    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS); // Barrier to synchronize team threads

    int np = amrex::ParallelDescriptor::NProcs();
    if (np==1) return;
    if(ntid==0)
      {
        pthread_mutex_lock(&(graph->rMap[f]->r_con.rcvLock));
        pthread_mutex_lock(&(graph->lMap[f]->r_con.rcvLock));
      }
    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS); // Barrier to synchronize team threads

    for(int i=0; i<graph->lMap[f]->r_con.nrcv; i++)
      if( (i%(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS)) == ntid)
        {
          Package *rcvMetaPackage = graph->rMap[f]->r_con.rcv[i].pQueue.dequeue(true);
          rcvMetaPackage->completed = false;
          rcvMetaPackage->served = false;
          rcvMetaPackage->request =0;// MPI_REQUEST_NULL;
          graph->rMap[f]->r_con.rcv[i].recycleQueue.enqueue(rcvMetaPackage,true);
          Package *rcvPackage = graph->lMap[f]->r_con.rcv[i].pQueue.dequeue(true);

          mf.m_fabs_v[f]->copyFromMem(graph->lMap[f]->r_con.rcv[i].dbx,0,nComp,rcvPackage->databuf);
          rcvPackage->completed = false;
          rcvPackage->notified = false;
          graph->lMap[f]->r_con.rcv[i].recycleQueue.enqueue(rcvPackage,true);
        }
    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS); // Barrier to synchronize team threads

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
    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS); // Barrier to synchronize team threads
  } // fillBoundaryPull


#if 0
#endif


#if 0
Array<const FabArrayBase::CopyComTagsContainer*> send_cctc;
Array<int> send_pr;
Array<const FabArrayBase::CopyComTagsContainer*> recv_cctc;
Array<int> recv_pr;


void Perilla::multifabExtractCopyAssoc(RegionGraph* gDst, RegionGraph* gSrc, const MultiFab& mfDst, const MultiFab& mfSrc, int nc, int ng, int ngSrc, const Periodicity& period)
{
#if 1
    int myProc = ParallelDescriptor::MyProc();
    int np = ParallelDescriptor::NProcs();
    try{
	if(true)//if(!(*mfSrc == *mfDst))
	{
#ifdef USE_PERILLA_PTHREADS
//	    if(perilla::isMasterThread())
#endif
            {
		if(ng > mfDst.nGrow()) cout <<"MULTIFAB_COPY_C: ng > mfDst.nGrow not supported in parallel copy"<< endl;
		if(ngSrc > mfSrc.nGrow()) cout <<"MULTIFAB_COPY_C: ngSrc > mfSrc.nGrow"<< endl;
		if(ngSrc > 0)
		{
		    // To be implemented
		    //do i = 1, nboxes(msrc%la)
		    //  call push_back(bl, grow(box_nodalize(get_box(msrc%la,i),msrc%nodal),lngsrc))
		    //end do
		    //call build(batmp, bl, sort = .false.)
		    //call destroy(bl)
		    //call build(lasrctmp, batmp, boxarray_bbox(batmp), explicit_mapping = get_proc(msrc%la))
		    //call destroy(batmp)
		    //call build(msrctmp, lasrctmp, nc = lnc, ng = 0)
		    //pmfsrc => msrctmp
		    assert(false);
		}
		if(np > 1)
		{
		    if(gSrc->sCopyMapHead == 0)
			gSrc->sCopyMapHead = new CopyMap();
		    else
		    {
			CopyMap *tmpCopyMap = new CopyMap();
			tmpCopyMap->next = gSrc->sCopyMapHead;
			gSrc->sCopyMapHead = tmpCopyMap;
		    }
		    if(gDst->rCopyMapHead == 0)
			gDst->rCopyMapHead = new CopyMap();
		    else
		    {
			CopyMap *tmpCopyMap = new CopyMap();
			tmpCopyMap->next = gDst->rCopyMapHead;
			gDst->rCopyMapHead = tmpCopyMap;
		    }
		    //gSrc->sCopyMapHead->map.reserve(mfSrc.size());
		    //gDst->rCopyMapHead->map.reserve(mfDst.size());
		    gSrc->sCopyMapHead->alloc_CopyMap(mfSrc);
		    gDst->rCopyMapHead->alloc_CopyMap(mfDst);
		}

		//if(gSrc->numTasks != mfSrc.IndexArray().size())
		//    std::cout<< "before " <<gSrc->numTasks << " now " <<mfSrc.size() << " at gID " << gSrc->graphID << std::endl;	

		gSrc->numFabs = mfSrc.size();
		gDst->numFabs = mfDst.size();	
		gSrc->numTasks = mfSrc.IndexArray().size();
		gDst->numTasks = mfDst.IndexArray().size();
	    }
#ifdef USE_PERILLA_PTHREADS
//	    perilla::syncAllThreads();
#endif
	    const FabArrayBase::CPC *TheCPC= &mfDst.getCPC(ng, mfSrc, ngSrc, period);;

	    int nfabsSrc = mfSrc.IndexArray().size();
	    int nfabsDst = mfDst.IndexArray().size();

	    const int nloc_cpAsc = TheCPC->m_LocTags->size();
	    const int nsnds_cpAsc = TheCPC->m_SndTags->size();
	    const int nrcvs_cpAsc = TheCPC->m_RcvTags->size();     
#ifdef USE_PERILLA_PTHREADS
//	    perilla::syncAllThreads();
#endif

	    if(np > 1){
#ifdef USE_PERILLA_PTHREADS
//		if(perilla::isMasterThread())
#endif
		{
		    send_cctc.reserve(nsnds_cpAsc);

		    for (FabArrayBase::MapOfCopyComTagContainers::const_iterator m_it = TheCPC->m_SndTags->begin(),
			    m_End = TheCPC->m_SndTags->end();
			    m_it != m_End;
			    ++m_it)
		    {
			if(m_it->first != myProc)      // Not destined to me.
			{
			    send_pr.push_back(m_it->first);
			    send_cctc.push_back(&(m_it->second));
			}
		    }

		    recv_cctc.reserve(nrcvs_cpAsc);

		    for (FabArrayBase::MapOfCopyComTagContainers::const_iterator m_it = TheCPC->m_RcvTags->begin(),
			    m_End = TheCPC->m_RcvTags->end();
			    m_it != m_End;
			    ++m_it)
		    {
			if(m_it->first != myProc)      // I am not the source for this receipt
			{
			    recv_pr.push_back(m_it->first);
			    recv_cctc.push_back(&(m_it->second));
			}
		    }
		}
	    }
#ifdef USE_PERILLA_PTHREADS
//	    perilla::syncAllThreads();
#endif

//#ifndef USE_PERILLA_PTHREADS
	    #pragma omp parallel shared(gSrc, gDst, mfSrc, mfDst, nfabsSrc, nfabsDst)
//#endif
	    {
#ifdef _OPENMP
		int tid = omp_get_thread_num();//perilla::tid();//omp_get_thread_num();	  
#else
		int tid=0;
#endif
		int tg = tid/perilla::NUM_THREADS_PER_TEAM;//perilla::wid();//WorkerThread::perilla_wid();
 		int nt= tid%perilla::NUM_THREADS_PER_TEAM;
		int fg;
		//std::cout<<"thread "<< tid<<"group "<<tg<< "Before parallel at gID " << gDst->graphID << " numTask " << gDst->numTasks << " numFabs " << gDst->numFabs <<std::endl;	

		for(int f=0; f<nfabsSrc; f++)
		{
		  //if(perilla::isMasterWorkerThread())
		  if(nt==0)
		    if(WorkerThread::isMyRegion(tg,f))
		    {
			int scnt = 0;
			FabCopyAssoc *cpSrc;
			//if(gDst->graphID > 25)
			//std::cout<< "Inside parallel Generating Send at tid " << tid << " f " << f << " gID " << gDst->graphID <<std::endl;	

			if(gSrc->task[f]->cpAsc_srcHead == 0)
			{
			    gSrc->task[f]->cpAsc_srcHead = new FabCopyAssoc();
			    cpSrc = gSrc->task[f]->cpAsc_srcHead;
			}
			else
			{
			    cpSrc = new FabCopyAssoc();
			    cpSrc->next = gSrc->task[f]->cpAsc_srcHead;
			    gSrc->task[f]->cpAsc_srcHead = cpSrc;
			}

			cpSrc->graphPartner = gDst;
			cpSrc->l_con.nscpy = 0;
			for(int i=0; i<nloc_cpAsc; i++)
			{
			    const FabArrayBase::CopyComTag& tag = (*TheCPC->m_LocTags)[i];
			    //if(f == tag.srcIndex)
			    if(mfSrc.IndexArray()[f] == tag.srcIndex)
				cpSrc->l_con.nscpy++;		  
			}
			cpSrc->l_con.scpy = new LocalCopyDescriptor[cpSrc->l_con.nscpy];		

			//if(gDst->graphID == 4 && tag.dstIndex == 60 )
			//std::cout<< "Inside parallel Generating Local Copy send at tid " << tid << " f " << f << " gID " << gDst->graphID <<" num local connections"<< nloc_cpAsc << std::endl;	

			for(int i=0; i<nloc_cpAsc; i++)
			{
			    const FabArrayBase::CopyComTag *tag = &(*TheCPC->m_LocTags)[i];
			    //if(f == tag.srcIndex)
			    if(mfSrc.IndexArray()[f] == tag->srcIndex)			
			    {
				cpSrc->l_con.scpy[scnt].ns = mfSrc.localindex(tag->srcIndex);
				cpSrc->l_con.scpy[scnt].nd = mfDst.localindex(tag->dstIndex);
				cpSrc->l_con.scpy[scnt].sbx = tag->sbox;
				cpSrc->l_con.scpy[scnt].dbx = tag->dbox;
				int psize = tag->sbox.numPts() * mfSrc.nComp(); //---------------------------------------------------------------????????????????
				//std::cout<< " gSrc ID "<< gSrc->graphID << " f "<<f<< " sndPkgsize " << psize <<std::endl;
				for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
				{
				    Package *tmpPkg = new Package(psize);
				    for(int j=0; j<psize; j++)
					tmpPkg->databuf[j] = 0;
				    cpSrc->l_con.scpy[scnt].pQueue.enqueue(tmpPkg);
				}
				for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
				    cpSrc->l_con.scpy[scnt].recycleQueue.enqueue(cpSrc->l_con.scpy[scnt].pQueue.dequeue());
				scnt++;
			    }
			}

			if(np > 1)
			{
			    cpSrc->r_con.nsnd = 0;
			    cpSrc->r_con.remotePushReady = false;
			    cpSrc->r_con.firingRuleCnt = 0;
			    for(int i=0; i<nsnds_cpAsc; i++)
			    {
				const FabArrayBase::CopyComTagsContainer& cctc = *send_cctc[i];
				for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
					it != cctc.end(); ++it)
				{
				    if(mfSrc.IndexArray()[f] == it->srcIndex)				
					cpSrc->r_con.nsnd++;
				}		      
			    } // for(i<nsnds_cpAsc)
			    cpSrc->r_con.snd = new RemoteCommDescriptor[cpSrc->r_con.nsnd];
			    scnt = 0;
			    for(int i=0; i<nsnds_cpAsc; i++)
			    {
				const FabArrayBase::CopyComTagsContainer& cctc = *send_cctc[i];
				for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
					it != cctc.end(); ++it)
				{
				    if(mfSrc.IndexArray()[f] == it->srcIndex)
				    {
					cpSrc->r_con.snd[scnt].ns = it->srcIndex;
					cpSrc->r_con.snd[scnt].nd = it->dstIndex;			      
					cpSrc->r_con.snd[scnt].lns = mfSrc.localindex(it->srcIndex);
					cpSrc->r_con.snd[scnt].lnd = mfDst.localindex(it->dstIndex);			      
					cpSrc->r_con.snd[scnt].sbx = it->sbox;
					cpSrc->r_con.snd[scnt].dbx = it->dbox;
					int psize = it->sbox.numPts() * mfSrc.nComp(); //---------------------------------------------------------------????????????????

					for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
					{
					    Package *tmpPkg = new Package(psize);
					    for(int j=0; j<psize; j++)
						tmpPkg->databuf[j] = 0;
					    cpSrc->r_con.snd[scnt].pQueue.enqueue(tmpPkg);
					}
					for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
					    cpSrc->r_con.snd[scnt].recycleQueue.enqueue(cpSrc->r_con.snd[scnt].pQueue.dequeue());			      
					scnt++;
				    }
				}		      		      
			    } // for(i<nsnds_cpAsc)		    		  
			} // if(np > 1)	      																	
		    } // if(fg==tg)

		    //perilla::syncAllThreads();
		    #pragma omp barrier
		    if(np > 1)
		    {
			//if(WorkerThread::perilla_isMasterWorkerThread() && tg==0)
		        if(tid==0)
			{

			    // std::cout<< "Inside parallel Generating Remote Send tg 0 at tid " << tid << " f " << f << " gID " << gDst->graphID <<std::endl;	

			    gSrc->sCopyMapHead->map[f]->r_con.nsnd = 0;
			    gSrc->sCopyMapHead->map[f]->r_con.firingRuleCnt = 0;
			    for(int i=0; i<nsnds_cpAsc; i++)
			    {
				const FabArrayBase::CopyComTagsContainer& cctc = *send_cctc[i];
				for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
					it != cctc.end(); ++it)
				{
				    if(mfSrc.IndexArray()[f] == it->srcIndex)				
					gSrc->sCopyMapHead->map[f]->r_con.nsnd++;
				}
			    } // for(i<nsnds_cpAsc)
			    gSrc->sCopyMapHead->map[f]->r_con.snd = new RemoteCommDescriptor[gSrc->sCopyMapHead->map[f]->r_con.nsnd];
			    int scnt = 0;
			    for(int i=0; i<nsnds_cpAsc; i++)
			    {
				const FabArrayBase::CopyComTagsContainer& cctc = *send_cctc[i];
				for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
					it != cctc.end(); ++it)
				{

				    if(mfSrc.IndexArray()[f] == it->srcIndex)
				    {

					//if(gDst->graphID == 31 && (it->dstIndex == 519))
					//std::cout <<"myP " <<myProc<< " Added in S Dep nd " << it->dstIndex << " ns "<< it->srcIndex << " f " << f << " i "<< scnt << " tg " <<tg << std::endl;

					gSrc->sCopyMapHead->map[f]->r_con.snd[scnt].ns = it->srcIndex;
					gSrc->sCopyMapHead->map[f]->r_con.snd[scnt].nd = it->dstIndex;
					gSrc->sCopyMapHead->map[f]->r_con.snd[scnt].r_gid = gDst->graphID-1;
					gSrc->sCopyMapHead->map[f]->r_con.snd[scnt].r_grids = (gDst->numFabs > gSrc->numFabs ? gDst->numFabs : gSrc->numFabs);
					gSrc->sCopyMapHead->map[f]->r_con.snd[scnt].lns = mfSrc.localindex(it->srcIndex);
					gSrc->sCopyMapHead->map[f]->r_con.snd[scnt].lnd = mfDst.localindex(it->dstIndex);				  
					gSrc->sCopyMapHead->map[f]->r_con.snd[scnt].sbx = it->sbox;
					gSrc->sCopyMapHead->map[f]->r_con.snd[scnt].dbx = it->dbox;

					int psize = it->sbox.numPts() * mfSrc.nComp(); //---------------------------------------------------------------????????????????

					gSrc->sCopyMapHead->map[f]->r_con.snd[scnt].sz = psize;
					gSrc->sCopyMapHead->map[f]->r_con.snd[scnt].pr = send_pr[i];

					for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
					{
					    //Package *tmpPkg = new Package(psize);
					    Package *tmpPkg = new Package();
					    //for(int j=0; j<psize; j++)
					    //tmpPkg->databuf[j] = 0;
					    gSrc->sCopyMapHead->map[f]->r_con.snd[scnt].pQueue.enqueue(tmpPkg);
					}
					for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
					    gSrc->sCopyMapHead->map[f]->r_con.snd[scnt].recycleQueue.enqueue(gSrc->sCopyMapHead->map[f]->r_con.snd[scnt].pQueue.dequeue());  
					scnt++;
				    }
				}
			    } // for(i<nsnds_cpAsc)
			} // if(tid==0)
		    } // if(np > 1)	  
		} // for(f<nfabsSrc)
		//	  std::cout<< "Barrier 2 " <<" tid "<<tid<<std::endl;	  
		//perilla::syncAllThreads();
		#pragma omp barrier
		for(int f=0; f<nfabsDst; f++)
		{
		  //if(perilla::isMasterWorkerThread())
		  if(nt==0)
		    if(WorkerThread::isMyRegion(tg,f))		
		    {
			//	  std::cout <<"tid: "<< tid << " f: "<< f << " is master "<<WorkerThread::isTeamMasterThread(tid) << " is my region "<<WorkerThread::isMyRegion(tg,f)<<std::endl;		  

			//if(gDst->graphID > 25)
			//std::cout<< "Inside parallel Generating Recive at tid " << tid << " f " << f << " gID " << gDst->graphID <<std::endl;	

			FabCopyAssoc *cpDst;
			if(gDst->task[f]->cpAsc_dstHead == 0)
			{
			    gDst->task[f]->cpAsc_dstHead = new FabCopyAssoc();
			    cpDst = gDst->task[f]->cpAsc_dstHead;
			}
			else
			{
			    cpDst = new FabCopyAssoc();
			    cpDst->next = gDst->task[f]->cpAsc_dstHead;
			    gDst->task[f]->cpAsc_dstHead = cpDst;
			}
			cpDst->graphPartner = gSrc;
			cpDst->l_con.ndcpy = 0;
			cpDst->l_con.firingRuleCnt = 0;
			cpDst->l_con.dcpyCnt = 0;
			for(int i=0; i<nloc_cpAsc; i++)
			{
			    const FabArrayBase::CopyComTag *tag = &(*TheCPC->m_LocTags)[i];
			    //if(f == tag.dstIndex)
			    if(mfDst.IndexArray()[f] == tag->dstIndex)
				cpDst->l_con.ndcpy++;		  
			}
			cpDst->l_con.dcpy = new LocalCopyDescriptor[cpDst->l_con.ndcpy];		
			int dcnt = 0;

			//if(gDst->graphID > 25)
			//std::cout<< "Inside parallel Generating Local copy recive at tid " << tid << " f " << f << " gID " << gDst->graphID <<std::endl;	


			//if(gDst->graphID ==27 && f == 633)
			//std::cout<< "tid " << tid << " f " << f << " gID " << gDst->graphID << " numReciv " << nloc_cpAsc << " ndcpy " << cpDst->l_con.ndcpy <<std::endl;	


			for(int i=0; i<nloc_cpAsc; i++)
			{
			    const FabArrayBase::CopyComTag *tag = &(*TheCPC->m_LocTags)[i];
			    //if(f == tag->dstIndex)
			    if(mfDst.IndexArray()[f] == tag->dstIndex)
			    {

				//if(gDst->graphID == 4 && (tag->dstIndex == 60 || tag->dstIndex == 59))
				//std::cout<< "dcpy tid " << tid << " f " << f << " i " << i << " dcnt " << dcnt << " ns "<<tag->srcIndex << " nd "<<tag->dstIndex << " lo " << tag->dbox.smallEnd() << " hi " << tag->dbox.bigEnd() <<std::endl;	

				cpDst->l_con.dcpy[dcnt].ns = mfSrc.localindex(tag->srcIndex);
				cpDst->l_con.dcpy[dcnt].nd = mfDst.localindex(tag->dstIndex);
				cpDst->l_con.dcpy[dcnt].sbx = tag->sbox;
				cpDst->l_con.dcpy[dcnt].dbx = tag->dbox;

				// if(gDst->graphID > 25 && f == 633)
				//std::cout<< " Generating Package tid " << tid << " i " << i <<std::endl;	

				int psize = tag->dbox.numPts() * mfSrc.nComp(); //---------------------------------------------------------------????????????????
				cpDst->l_con.dcpy[dcnt].sz = psize;

				if(!gDst->isDepGraph)
				{
				    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
				    {
					Package *tmpPkg = new  Package(psize);

					// if(tmpPkg == nullptr)
					//std::cout<<"Found the culprit tid " << tid << " f " << f << " i " << i << std::endl;

					for(int j=0; j<psize; j++)
					    tmpPkg->databuf[j] = 0;
					cpDst->l_con.dcpy[dcnt].pQueue.enqueue(tmpPkg);
				    }

				    // if(gDst->graphID > 25 && f == 633)
				    //std::cout<< " Generating  now in reQ Package tid " << tid << " i " << i <<std::endl;	

				    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
					cpDst->l_con.dcpy[dcnt].recycleQueue.enqueue(cpDst->l_con.dcpy[dcnt].pQueue.dequeue());		      

				    //if(gDst->graphID > 25 && f == 633)
				    //  std::cout<< " Generated Package tid " << tid << " i " << i <<std::endl;	
				}

				dcnt++;
			    }
			}

			// if(gDst->graphID > 25 && f > 630)
			//std::cout<< "Safe now tid " << tid << " f " << f << " gID " << gDst->graphID << " numReciv " << nloc_cpAsc <<std::endl;	

			RegionGraph* depGraph = gDst->srcLinkGraph;
			for(int df=0; df < gDst->task[f]->depTaskIDs.size(); df++)
			{
			    int dfi = gDst->task[f]->depTaskIDs[df];
			    FabCopyAssoc *cpdDst = depGraph->task[dfi]->cpAsc_dstHead;
			    for(int i=0; i<cpdDst->l_con.ndcpy ; i++)
			    {
				for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
				{
				    int psize = cpdDst->l_con.dcpy[i].sz;
				    Package *tmpPkg = new  Package(psize);       
				    for(int j=0; j<psize; j++)
					tmpPkg->databuf[j] = 0;
				    cpdDst->l_con.dcpy[i].pQueue.enqueue(tmpPkg);
				}			      
				for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
				    cpdDst->l_con.dcpy[i].recycleQueue.enqueue(cpdDst->l_con.dcpy[i].pQueue.dequeue());		      
			    }
			}		  

			if(np > 1)
			{
			    cpDst->r_con.nrcv = 0;
			    cpDst->r_con.remotePullDone = false;
			    cpDst->r_con.firingRuleCnt = 0;
			    for(int i=0; i<nrcvs_cpAsc; i++)
			    {
				const FabArrayBase::CopyComTagsContainer& cctc = *recv_cctc[i];
				for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
					it != cctc.end(); ++it)
				{
				    if(mfDst.IndexArray()[f] == it->dstIndex)
					cpDst->r_con.nrcv++;
				}		      
			    } // for(i<nrcvs_cpAsc)
			    cpDst->r_con.rcv = new RemoteCommDescriptor[cpDst->r_con.nrcv];
			    dcnt = 0;
			    for(int i=0; i<nrcvs_cpAsc; i++)
			    {
				const FabArrayBase::CopyComTagsContainer& cctc = *recv_cctc[i];
				for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
					it != cctc.end(); ++it)
				{
				    //if(f == it->dstIndex)
				    if(mfDst.IndexArray()[f] == it->dstIndex)
				    {
					cpDst->r_con.rcv[dcnt].nd = it->dstIndex;
					cpDst->r_con.rcv[dcnt].ns = it->srcIndex;
					cpDst->r_con.rcv[dcnt].lnd = mfDst.localindex(it->dstIndex);
					cpDst->r_con.rcv[dcnt].lns = mfSrc.localindex(it->srcIndex);				  
					cpDst->r_con.rcv[dcnt].sbx = it->sbox;
					cpDst->r_con.rcv[dcnt].dbx = it->dbox;
					int psize = it->dbox.numPts() * mfDst.nComp(); //---------------------------------------------------------------????????????????
					cpDst->r_con.rcv[dcnt].sz = psize;

					if(!gDst->isDepGraph)
					{
					    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
					    {
						Package *tmpPkg = new Package(psize);
						for(int j=0; j<psize; j++)
						    tmpPkg->databuf[j] = 0;
						cpDst->r_con.rcv[dcnt].pQueue.enqueue(tmpPkg);
					    }
					    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
						cpDst->r_con.rcv[dcnt].recycleQueue.enqueue(cpDst->r_con.rcv[dcnt].pQueue.dequeue());		      			      
					}

					dcnt++;
				    }
				}
			    }// for(i<nrcvs_cpAsc)

			    RegionGraph* depGraph = gDst->srcLinkGraph;
			    for(int df=0; df < gDst->task[f]->depTaskIDs.size(); df++)
			    {
				int dfi = gDst->task[f]->depTaskIDs[df];
				FabCopyAssoc *cpdDst = depGraph->task[dfi]->cpAsc_dstHead;
				for(int i=0; i<cpdDst->r_con.nrcv ; i++)
				{
				    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
				    {
					int psize = cpdDst->r_con.rcv[i].sz;
					Package *tmpPkg = new  Package(psize);       
					for(int j=0; j<psize; j++)
					    tmpPkg->databuf[j] = 0;
					cpdDst->r_con.rcv[i].pQueue.enqueue(tmpPkg);
				    }			      
				    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
					cpdDst->r_con.rcv[i].recycleQueue.enqueue(cpdDst->r_con.rcv[i].pQueue.dequeue());		      
				}
			    }		  


			} // if(np > 1)
		    }// if(fg==tg)

		    //perilla::syncAllThreads();
		    #pragma omp barrier

		    if(np > 1)
		    {
			//if(WorkerThread::perilla_isMasterWorkerThread() && tg==0)			
			if(tid==0)
			{
			    //  std::cout<< "Inside parallel Generating Remote Recive tg 0 at tid " << tid << " f " << f << " gID " << gDst->graphID <<std::endl;	
			    gDst->rCopyMapHead->map[f]->r_con.nrcv = 0;
			    gDst->rCopyMapHead->map[f]->r_con.firingRuleCnt = 0;
			    for(int i=0; i<nrcvs_cpAsc; i++)
			    {
				const FabArrayBase::CopyComTagsContainer& cctc = *recv_cctc[i];
				for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
					it != cctc.end(); ++it)
				{
				    //if(f == it->dstIndex)
				    if(mfDst.IndexArray()[f] == it->dstIndex)
					gDst->rCopyMapHead->map[f]->r_con.nrcv++;
				}
			    }
			    gDst->rCopyMapHead->map[f]->r_con.rcv = new RemoteCommDescriptor[gDst->rCopyMapHead->map[f]->r_con.nrcv];
			    int dcnt = 0;
			    for(int i=0; i<nrcvs_cpAsc; i++)
			    {
				const FabArrayBase::CopyComTagsContainer& cctc = *recv_cctc[i];
				for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
					it != cctc.end(); ++it)
				{
				    //if(f == it->dstIndex)
				    if(mfDst.IndexArray()[f] == it->dstIndex)
				    {

					// if(myProc==54 && gDst->graphID == 25 && f == 10)
					// std::cout <<"myP " <<myProc<<" Dep n R Added nd " << it->dstIndex << " ns "<< it->srcIndex << " f " << f << " sgID "<< gSrc->graphID <<" tg "<<tg<< " from P " << recv_pr[i] <<std::endl;

					gDst->rCopyMapHead->map[f]->r_con.rcv[dcnt].nd = it->dstIndex;
					gDst->rCopyMapHead->map[f]->r_con.rcv[dcnt].ns = it->srcIndex;
					gDst->rCopyMapHead->map[f]->r_con.rcv[dcnt].lnd = mfDst.localindex(it->dstIndex);
					gDst->rCopyMapHead->map[f]->r_con.rcv[dcnt].lns = mfSrc.localindex(it->srcIndex);
					gDst->rCopyMapHead->map[f]->r_con.rcv[dcnt].r_grids = (gDst->numFabs > gSrc->numFabs ? gDst->numFabs : gSrc->numFabs);
					gDst->rCopyMapHead->map[f]->r_con.rcv[dcnt].sbx = it->sbox;
					gDst->rCopyMapHead->map[f]->r_con.rcv[dcnt].dbx = it->dbox;

					int psize = it->dbox.numPts() * mfDst.nComp(); //---------------------------------------------------------------????????????????

					gDst->rCopyMapHead->map[f]->r_con.rcv[dcnt].sz = psize;
					gDst->rCopyMapHead->map[f]->r_con.rcv[dcnt].pr = recv_pr[i];

					BL_ASSERT(gDst->rCopyMapHead->map[f]->r_con.rcv[dcnt].lnd == f);

					if(Perilla::genTags)
					{
					    try{
						std::map<int,int>::iterator itr = tagMap[recv_pr[i]][gDst->graphID-1][it->dstIndex][it->srcIndex].find(psize);
						if( itr != tagMap[recv_pr[i]][gDst->graphID-1][it->dstIndex][it->srcIndex].end())
						{
						    //gDst->rCopyMapHead->map[f]->r_con.rcv[dcnt].lnd = itr->second;
						}
						else
						{
						    tagMap[recv_pr[i]][gDst->graphID-1][it->dstIndex][it->srcIndex][psize] = Perilla::uTags++;
						    //gDst->rCopyMapHead->map[f]->r_con.rcv[dcnt].lnd = Perilla::uTags++;
						    std::map<int,int>::iterator itr2 = pTagCnt[recv_pr[i]].find(gDst->graphID-1);
						    if(itr2 != pTagCnt[recv_pr[i]].end())
							pTagCnt[recv_pr[i]][gDst->graphID-1] = pTagCnt[recv_pr[i]][gDst->graphID-1] + 1;
						    else
							pTagCnt[recv_pr[i]][gDst->graphID-1] = 1;									     									      
						}
					    }
					    catch(std::exception& e)
					    {
						std::cout <<"Inside tagGeneration gID "<< gDst->graphID <<" "<< e.what() << '\n';
					    }
					}
					//tagMap[recv_pr[i]][gDst->graphID][it->dstIndex][it->srcIndex] = pTagCnt[recv_pr[i]];				  


					for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
					{
					    //Package *tmpPkg = new Package(psize);
					    Package *tmpPkg = new Package();
					    //for(int j=0; j<psize; j++)
					    //tmpPkg->databuf[j] = 0;
					    gDst->rCopyMapHead->map[f]->r_con.rcv[dcnt].pQueue.enqueue(tmpPkg);
					}
					for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
					    gDst->rCopyMapHead->map[f]->r_con.rcv[dcnt].recycleQueue.enqueue(gDst->rCopyMapHead->map[f]->r_con.rcv[dcnt].pQueue.dequeue());
					dcnt++;
				    }
				}
			    } // for(i<nrcvs_cpAsc)

			} // if(tid==0)
		    } // if(np > 1)
		    //} //if(fg==tg)
	    } // for(f<nfabsDst)
	    // std::cout<< "Barrier 4" <<" tid "<<tid <<std::endl;	      	  	  
	    //perilla::syncAllThreads();
	    #pragma omp for
	    for(int f=0; f<nfabsSrc; f++)
	    {
	      //if(perilla::isMasterWorkerThread())
              if(nt==0)
		if(WorkerThread::isMyRegion(tg,f))	      
		{	

		    //if(gDst->graphID > 25)
		    //std::cout<< "Inside parallel Generating Send partners at tid " << tid << " f " << f << " gID " << gDst->graphID <<std::endl;		  

		    for(int i=0; i<gSrc->task[f]->cpAsc_srcHead->l_con.nscpy; i++)
		    {
			int lnd = gSrc->task[f]->cpAsc_srcHead->l_con.scpy[i].nd;
			for(int j=0; j<gDst->task[ lnd ]->cpAsc_dstHead->l_con.ndcpy; j++)
			    if(gSrc->task[f]->cpAsc_srcHead->l_con.scpy[i].dbx == gDst->task[ lnd ]->cpAsc_dstHead->l_con.dcpy[j].dbx)
				gSrc->task[f]->cpAsc_srcHead->l_con.scpy[i].dPartner = j;
		    }
		}
	    } // for(f<nfabsSrc)
	    //std::cout<< "Barrier 5" <<" tid "<<tid<<std::endl;	      	  	  
	    //perilla::syncAllThreads();
	    #pragma omp for
	    for(int f=0; f<nfabsDst; f++)
	    {
	      //if(perilla::isMasterWorkerThread())
	      if(nt==0)
		if(WorkerThread::isMyRegion(tg,f))
		{
		    //if(gDst->graphID > 25)
		    //std::cout<< "Inside parallel Generating Recive partners at tid " << tid << " f " << f << " gID " << gDst->graphID <<std::endl;	

		    for(int i=0; i<gDst->task[f]->cpAsc_dstHead->l_con.ndcpy; i++)
		    {
			int lns = gDst->task[f]->cpAsc_dstHead->l_con.dcpy[i].ns;
			for(int j=0; j<gSrc->task[ lns ]->cpAsc_srcHead->l_con.nscpy; j++)
			    if(gDst->task[f]->cpAsc_dstHead->l_con.dcpy[i].dbx == gSrc->task[ lns ]->cpAsc_srcHead->l_con.scpy[j].dbx)
				gDst->task[f]->cpAsc_dstHead->l_con.dcpy[i].sPartner = j;
		    }
		}
	    } // for(f<nfabsDst)																 
	} // omp parallel
    } // if(!(*mfSrc == *mfDst))    
}
catch(std::exception& e)
{
    std::cout <<"Inside MFcopyAssoc gID "<< gDst->graphID <<" "<< e.what() << '\n';
}


//std::cout<< "All done safely at gID " << gDst->graphID <<std::endl;	

#endif

} // multifabExtractCopyAssoc

#endif



#if 0
void Perilla::multifabCopyPull(RegionGraph* destGraph, RegionGraph* srcGraph, MultiFab* mfDst, MultiFab* mfSrc, int f, int dstcomp, int srccomp, int nc, int ng, int ngsrc, bool singleT)
{
    int myProc = ParallelDescriptor::MyProc();

    int ntid = WorkerThread::perilla_wtid();
    int tg = WorkerThread::perilla_wid();
    //MultiFab* mfDst = destGraph->assocMF;
    //MultiFab* mfSrc = srcGraph->assocMF;
    if(nc<1) cout <<"MULTIFAB_COPY_C: nc must be >= 1"<< endl;
    if(mfDst->nComp() < (dstcomp-1)) cout <<"MULTIFAB_COPY_C: nc too large for dst multifab"<< endl;
    //if(mfSrc->nComp() < (srccomp-1)) cout <<"MULTIFAB_COPY_C: nc too large for src multifab"<< endl;

    if(true)//if(!(*mfDst == *mfSrc))
    {
	if(ng > mfDst->nGrow()) cout <<"MULTIFAB_COPY_C: ng > 0 not supported in parallel copy"<< endl;
	//if(ngsrc > mfSrc->nGrow()) cout <<"MULTIFAB_COPY_C: ngsrc > msrc%ng"<< endl;
	FabCopyAssoc* cpDst = destGraph->task[f]->cpAsc_dstHead;
	while(cpDst != 0)
	{
	    if(cpDst->graphPartner == srcGraph)
		break;
	    cpDst = cpDst->next;
	} 
	if(cpDst == 0) cout <<"Metadata for across grid copy not found"<< endl;
	//destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);

	if(singleT)
	{
	    pthread_mutex_lock(&(cpDst->l_con.dLock));
	    for(int i=0; i<cpDst->l_con.ndcpy; i++)
	    {
		Package* rcvPackage = cpDst->l_con.dcpy[i].pQueue.getFront(true); // corrected from recycleQ to pQ
		mfDst->m_fabs_v[f]->copyFromMem(cpDst->l_con.dcpy[i].dbx,dstcomp,nc,rcvPackage->databuf);
	    }	
	    for(int i=0; i<cpDst->l_con.ndcpy; i++)
		cpDst->l_con.dcpy[i].recycleQueue.enqueue(cpDst->l_con.dcpy[i].pQueue.dequeue(true),true); // corrected from pQ to recycleQ and from recycleQ to pQ
	    cpDst->l_con.firingRuleCnt = cpDst->l_con.firingRuleCnt - cpDst->l_con.ndcpy;
	    pthread_mutex_unlock(&(cpDst->l_con.dLock));
	}
	else
	{
	    if(ntid==0)
		pthread_mutex_lock(&(cpDst->l_con.dLock));
	    destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);

	    for(int i=0; i<cpDst->l_con.ndcpy; i++)
		if((i%(perilla::NUM_THREADS_PER_TEAM-1)) == ntid)
		{
		    Package* rcvPackage = cpDst->l_con.dcpy[i].pQueue.getFront(true); // corrected from recycleQ to pQ
		    mfDst->m_fabs_v[f]->copyFromMem(cpDst->l_con.dcpy[i].dbx,dstcomp,nc,rcvPackage->databuf);
		}	
	    destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);       	

	    if(ntid == 0)
	    {
		for(int i=0; i<cpDst->l_con.ndcpy; i++)
		    cpDst->l_con.dcpy[i].recycleQueue.enqueue(cpDst->l_con.dcpy[i].pQueue.dequeue(true),true); // corrected from pQ to recycleQ and from recycleQ to pQ
		cpDst->l_con.firingRuleCnt = cpDst->l_con.firingRuleCnt - cpDst->l_con.ndcpy;
		pthread_mutex_unlock(&(cpDst->l_con.dLock));
	    }
	    destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
	}

	int np = ParallelDescriptor::NProcs();
	if(np == 1)
	    return;

	if(singleT)
	{
	    pthread_mutex_lock(&(destGraph->rCopyMapHead->map[f]->r_con.rcvLock));
	    pthread_mutex_lock(&(cpDst->r_con.rcvLock));
	    for(int i=0; i<cpDst->r_con.nrcv; i++)
	    {
		///*
		Package *rcvMetaPackage = destGraph->rCopyMapHead->map[f]->r_con.rcv[i].pQueue.dequeue(true);
		rcvMetaPackage->completed = false;
		rcvMetaPackage->served = false;
		rcvMetaPackage->request = MPI_REQUEST_NULL;	  
		destGraph->rCopyMapHead->map[f]->r_con.rcv[i].recycleQueue.enqueue(rcvMetaPackage,true);

		Package* rcvPackage = cpDst->r_con.rcv[i].pQueue.dequeue(true);                               // corrected from recycleQ to pQ
		mfDst->m_fabs_v[f]->copyFromMem(cpDst->r_con.rcv[i].dbx,dstcomp,nc,rcvPackage->databuf);	      
		rcvPackage->notified = false;
		rcvPackage->completed = false;
		cpDst->r_con.rcv[i].recycleQueue.enqueue(rcvPackage,true);                         // corrected from pQ to recycleQ	      
		//*/

		//Package* rcvPackage = cpDst->r_con.rcv[i].pQueue.getFront(true);                               // corrected from recycleQ to pQ
		//mfDst->m_fabs_v[f]->copyFromMem(cpDst->r_con.rcv[i].dbx,dstcomp,nc,rcvPackage->databuf);
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
	    pthread_mutex_unlock(&(destGraph->rCopyMapHead->map[f]->r_con.rcvLock));

	}
	else
	{	
	    if(ntid==0)
	    {
		pthread_mutex_lock(&(destGraph->rCopyMapHead->map[f]->r_con.rcvLock));
		pthread_mutex_lock(&(cpDst->r_con.rcvLock));
	    }
	    destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);

	    for(int i=0; i<cpDst->r_con.nrcv; i++)
		if((i%(perilla::NUM_THREADS_PER_TEAM-1)) == ntid)
		{
		    ///*
		    Package *rcvMetaPackage = destGraph->rCopyMapHead->map[f]->r_con.rcv[i].pQueue.dequeue(true);
		    rcvMetaPackage->completed = false;
		    rcvMetaPackage->served = false;
		    rcvMetaPackage->request = MPI_REQUEST_NULL;	  
		    destGraph->rCopyMapHead->map[f]->r_con.rcv[i].recycleQueue.enqueue(rcvMetaPackage,true);

		    Package* rcvPackage = cpDst->r_con.rcv[i].pQueue.dequeue(true);                               // corrected from recycleQ to pQ
		    mfDst->m_fabs_v[f]->copyFromMem(cpDst->r_con.rcv[i].dbx,dstcomp,nc,rcvPackage->databuf);	      
		    rcvPackage->notified = false;
		    rcvPackage->completed = false;
		    cpDst->r_con.rcv[i].recycleQueue.enqueue(rcvPackage,true);                         // corrected from pQ to recycleQ	      
		    //*/

		    //Package* rcvPackage = cpDst->r_con.rcv[i].pQueue.getFront(true);                               // corrected from recycleQ to pQ
		    //mfDst->m_fabs_v[f]->copyFromMem(cpDst->r_con.rcv[i].dbx,dstcomp,nc,rcvPackage->databuf);

		}
	    destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);

	    if(ntid==0)
	    {
		cpDst->r_con.firingRuleCnt = cpDst->r_con.firingRuleCnt - cpDst->r_con.nrcv;

		cpDst->r_con.remotePullDone = true;
		///*
		for(int i=0; i<cpDst->r_con.nrcv; i++)
		    if(cpDst->r_con.rcv[i].pQueue.queueSize(true) >= 1)
			if(cpDst->r_con.rcv[i].pQueue.getFront(true)->checkRequest())
			    cpDst->r_con.firingRuleCnt++;
		//*/
		pthread_mutex_unlock(&(cpDst->r_con.rcvLock));
		pthread_mutex_unlock(&(destGraph->rCopyMapHead->map[f]->r_con.rcvLock));
	    }
	    destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
	}
    } // if(!(*mfDst == *mfSrc))

} // multifabCopyPull
#endif



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




  void Perilla::fillBoundaryPull(amrex::RGIter& rgi, RegionGraph* rg, amrex::MultiFab& mf)
  {
    if(rgi.currentItr != 1)
      return;

    int f = rgi.currentRegion;
    fillBoundaryPull(rg, &mf, f);
  }

  void Perilla::fillBoundaryPull(amrex::RGIter& rgi, amrex::MultiFab& mf)
  {
    if(rgi.currentItr != 1)
      return;

    int f = rgi.currentRegion;
    fillBoundaryPull(rgi.itrGraph, &mf, f);
  }


