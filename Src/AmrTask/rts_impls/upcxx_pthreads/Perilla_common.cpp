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


#if 0
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


#endif
