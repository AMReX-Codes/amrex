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
using namespace upcxx;


struct sMsgMap_t{
    std::map< int, std::map< int,  std::list< Package* > > > map; 
    volatile int size=0;
    pthread_mutex_t lock= PTHREAD_MUTEX_INITIALIZER;
}sMsgMap;

struct rMsgMap_t{
    std::map< int, std::map< int,  std::list< Package* > > > map; 
    volatile int size=0;
    pthread_mutex_t lock= PTHREAD_MUTEX_INITIALIZER;
}rMsgMap;

struct getReq_t{
    int src;
    int tag;
    upcxx::global_ptr<double> sbuf;
    int size;
    getReq_t(int _src, int _tag, upcxx::global_ptr<double> _sbuf, int _size):src(_src), tag(_tag), sbuf(_sbuf), size(_size){}
};

struct pendingGetList_t{
    std::list< getReq_t* > _pendingGets;
    pthread_mutex_t lock= PTHREAD_MUTEX_INITIALIZER;
    void add(getReq_t* req){
         pthread_mutex_lock(&lock);
	_pendingGets.push_back(req);
         pthread_mutex_unlock(&lock);
    }
    void process(){
        if(_pendingGets.size()==0) return;
        pthread_mutex_lock(&(rMsgMap.lock));
        pthread_mutex_lock(&lock);
        std::list< getReq_t* >::iterator it= _pendingGets.begin();
        while(it != _pendingGets.end()){
            double* localbuf= NULL;
            int src= (*it)->src; 
	    int tag= (*it)->tag;
            if(rMsgMap.map.find(src) != rMsgMap.map.end()){
	        if(rMsgMap.map[src].find(tag) != rMsgMap.map[src].end()){
                   if(rMsgMap.map[src][tag].size() >0){
                       rMsgMap.map[src][tag].front()->tag= tag;
                       localbuf= (rMsgMap.map[src][tag].front()->databuf).local(); //(double*) (static_cast<upcxx::global_ptr<void> > (rMsgMap.map[src][tag].front()->databuf).local());
                       *(rMsgMap.map[src][tag].front()->request)= upcxx::rget((*it)->sbuf, localbuf, (*it)->size);
                       rMsgMap.map[src][tag].pop_front();
                       rMsgMap.size--;
                       std::list< getReq_t* >::iterator it1= it;
                       it++;
		       delete (*it);
                       _pendingGets.erase(it1);
	           }else it++;
		}else it++;
            }else it++;
        }
        pthread_mutex_unlock(&lock);
        pthread_mutex_unlock(&(rMsgMap.lock));
    }
} pendingGetList;



#if 0

int Perilla::tagGen(int src, int dest, int channelID, int nFabs, int nChannels)
{
    int maxRange;

    maxRange = nFabs;
    if(nFabs*nChannels > perilla::MAX_SQRT_TAG) maxRange= 1024;    
    return (src%maxRange)*maxRange + (dest%maxRange) + channelID*(perilla::MAX_SQRT_TAG*(perilla::MAX_SQRT_TAG+1)/nChannels);

    //int nfabs = 256;

    //if(src >= nfabs || dest>=nfabs)
    // std::cout<<"Warnig Tag" << src << " " << dest << " "<<nFabs <<std::endl;

    //maxRange = nfabs*nfabs;
    //return src*nfabs + dest + channelID*(maxRange+2);


    /*int maxSR = *nChannels;

      std::cout << "tag " << (src%maxRange)*maxRange + (dest%maxRange) + channelID*(maxSR*(maxSR+1)/nChannels) << " " <<MPI_TAG_UB <<std::endl;

      if( (src%maxRange)*maxRange + (dest%maxRange) + channelID*(maxSR*(maxSR+1)/nChannels) >= MPI_TAG_UB )
      std::cout << "Out of Bound tag " << (src%maxRange)*maxRange + (dest%maxRange) + channelID*(maxSR*(maxSR+1)/nChannels) << " " <<MPI_TAG_UB <<std::endl;
      return (src%maxRange)*maxRange + (dest%maxRange) + channelID*(maxSR*(maxSR+1)/nChannels);
     */
}
#endif


void Perilla::syncProcesses(){
    upcxx::barrier();
}


#if 0
void Perilla::multifabBuildFabCon(RegionGraph* rg, const MultiFab& mf, const Periodicity& period)
{
    int np = ParallelDescriptor::NProcs();
    int myProc = ParallelDescriptor::MyProc();
    int numfabs = mf.IndexArray().size();
    bool cross = false;
    const FabArrayBase::FB& TheFB = mf.getFB(mf.nGrowVect(), period, false, false);
    const int n_loc_mf = TheFB.m_LocTags->size();
    const int n_snds_mf = TheFB.m_SndTags->size();
    const int n_rcvs_mf = TheFB.m_RcvTags->size();

    Vector<const FabArrayBase::CopyComTagsContainer*> send_cctc;
    Vector<int> send_pr;
    send_cctc.reserve(n_snds_mf);

    for (FabArrayBase::MapOfCopyComTagContainers::const_iterator m_it = TheFB.m_SndTags->begin(),
            m_End = TheFB.m_SndTags->end();
            m_it != m_End;
            ++m_it)
    {
        if(m_it->first != myProc)      // Not destined to me.
        {
            send_pr.push_back(m_it->first);
            send_cctc.push_back(&(m_it->second));
        }
    }

    Vector<const FabArrayBase::CopyComTagsContainer*> recv_cctc;
    Vector<int> recv_pr;
    recv_cctc.reserve(n_rcvs_mf);

    for (FabArrayBase::MapOfCopyComTagContainers::const_iterator m_it = TheFB.m_RcvTags->begin(),
            m_End = TheFB.m_RcvTags->end();
            m_it != m_End;
            ++m_it)
    {
        if(m_it->first != myProc)      // I am not the source for this receipt
        {
            recv_pr.push_back(m_it->first);
            recv_cctc.push_back(&(m_it->second));
        }
    }

//#pragma omp parallel shared(rg, mf, numfabs, np, TheFB, recv_cctc, send_cctc)
    {
        //int tg = omp_get_thread_num();
        int fg;
//        if(WorkerThread::perilla_isCommunicationThread())
//#pragma omp single
        {
            //bool cc = !mf->is_nodal(); //  cc = multifab_cell_centered_q(mf)
            //mf->sMap.reserve(numfabs);
            //mf->rMap.reserve(numfabs);
            //std::cout<< "Allocating sMap and rMap" <<std::endl;
            rg->alloc_lMap(mf);
            rg->alloc_sMap(mf);
            rg->alloc_rMap(mf);
        }
//#pragma omp barrier      
        //if(tid==0)                              
        {
            //bool cc = !mf->is_nodal(); //  cc = multifab_cell_centered_q(mf)
            //mf->sMap.reserve(numfabs);
            //mf->rMap.reserve(numfabs);
//#pragma omp for
            for(int f=0; f<numfabs; f++) //        !create local communication metadata for each fab
            {
//                if(WorkerThread::isMyRegion(tg,f) && WorkerThread::perilla_isMasterWorkerThread())
                {
                    rg->lMap[f]->l_con.nscpy = 0;

                    //for(int i=0; i<bxasc->l_con.ncpy; i++)
                    for(int i=0; i<n_loc_mf; i++)
                    {
                        const FabArrayBase::CopyComTag& tag = (*TheFB.m_LocTags)[i];

                        //std::cout << tag.srcIndex << " " << tag.dstIndex << " " <<tag.sbox.smallEnd() <<" "<< tag.sbox.bigEnd() << std::endl;

                        BL_ASSERT(mf.distributionMap[tag.dstIndex] == myProc);
                        BL_ASSERT(mf.distributionMap[tag.srcIndex] == myProc);
                        //get(tag.dstIndex).copy(get(tag.srcIndex),tag.box,scomp,tag.box,scomp,ncomp);
                        //if(f == local_index(mf,bxasc->l_con.cpy[i].ns)) //LocalIndex
                        if(mf.IndexArray()[f] == tag.srcIndex)
                            rg->lMap[f]->l_con.nscpy++;
                        //if(f == local_index(mf,bxasc->l_con.cpy[i].nd)) //LocalIndex
                        if(mf.IndexArray()[f] == tag.dstIndex)
                            rg->lMap[f]->l_con.ndcpy++;
                    }
                    /*
                       if(rg->lMap[f]->l_con.nscpy+rg->lMap[f]->l_con.ndcpy != n_loc_mf)
                       std::cout<< "Diff in Sum " << rg->lMap[f]->l_con.nscpy << " " <<rg->lMap[f]->l_con.ndcpy << " " << n_loc_mf <<std::endl;
                       BL_ASSERT(rg->lMap[f]->l_con.nscpy+rg->lMap[f]->l_con.ndcpy == n_loc_mf);
                     */
                }
            }
        }
//#pragma omp barrier
        //now we know how many copying segments each fab owns as source and destination allocate memory for metadata   
//#pragma omp for
        for(int f=0; f<numfabs; f++)
        {
            //fg = f % (omp_get_num_threads()/perilla::NUM_THREADS_PER_TEAM);   /// need to check if computing correct ???????
            //if((fg == tg) && ((tid%perilla::NUM_THREADS_PER_TEAM)==1))
//            if(WorkerThread::isMyRegion(tg,f) && WorkerThread::perilla_isMasterWorkerThread())
            {
                //omp_init_lock(&(rg->lMap[f]->l_con.sLock));
                //omp_init_lock(&(rg->lMap[f]->l_con.dLock));
                //omp_init_lock(&(rg->lMap[f]->l_con.ghostLock));

                //std::cout<< "MF l_con nscpy " <<rg->lMap[f]->l_con.nscpy << " ndcpy " << rg->lMap[f]->l_con.ndcpy <<std::endl;

                rg->lMap[f]->l_con.scpy = new LocalCopyDescriptor[rg->lMap[f]->l_con.nscpy];
                rg->lMap[f]->l_con.dcpy = new LocalCopyDescriptor[rg->lMap[f]->l_con.ndcpy];
                rg->lMap[f]->l_con.scpyCnt = 0;
                rg->lMap[f]->l_con.dcpyCnt = 0;
            }
        }
//#pragma omp barrier
        if(np > 1)
        {
//#pragma omp for
            for(int f=0; f<numfabs; f++)
            {
//                if(WorkerThread::perilla_isMasterWorkerThread() && WorkerThread::isMyRegion(tg,f))
                {
                    rg->lMap[f]->r_con.nrcv = 0;
                    rg->lMap[f]->r_con.nsnd = 0;
                    rg->lMap[f]->r_con.firingRuleCnt = 0;

                    //for(int i=0; i<bxasc->r_con.nsnd; i++)
                    for(int i=0; i<n_snds_mf; i++)
                    {
                        const FabArrayBase::CopyComTagsContainer& cctc = *send_cctc[i];
                        for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
                                it != cctc.end(); ++it)
                        {
                            //if(f == local_index(mf,bxasc->r_con.snd[i].ns)) //LocalIndex
                            if(mf.IndexArray()[f] == it->srcIndex)
                            {
                                rg->lMap[f]->r_con.nsnd++;
                            }
                        }
                    }
                    //for(int i=0; i<bxasc->r_con.nrcv; i++)
                    for(int i=0; i<n_rcvs_mf; i++)
                    {
                        const FabArrayBase::CopyComTagsContainer& cctc = *recv_cctc[i];
                        for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
                                it != cctc.end(); ++it)
                        {
                            //if(f == local_index(mf,bxasc->r_con.rcv[i].nd)) //LocalIndex
                            if(mf.IndexArray()[f] == it->dstIndex)
                            {
                                rg->lMap[f]->r_con.nrcv++;
                            }
                        }
                    }
                    //rg->sMap[f]->r_con.sndLock = new omp_lock_t;
                    //rg->rMap[f]->r_con.rcvLock = new omp_lock_t;
                    //omp_init_lock(rg->sMap[f]->r_con.sndLock);
                    //omp_init_lock(rg->rMap[f]->r_con.rcvLock);
                    rg->lMap[f]->r_con.snd = new RemoteCommDescriptor[rg->lMap[f]->r_con.nsnd];
                    rg->lMap[f]->r_con.rcv = new RemoteCommDescriptor[rg->lMap[f]->r_con.nrcv];
                }
            }
 //           if(WorkerThread::perilla_isMasterWorkerThread() && tg==0)
            {
//#pragma omp for
                for(int f=0; f<numfabs; f++)
                {
                    rg->rMap[f]->r_con.nrcv = 0;
                    rg->sMap[f]->r_con.nsnd = 0;

                    //for(int i=0; i<bxasc->r_con.nsnd; i++)
                    for(int i=0; i<n_snds_mf; i++)
                    {
                        const FabArrayBase::CopyComTagsContainer& cctc = *send_cctc[i];
                        for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
                                it != cctc.end(); ++it)
                        {
                            //if(f == local_index(mf,bxasc->r_con.snd[i].ns)) //LocalIndex
                            if(mf.IndexArray()[f] == it->srcIndex)
                            {
                                rg->sMap[f]->r_con.nsnd++;
                            }
                        }
                    }
                    //for(int i=0; i<bxasc->r_con.nrcv; i++)
                    for(int i=0; i<n_rcvs_mf; i++)
                    {
                        const FabArrayBase::CopyComTagsContainer& cctc = *recv_cctc[i];
                        for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
                                it != cctc.end(); ++it)
                        {
                            //if(f == local_index(mf,bxasc->r_con.rcv[i].nd)) //LocalIndex
                            if(mf.IndexArray()[f] == it->dstIndex)
                            {
                                rg->rMap[f]->r_con.nrcv++;
                            }
                        }
                    }
                    //rg->sMap[f]->r_con.sndLock = new omp_lock_t;
                    //rg->rMap[f]->r_con.rcvLock = new omp_lock_t;
                    //omp_init_lock(rg->sMap[f]->r_con.sndLock);
                    //omp_init_lock(rg->rMap[f]->r_con.rcvLock);
                    rg->sMap[f]->r_con.snd = new RemoteCommDescriptor[rg->sMap[f]->r_con.nsnd];
                    rg->rMap[f]->r_con.rcv = new RemoteCommDescriptor[rg->rMap[f]->r_con.nrcv];
                }
            }
        }
    } // omp parallel
    //std::cout<< "counting done " <<std::endl;
    //    !!touch data to bind pages to the NUMA node
//#pragma omp parallel shared(mf, numfabs, TheFB, recv_cctc, send_cctc)
    {
//        int tg = WorkerThread::perilla_wid();

        //      std::cout<< "Barr 4- "<< tid <<" "<< tg << " " << WorkerThread::isTeamMasterThread(tid) << std::endl;

        //      std::cout<< "Barr 5" <<std::endl;
        int fg, scnt, dcnt;

//#pragma omp for
        for(int f=0; f<numfabs; f++)
        {
            //fg = f % (omp_get_num_threads()/perilla::NUM_THREADS_PER_TEAM);
            //if((fg == tg) && ((tid%perilla::NUM_THREADS_PER_TEAM)==1))

            //if((fg == tg) && ((tid%perilla::NUM_THREADS_PER_TEAM)==0))
 //           if(WorkerThread::isMyRegion(tg,f) && WorkerThread::perilla_isMasterWorkerThread())
            {
                rg->lMap[f]->l_con.localBarrier = new Barrier(perilla::NUM_THREADS_PER_TEAM-1);
                // !create local communication meta data for sources and destinations
                scnt = -1;
                dcnt = -1;
                //for(int i=0; i<bxasc->l_con.ncpy; i++)
                for(int i=0; i<n_loc_mf; i++)
                {
                    const FabArrayBase::CopyComTag& tag = (*TheFB.m_LocTags)[i];
                    BL_ASSERT(mf.distributionMap[tag.dstIndex] == myProc);
                    BL_ASSERT(mf.distributionMap[tag.srcIndex] == myProc);

                    //if(f == local_index(mf,bxasc->l_con.cpy[i].ns)) //LocalIndex
                    if(mf.IndexArray()[f] == tag.srcIndex)
                    {
                        scnt++;
                        //omp_init_lock(&(rg->lMap[f]->l_con.scpy[scnt].ghostLock));
                        rg->lMap[f]->l_con.scpy[scnt].ns = mf.localindex(tag.srcIndex); //local_index(mf,bxasc->l_con.cpy[i].ns); //LocalIndex
                        rg->lMap[f]->l_con.scpy[scnt].nd = mf.localindex(tag.dstIndex); //local_index(mf,bxasc->l_con.cpy[i].nd); //LocalIndex
                        rg->lMap[f]->l_con.scpy[scnt].sbx = tag.sbox; //bxasc->l_con.cpy[i].sbx;
                        rg->lMap[f]->l_con.scpy[scnt].dbx = tag.dbox; //bxasc->l_con.cpy[i].dbx;                    
                        // !create queues for ghost cells
                        //call queue_init(mf%fbs(f)%l_con%scpy(scnt)%pQueue)
                        //call queue_init(mf%fbs(f)%l_con%scpy(scnt)%recycleQueue)

                        int psize = tag.sbox.numPts() * mf.nComp(); //---------------------------------------------------------------????????????????
                        /*
                           p => dataptr(mf%fbs(f), mf%fbs(f)%l_con%scpy(scnt)%sbx, 1, mf%nc)
                           s1= size(p,1)
                           s2= size(p,2)
                           s3= size(p,3)
                           s4= size(p,4)
                           s1*s2*s3*s4
                         */
                        for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
                        {
                            Package *tmpPkg = new Package(psize);
                            void* local_ptr= tmpPkg->databuf.local();//(static_cast<upcxx::global_ptr<void> >(tmpPkg->databuf)).local();
                            for(int j=0; j<psize; j++)
                                ((double*)local_ptr)[j]= 0;
                            rg->lMap[f]->l_con.scpy[scnt].pQueue.enqueue(tmpPkg);
                        }
                        for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
                            rg->lMap[f]->l_con.scpy[scnt].recycleQueue.enqueue(rg->lMap[f]->l_con.scpy[scnt].pQueue.dequeue());
                    }
                    //if(f == local_index(mf,bxasc->l_con.cpy[i].nd)) //LocalIndex
                    if(mf.IndexArray()[f] == tag.dstIndex)
                    {
                        dcnt++;
                        rg->lMap[f]->l_con.dcpy[dcnt].ns = mf.localindex(tag.srcIndex); //local_index(mf,bxasc->l_con.cpy[i].ns); //LocalIndex
                        rg->lMap[f]->l_con.dcpy[dcnt].nd = mf.localindex(tag.dstIndex); //local_index(mf,bxasc->l_con.cpy[i].nd); //LocalIndex
                        rg->lMap[f]->l_con.dcpy[dcnt].sbx = tag.sbox; //bxasc->l_con.cpy[i].sbx;
                        rg->lMap[f]->l_con.dcpy[dcnt].dbx = tag.dbox; //bxasc->l_con.cpy[i].dbx;                    
                        //call queue_init(mf%fbs(f)%l_con%dcpy(dcnt)%pQueue)
                        //call queue_init(mf%fbs(f)%l_con%dcpy(dcnt)%recycleQueue)
                        int psize = tag.dbox.numPts() * mf.nComp(); //---------------------------------------------------------------????????????????
                        /*
                           p => dataptr(mf%fbs(f), mf%fbs(f)%l_con%dcpy(dcnt)%dbx, 1, mf%nc)
                           s1= size(p,1)
                           s2= size(p,2)
                           s3= size(p,3)
                           s4= size(p,4)
                           s1*s2*s3*s4
                         */
                        for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
                        {
                            Package *tmpPkg = new Package(psize);
                            void* local_ptr= tmpPkg->databuf.local();//(static_cast<upcxx::global_ptr<void> >(tmpPkg->databuf)).local();
                            for(int j=0; j<psize; j++)
                                ((double*)local_ptr)[j] =0;
                            rg->lMap[f]->l_con.dcpy[dcnt].pQueue.enqueue(tmpPkg);
                        }
                        for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
                            rg->lMap[f]->l_con.dcpy[dcnt].recycleQueue.enqueue(rg->lMap[f]->l_con.dcpy[dcnt].pQueue.dequeue());
                    }
                } // for(i<n_loc_mf)
                //std::cout<< scnt << " " << dcnt << std::endl;
            }
        }// for(f<numfabs)

//#pragma omp barrier       

 //       if(WorkerThread::perilla_isMasterWorkerThread() && tg==0)
//#pragma omp for
            for(int f=0; f<numfabs; f++)
            {
                for(int i=0; i<rg->lMap[f]->l_con.nscpy; i++)
                    for(int j=0; j<rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.ndcpy; j++)
                        if(rg->lMap[f]->l_con.scpy[i].dbx == rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dcpy[j].dbx)
                            rg->lMap[f]->l_con.scpy[i].dPartner = j;

                for(int i=0; i<rg->lMap[f]->l_con.ndcpy; i++)
                    for(int j=0; j<rg->lMap[rg->lMap[f]->l_con.dcpy[i].ns]->l_con.nscpy; j++)
                        if(rg->lMap[f]->l_con.dcpy[i].dbx == rg->lMap[rg->lMap[f]->l_con.dcpy[i].ns]->l_con.scpy[j].dbx)
                            rg->lMap[f]->l_con.dcpy[i].sPartner = j;
            }
    }
    if(np == 1) return;

    //std::cout<< "local init done" <<std::endl;

//#pragma omp parallel shared(rg, mf, numfabs)
    {
        int tg = WorkerThread::perilla_wid();
        int fg, nsnd, nrcv;

        for(int f=0; f<numfabs; f++)
        {
            //fg = f % (omp_get_num_threads()/perilla::NUM_THREADS_PER_TEAM);
            //if((fg == tg) && ((tid%perilla::NUM_THREADS_PER_TEAM)==0))
 //           if(WorkerThread::isMyRegion(tg,f) && WorkerThread::perilla_isMasterWorkerThread())
            {
                //rg->lMap[f]->r_con.sndLock = new omp_lock_t;
                //rg->lMap[f]->r_con.rcvLock = new omp_lock_t;
                //omp_init_lock(rg->lMap[f]->r_con.sndLock);
                //omp_init_lock(rg->lMap[f]->r_con.rcvLock);
                //rg->lMap[f]->r_con.snd = new RemoteCommDescriptor[rg->lMap[f]->r_con.nsnd];
                //rg->lMap[f]->r_con.rcv = new RemoteCommDescriptor[rg->lMap[f]->r_con.nrcv];
                nrcv= -1;
                //for(int i=0; i<bxasc->r_con.nrcv; i++)
                for(int i=0; i<n_rcvs_mf; i++)
                {
                    const FabArrayBase::CopyComTagsContainer& cctc = *recv_cctc[i];
                    int pr = recv_pr[i];
                    for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
                            it != cctc.end(); ++it)
                    {
                        //if(f == local_index(mf,bxasc->r_con.rcv[i].nd)) //LocalIndex
                        if(mf.IndexArray()[f] == it->dstIndex)
                        {
                            nrcv++;
                            rg->lMap[f]->r_con.rcv[nrcv].ns = it->srcIndex; //bxasc->r_con.rcv[i].ns;
                            //rg->lMap[f]->r_con.rcv[nrcv].lnd = ; //local_index(mf,bxasc->r_con.rcv[i].nd); // not used anywhere so deferred ---------????????
                            //rg->lMap[f]->r_con.rcv[nrcv].lns = -1; //undefined
                            rg->lMap[f]->r_con.rcv[nrcv].nd = it->dstIndex; //bxasc->r_con.rcv[i].nd;
                            rg->lMap[f]->r_con.rcv[nrcv].lnd = mf.localindex(it->dstIndex);
                            rg->lMap[f]->r_con.rcv[nrcv].lns = mf.localindex(it->srcIndex);
                            rg->lMap[f]->r_con.rcv[nrcv].sbx = it->sbox; //bxasc->r_con.rcv[i].sbx;
                            rg->lMap[f]->r_con.rcv[nrcv].dbx = it->dbox; //bxasc->r_con.rcv[i].dbx;
                            rg->lMap[f]->r_con.rcv[nrcv].pr = pr; //bxasc->r_con.rcv[i].pr;
                            rg->lMap[f]->r_con.rcv[nrcv].cnt = 0;
                            //!create queues for ghost cells
                            //call queue_init(mf%fbs(f)%r_con%rcv(nrcv)%pQueue)
                            //call queue_init(mf%fbs(f)%r_con%rcv(nrcv)%recycleQueue)
                            int psize = it->sbox.numPts() * mf.nComp(); //-----------------------------------------------------------????????????????
                            /*
                               p => dataptr(mf%fbs(f), mf%fbs(f)%r_con%rcv(nrcv)%dbx, 1, mf%nc)
                               s1= size(p,1)
                               s2= size(p,2)
                               s3= size(p,3)
                               s4= size(p,4)
                               s1*s2*s3*s4
                             */
                            rg->lMap[f]->r_con.rcv[nrcv].sz = psize;
                            for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
                            {
                                Package *tmpPkg = new Package(psize);
                                void* local_ptr= tmpPkg->databuf.local();//(static_cast<upcxx::global_ptr<void> >(tmpPkg->databuf)).local();
                                for(int j=0; j<psize; j++)
                                    ((double*)local_ptr)[j]=0;
                                rg->lMap[f]->r_con.rcv[nrcv].pQueue.enqueue(tmpPkg);
                            }
                            for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
                                rg->lMap[f]->r_con.rcv[nrcv].recycleQueue.enqueue(rg->lMap[f]->r_con.rcv[nrcv].pQueue.dequeue());
                        }
                    }
                } // for(i<n_rcvs_mf)
                nsnd = -1;
                //for(int i=0; i<bxasc->r_con.nsnd; i++)
                for(int i=0; i<n_snds_mf; i++)
                {
                    const FabArrayBase::CopyComTagsContainer& cctc = *send_cctc[i];
                    int pr = send_pr[i];
                    for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
                            it != cctc.end(); ++it)
                    {
                        //if(f == local_index(mf,bxasc->r_con.snd[i].ns)) //LocalIndex
                        if(mf.IndexArray()[f] == it->srcIndex )
                        {
                            nsnd++;
                            rg->lMap[f]->r_con.snd[nsnd].ns = it->srcIndex; //bxasc->r_con.snd[i].ns;
                            rg->lMap[f]->r_con.snd[nsnd].nd = it->dstIndex; //bxasc->r_con.snd[i].nd;
                            //rg->lMap[f]->r_con.snd[nsnd].lns = ; //local_index(mf,bxasc->r_con.snd[i].ns); //not used anywhere so deferred ------?????????
                            //rg->lMap[f]->r_con.snd[nsnd].lnd = -1; //undefined
                            rg->lMap[f]->r_con.snd[nsnd].lns = mf.localindex(it->srcIndex);
                            rg->lMap[f]->r_con.snd[nsnd].lnd = mf.localindex(it->dstIndex);
                            rg->lMap[f]->r_con.snd[nsnd].sbx = it->sbox; //bxasc->r_con.snd[i].sbx;
                            rg->lMap[f]->r_con.snd[nsnd].dbx = it->dbox; //bxasc->r_con.snd[i].dbx;
                            rg->lMap[f]->r_con.snd[nsnd].pr = pr; //bxasc->r_con.snd[i].pr;
                            rg->lMap[f]->r_con.snd[nsnd].cnt = 0;
                            //!create queues for ghost cells            
                            //call queue_init(mf%fbs(f)%r_con%snd(nsnd)%pQueue)
                            //call queue_init(mf%fbs(f)%r_con%snd(nsnd)%recycleQueue)
                            int psize = it->sbox.numPts() * mf.nComp(); //---------------------------------------------------------------????????????????
                            /*
                               p => dataptr(mf%fbs(f), mf%fbs(f)%r_con%snd(nsnd)%sbx, 1, mf%nc)
                               s1= size(p,1)
                               s2= size(p,2)
                               s3= size(p,3)
                               s4= size(p,4)
                               s1*s2*s3*s4
                             */
                            rg->lMap[f]->r_con.snd[nsnd].sz = psize;
                            for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
                            {
                                Package *tmpPkg = new Package(psize);
                                void* local_ptr= tmpPkg->databuf.local(); //(static_cast<upcxx::global_ptr<void> >(tmpPkg->databuf)).local();                  
                                for(int j=0; j<psize; j++)
                                    ((double*)local_ptr)[j]=0;
                                rg->lMap[f]->r_con.snd[nsnd].pQueue.enqueue(tmpPkg);
                            }
                            for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
                                rg->lMap[f]->r_con.snd[nsnd].recycleQueue.enqueue(rg->lMap[f]->r_con.snd[nsnd].pQueue.dequeue());

                            //std::cout<< "RQ f "<< f << " i "<< nsnd <<std::endl;
                        }
                    }
                } // for(i<n_snds_mf)
                //std::cout<< "tid "<< tid << " f "<< f << " nfabs "<< numfabs <<std::endl;
            }// if(fg==tg...)
            //#pragma omp barrier
        }//for(f<numfabs)

        //std::cout<< "Barr 1 tid " << tid <<std::endl;

//#pragma omp barrier      //----------------------------------- Barrier ------------------------------------------      

        //if(tid == 0)
 //       if(WorkerThread::perilla_isMasterWorkerThread() && tg==0)
        {
            for(int f=0; f<numfabs; f++)
            {
                for(int i=0; i<rg->lMap[f]->r_con.nsnd; i++)
                {
                    rg->sMap[f]->r_con.snd[i].ns = rg->lMap[f]->r_con.snd[i].ns;
                    rg->sMap[f]->r_con.snd[i].nd = rg->lMap[f]->r_con.snd[i].nd;
                    rg->sMap[f]->r_con.snd[i].lns = rg->lMap[f]->r_con.snd[i].lns;
                    rg->sMap[f]->r_con.snd[i].lnd = rg->lMap[f]->r_con.snd[i].lnd;
                    rg->sMap[f]->r_con.snd[i].r_gid = rg->graphID-1;
                    rg->sMap[f]->r_con.snd[i].r_grids = rg->numFabs;
                    rg->sMap[f]->r_con.snd[i].sbx = rg->lMap[f]->r_con.snd[i].sbx;
                    rg->sMap[f]->r_con.snd[i].dbx = rg->lMap[f]->r_con.snd[i].dbx;
                    rg->sMap[f]->r_con.snd[i].pr = rg->lMap[f]->r_con.snd[i].pr;
                    rg->sMap[f]->r_con.snd[i].sz = rg->lMap[f]->r_con.snd[i].sz;
                    rg->sMap[f]->r_con.snd[i].cnt = 0;
                    rg->lMap[f]->r_con.snd[i].cnt = 0;

                    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
                    {
                        Package *tmpPkg = new Package(rg->lMap[f]->r_con.snd[i].sz);
                        void* local_ptr= tmpPkg->databuf.local(); //(static_cast<upcxx::global_ptr<void> >(tmpPkg->databuf)).local();
                        for(int j=0; j<rg->lMap[f]->r_con.snd[i].sz; j++)
                            ((double*)local_ptr)[j]=0;
                        rg->sMap[f]->r_con.snd[i].pQueue.enqueue(tmpPkg);
                    }
                    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
                        rg->sMap[f]->r_con.snd[i].recycleQueue.enqueue(rg->sMap[f]->r_con.snd[i].pQueue.dequeue());
                }
                for(int i=0; i<rg->lMap[f]->r_con.nrcv; i++)
                {
                    rg->rMap[f]->r_con.rcv[i].ns = rg->lMap[f]->r_con.rcv[i].ns;
                    rg->rMap[f]->r_con.rcv[i].nd = rg->lMap[f]->r_con.rcv[i].nd;
                    rg->rMap[f]->r_con.rcv[i].lns = rg->lMap[f]->r_con.rcv[i].lns;
                    rg->rMap[f]->r_con.rcv[i].lnd = rg->lMap[f]->r_con.rcv[i].lnd;
                    rg->rMap[f]->r_con.rcv[i].r_gid = rg->graphID-1;
                    rg->rMap[f]->r_con.rcv[i].r_grids = rg->numFabs;
                    rg->rMap[f]->r_con.rcv[i].sbx = rg->lMap[f]->r_con.rcv[i].sbx;
                    rg->rMap[f]->r_con.rcv[i].dbx = rg->lMap[f]->r_con.rcv[i].dbx;
                    rg->rMap[f]->r_con.rcv[i].pr = rg->lMap[f]->r_con.rcv[i].pr;
                    rg->rMap[f]->r_con.rcv[i].sz = rg->lMap[f]->r_con.rcv[i].sz;
                    rg->rMap[f]->r_con.rcv[i].cnt = 0;
                    rg->lMap[f]->r_con.rcv[i].cnt = 0;

                    if(Perilla::genTags)
                    {
                        try
                        {
                            int rcv_pr = rg->rMap[f]->r_con.rcv[i].pr;
                            int dstIndex = rg->rMap[f]->r_con.rcv[i].nd;
                            int srcIndex = rg->rMap[f]->r_con.rcv[i].ns;
                            int psize = rg->rMap[f]->r_con.rcv[i].sz;
                            std::map<int,int>::iterator itr = tagMap[rcv_pr][rg->graphID-1][dstIndex][srcIndex].find(psize);
                            if( itr != tagMap[rcv_pr][rg->graphID-1][dstIndex][srcIndex].end())
                            {
                                //rg->rCopyMapHead->map[f]->r_con.rcv[dcnt].lnd = itr->second;
                            }
                            else
                            {
                                tagMap[rcv_pr][rg->graphID-1][dstIndex][srcIndex][psize] = Perilla::uTags++;
                                //rg->rCopyMapHead->map[f]->r_con.rcv[dcnt].lnd = Perilla::uTags++;
                                std::map<int,int>::iterator itr2 = pTagCnt[rcv_pr].find(rg->graphID-1);
                                if(itr2 != pTagCnt[rcv_pr].end())
                                    pTagCnt[rcv_pr][rg->graphID-1] = pTagCnt[rcv_pr][rg->graphID-1] + 1;
                                else
                                    pTagCnt[rcv_pr][rg->graphID-1] = 1;
                            }
                        }
                        catch(std::exception& e)
                        {
                            std::cout <<"Inside tagGeneration gID "<< rg->graphID <<" "<< e.what() << '\n';
                        }
                    }
                    //tagMap[rcv_pr][rg->graphID][it->dstIndex][it->srcIndex] = pTagCnt[rcv_pr];                                  

                    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
                    {
                        Package *tmpPkg = new Package(rg->lMap[f]->r_con.rcv[i].sz);
                        void* local_ptr= tmpPkg->databuf.local(); //(static_cast<upcxx::global_ptr<void> >(tmpPkg->databuf)).local();
                        for(int j=0; j<rg->lMap[f]->r_con.rcv[i].sz; j++)
                            ((double*)local_ptr)[j]=0;
                        rg->rMap[f]->r_con.rcv[i].pQueue.enqueue(tmpPkg);
                    }
                    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
                        rg->rMap[f]->r_con.rcv[i].recycleQueue.enqueue(rg->rMap[f]->r_con.rcv[i].pQueue.dequeue());
                }
            }
        }// if(tid==0)

    }// omp parallel    
}// multifabBuildFabCon


#endif

#if 0
#ifdef USE_PERILLA_PTHREADS
Array<const FabArrayBase::CopyComTagsContainer*> send_cctc1;
Array<int> send_pr1;
Array<const FabArrayBase::CopyComTagsContainer*> recv_cctc1;
Array<int> recv_pr1;
#endif

void Perilla::multifabBuildFabCon(RegionGraph* rg, const MultiFab& mf, const Periodicity& period)
{
    int np = ParallelDescriptor::NProcs();
    int myProc = ParallelDescriptor::MyProc();
    int numfabs = mf.IndexArray().size();
    bool cross = false;
    const FabArrayBase::FB& TheFB = mf.getFB(period, false, false);
    const int n_loc_mf = TheFB.m_LocTags->size();
    const int n_snds_mf = TheFB.m_SndTags->size();
    const int n_rcvs_mf = TheFB.m_RcvTags->size();


#ifndef USE_PERILLA_PTHREADS
    Array<const FabArrayBase::CopyComTagsContainer*> send_cctc1;
    Array<int> send_pr1;
    Array<const FabArrayBase::CopyComTagsContainer*> recv_cctc1;
    Array<int> recv_pr1;
#endif


#ifdef USE_PERILLA_PTHREADS
    if(perilla::isMasterThread){
#endif
	send_cctc1.reserve(n_snds_mf);

	for (FabArrayBase::MapOfCopyComTagContainers::const_iterator m_it = TheFB.m_SndTags->begin(),
		m_End = TheFB.m_SndTags->end();
		m_it != m_End;
		++m_it)
	{
	    if(m_it->first != myProc)      // Not destined to me.
	    {
		send_pr1.push_back(m_it->first);
		send_cctc1.push_back(&(m_it->second));
	    }
	}

	recv_cctc1.reserve(n_rcvs_mf);

	for (FabArrayBase::MapOfCopyComTagContainers::const_iterator m_it = TheFB.m_RcvTags->begin(),
		m_End = TheFB.m_RcvTags->end();
		m_it != m_End;
		++m_it)
	{
	    if(m_it->first != myProc)      // I am not the source for this receipt
	    {
		recv_pr1.push_back(m_it->first);
		recv_cctc1.push_back(&(m_it->second));
	    }
	}
#ifdef USE_PERILLA_PTHREADS
    }
#endif

#ifndef USE_PERILLA_PTHREADS
#pragma omp parallel shared(rg, mf, numfabs, np, TheFB, recv_cctc1, send_cctc1)
#endif
    {
	int tg = WorkerThread::perilla_wid();
	int fg;
	if(WorkerThread::perilla_isCommunicationThread())	
	{	  
	    //bool cc = !mf->is_nodal(); //  cc = multifab_cell_centered_q(mf)
	    //mf->sMap.reserve(numfabs);
	    //mf->rMap.reserve(numfabs);
	    //std::cout<< "Allocating sMap and rMap" <<std::endl;
	    rg->alloc_lMap(mf);	  
	    rg->alloc_sMap(mf);
	    rg->alloc_rMap(mf);
	}
	perilla::syncAllThreads();
	//#pragma omp barrier      
	//if(tid==0)            	          
	{	  
	    //bool cc = !mf->is_nodal(); //  cc = multifab_cell_centered_q(mf)
	    //mf->sMap.reserve(numfabs);
	    //mf->rMap.reserve(numfabs);
	    for(int f=0; f<numfabs; f++) //	   !create local communication metadata for each fab
	    {
		if(WorkerThread::isMyRegion(tg,f) && WorkerThread::perilla_isMasterWorkerThread())
		{	      
		    rg->lMap[f]->l_con.nscpy = 0;

		    //for(int i=0; i<bxasc->l_con.ncpy; i++)
		    for(int i=0; i<n_loc_mf; i++)
		    {
			const FabArrayBase::CopyComTag& tag = (*TheFB.m_LocTags)[i];

			//std::cout << tag.srcIndex << " " << tag.dstIndex << " " <<tag.sbox.smallEnd() <<" "<< tag.sbox.bigEnd() << std::endl;

			BL_ASSERT(mf.distributionMap[tag.dstIndex] == myProc);
			BL_ASSERT(mf.distributionMap[tag.srcIndex] == myProc);
			//get(tag.dstIndex).copy(get(tag.srcIndex),tag.box,scomp,tag.box,scomp,ncomp);
			//if(f == local_index(mf,bxasc->l_con.cpy[i].ns)) //LocalIndex
			if(mf.IndexArray()[f] == tag.srcIndex)
			    rg->lMap[f]->l_con.nscpy++;		  
			//if(f == local_index(mf,bxasc->l_con.cpy[i].nd)) //LocalIndex
			if(mf.IndexArray()[f] == tag.dstIndex)
			    rg->lMap[f]->l_con.ndcpy++;
		    }
		    /*
		       if(rg->lMap[f]->l_con.nscpy+rg->lMap[f]->l_con.ndcpy != n_loc_mf)
		       std::cout<< "Diff in Sum " << rg->lMap[f]->l_con.nscpy << " " <<rg->lMap[f]->l_con.ndcpy << " " << n_loc_mf <<std::endl;
		       BL_ASSERT(rg->lMap[f]->l_con.nscpy+rg->lMap[f]->l_con.ndcpy == n_loc_mf);
		     */
		}
	    }
	}
	perilla::syncAllThreads();
	//now we know how many copying segments each fab owns as source and destination allocate memory for metadata   
	for(int f=0; f<numfabs; f++)
	{
	    //fg = f % (omp_get_num_threads()/perilla::NUM_THREADS_PER_TEAM);   /// need to check if computing correct ???????
	    //if((fg == tg) && ((tid%perilla::NUM_THREADS_PER_TEAM)==1))
	    if(WorkerThread::isMyRegion(tg,f) && WorkerThread::perilla_isMasterWorkerThread())	    
	    {
		//omp_init_lock(&(rg->lMap[f]->l_con.sLock));
		//omp_init_lock(&(rg->lMap[f]->l_con.dLock));
		//omp_init_lock(&(rg->lMap[f]->l_con.ghostLock));

		//std::cout<< "MF l_con nscpy " <<rg->lMap[f]->l_con.nscpy << " ndcpy " << rg->lMap[f]->l_con.ndcpy <<std::endl;

		rg->lMap[f]->l_con.scpy = new LocalCopyDescriptor[rg->lMap[f]->l_con.nscpy];
		rg->lMap[f]->l_con.dcpy = new LocalCopyDescriptor[rg->lMap[f]->l_con.ndcpy];
		rg->lMap[f]->l_con.scpyCnt = 0;
		rg->lMap[f]->l_con.dcpyCnt = 0;
	    }
	}
	perilla::syncAllThreads();
	if(np > 1)
	{
	    for(int f=0; f<numfabs; f++)
	    {
		if(WorkerThread::perilla_isMasterWorkerThread() && WorkerThread::isMyRegion(tg,f))      
		{
		    rg->lMap[f]->r_con.nrcv = 0;
		    rg->lMap[f]->r_con.nsnd = 0;
		    rg->lMap[f]->r_con.firingRuleCnt = 0;

		    //for(int i=0; i<bxasc->r_con.nsnd; i++)
		    for(int i=0; i<n_snds_mf; i++)
		    {
			const FabArrayBase::CopyComTagsContainer& cctc = *send_cctc1[i];
			for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
				it != cctc.end(); ++it)
			{
			    //if(f == local_index(mf,bxasc->r_con.snd[i].ns)) //LocalIndex
			    if(mf.IndexArray()[f] == it->srcIndex)
			    {
				rg->lMap[f]->r_con.nsnd++;
			    }
			}
		    }
		    //for(int i=0; i<bxasc->r_con.nrcv; i++)
		    for(int i=0; i<n_rcvs_mf; i++)
		    {
			const FabArrayBase::CopyComTagsContainer& cctc = *recv_cctc1[i];
			for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
				it != cctc.end(); ++it)
			{
			    //if(f == local_index(mf,bxasc->r_con.rcv[i].nd)) //LocalIndex
			    if(mf.IndexArray()[f] == it->dstIndex)
			    {
				rg->lMap[f]->r_con.nrcv++;
			    }
			}
		    }
		    //rg->sMap[f]->r_con.sndLock = new omp_lock_t;
		    //rg->rMap[f]->r_con.rcvLock = new omp_lock_t;
		    //omp_init_lock(rg->sMap[f]->r_con.sndLock);
		    //omp_init_lock(rg->rMap[f]->r_con.rcvLock);
		    rg->lMap[f]->r_con.snd = new RemoteCommDescriptor[rg->lMap[f]->r_con.nsnd];
		    rg->lMap[f]->r_con.rcv = new RemoteCommDescriptor[rg->lMap[f]->r_con.nrcv];
		}		
	    }	
	    if(WorkerThread::perilla_isMasterWorkerThread() && tg==0)      
	    {
		for(int f=0; f<numfabs; f++)
		{
		    rg->rMap[f]->r_con.nrcv = 0;
		    rg->sMap[f]->r_con.nsnd = 0;

		    //for(int i=0; i<bxasc->r_con.nsnd; i++)
		    for(int i=0; i<n_snds_mf; i++)
		    {
			const FabArrayBase::CopyComTagsContainer& cctc = *send_cctc1[i];
			for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
				it != cctc.end(); ++it)
			{
			    //if(f == local_index(mf,bxasc->r_con.snd[i].ns)) //LocalIndex
			    if(mf.IndexArray()[f] == it->srcIndex)
			    {
				rg->sMap[f]->r_con.nsnd++;
			    }
			}
		    }
		    //for(int i=0; i<bxasc->r_con.nrcv; i++)
		    for(int i=0; i<n_rcvs_mf; i++)
		    {
			const FabArrayBase::CopyComTagsContainer& cctc = *recv_cctc1[i];
			for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
				it != cctc.end(); ++it)
			{
			    //if(f == local_index(mf,bxasc->r_con.rcv[i].nd)) //LocalIndex
			    if(mf.IndexArray()[f] == it->dstIndex)
			    {
				rg->rMap[f]->r_con.nrcv++;
			    }
			}
		    }
		    //rg->sMap[f]->r_con.sndLock = new omp_lock_t;
		    //rg->rMap[f]->r_con.rcvLock = new omp_lock_t;
		    //omp_init_lock(rg->sMap[f]->r_con.sndLock);
		    //omp_init_lock(rg->rMap[f]->r_con.rcvLock);
		    rg->sMap[f]->r_con.snd = new RemoteCommDescriptor[rg->sMap[f]->r_con.nsnd];
		    rg->rMap[f]->r_con.rcv = new RemoteCommDescriptor[rg->rMap[f]->r_con.nrcv];
		}
	    }
	}
    } // omp parallel
    //std::cout<< "counting done " <<std::endl;
    //    !!touch data to bind pages to the NUMA node

#ifndef USE_PERILLA_PTHREADS
#pragma omp parallel shared(mf, numfabs, TheFB, recv_cctc1, send_cctc1)
#endif
    {
	int tg = WorkerThread::perilla_wid();

	//      std::cout<< "Barr 4- "<< tid <<" "<< tg << " " << WorkerThread::isTeamMasterThread(tid) << std::endl;

	//      std::cout<< "Barr 5" <<std::endl;
	int fg, scnt, dcnt;

	for(int f=0; f<numfabs; f++)
	{
	    //fg = f % (omp_get_num_threads()/perilla::NUM_THREADS_PER_TEAM);
	    //if((fg == tg) && ((tid%perilla::NUM_THREADS_PER_TEAM)==1))

	    //if((fg == tg) && ((tid%perilla::NUM_THREADS_PER_TEAM)==0))
	    if(WorkerThread::isMyRegion(tg,f) && WorkerThread::perilla_isMasterWorkerThread())
	    {
		rg->lMap[f]->l_con.localBarrier = new Barrier(perilla::NUM_THREADS_PER_TEAM-1);
		// !create local communication meta data for sources and destinations
		scnt = -1;
		dcnt = -1;
		//for(int i=0; i<bxasc->l_con.ncpy; i++)
		for(int i=0; i<n_loc_mf; i++)
		{
		    const FabArrayBase::CopyComTag& tag = (*TheFB.m_LocTags)[i];
		    BL_ASSERT(mf.distributionMap[tag.dstIndex] == myProc);
		    BL_ASSERT(mf.distributionMap[tag.srcIndex] == myProc);

		    //if(f == local_index(mf,bxasc->l_con.cpy[i].ns)) //LocalIndex
		    if(mf.IndexArray()[f] == tag.srcIndex)
		    {
			scnt++;
			//omp_init_lock(&(rg->lMap[f]->l_con.scpy[scnt].ghostLock));
			rg->lMap[f]->l_con.scpy[scnt].ns = mf.localindex(tag.srcIndex); //local_index(mf,bxasc->l_con.cpy[i].ns); //LocalIndex
			rg->lMap[f]->l_con.scpy[scnt].nd = mf.localindex(tag.dstIndex); //local_index(mf,bxasc->l_con.cpy[i].nd); //LocalIndex
			rg->lMap[f]->l_con.scpy[scnt].sbx = tag.sbox; //bxasc->l_con.cpy[i].sbx;
			rg->lMap[f]->l_con.scpy[scnt].dbx = tag.dbox; //bxasc->l_con.cpy[i].dbx;		    
			// !create queues for ghost cells
			//call queue_init(mf%fbs(f)%l_con%scpy(scnt)%pQueue)
			//call queue_init(mf%fbs(f)%l_con%scpy(scnt)%recycleQueue)
			int psize = tag.sbox.numPts() * mf.nComp(); //---------------------------------------------------------------????????????????
			/*
			   p => dataptr(mf%fbs(f), mf%fbs(f)%l_con%scpy(scnt)%sbx, 1, mf%nc)
			   s1= size(p,1)
			   s2= size(p,2)
			   s3= size(p,3)
			   s4= size(p,4)
			   s1*s2*s3*s4
			 */
			for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
			{
			    Package *tmpPkg = new Package(psize);
			    for(int j=0; j<psize; j++)
				tmpPkg->databuf[j] = 0;				  
			    rg->lMap[f]->l_con.scpy[scnt].pQueue.enqueue(tmpPkg);			
			}
			for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
			    rg->lMap[f]->l_con.scpy[scnt].recycleQueue.enqueue(rg->lMap[f]->l_con.scpy[scnt].pQueue.dequeue());
		    }	      
		    //if(f == local_index(mf,bxasc->l_con.cpy[i].nd)) //LocalIndex
		    if(mf.IndexArray()[f] == tag.dstIndex)		    
		    {
			dcnt++;		      		      
			rg->lMap[f]->l_con.dcpy[dcnt].ns = mf.localindex(tag.srcIndex); //local_index(mf,bxasc->l_con.cpy[i].ns); //LocalIndex
			rg->lMap[f]->l_con.dcpy[dcnt].nd = mf.localindex(tag.dstIndex); //local_index(mf,bxasc->l_con.cpy[i].nd); //LocalIndex
			rg->lMap[f]->l_con.dcpy[dcnt].sbx = tag.sbox; //bxasc->l_con.cpy[i].sbx;
			rg->lMap[f]->l_con.dcpy[dcnt].dbx = tag.dbox; //bxasc->l_con.cpy[i].dbx;		    
			//call queue_init(mf%fbs(f)%l_con%dcpy(dcnt)%pQueue)
			//call queue_init(mf%fbs(f)%l_con%dcpy(dcnt)%recycleQueue)
			int psize = tag.dbox.numPts() * mf.nComp(); //---------------------------------------------------------------????????????????
			/*
			   p => dataptr(mf%fbs(f), mf%fbs(f)%l_con%dcpy(dcnt)%dbx, 1, mf%nc)
			   s1= size(p,1)
			   s2= size(p,2)
			   s3= size(p,3)
			   s4= size(p,4)
			   s1*s2*s3*s4
			 */

			for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
			{
			    Package *tmpPkg = new Package(psize);
			    for(int j=0; j<psize; j++)
				tmpPkg->databuf[j] = 0;				  
			    rg->lMap[f]->l_con.dcpy[dcnt].pQueue.enqueue(tmpPkg);			
			}
			for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
			    rg->lMap[f]->l_con.dcpy[dcnt].recycleQueue.enqueue(rg->lMap[f]->l_con.dcpy[dcnt].pQueue.dequeue());
		    }
		} // for(i<n_loc_mf)
		//std::cout<< scnt << " " << dcnt << std::endl;
	    }
	}// for(f<numfabs)

	//#pragma omp barrier	  
	perilla::syncAllThreads();

	if(WorkerThread::perilla_isMasterWorkerThread() && tg==0)    
	    for(int f=0; f<numfabs; f++)
	    {
		for(int i=0; i<rg->lMap[f]->l_con.nscpy; i++)
		    for(int j=0; j<rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.ndcpy; j++)
			if(rg->lMap[f]->l_con.scpy[i].dbx == rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dcpy[j].dbx)
			    rg->lMap[f]->l_con.scpy[i].dPartner = j;

		for(int i=0; i<rg->lMap[f]->l_con.ndcpy; i++)
		    for(int j=0; j<rg->lMap[rg->lMap[f]->l_con.dcpy[i].ns]->l_con.nscpy; j++)
			if(rg->lMap[f]->l_con.dcpy[i].dbx == rg->lMap[rg->lMap[f]->l_con.dcpy[i].ns]->l_con.scpy[j].dbx)
			    rg->lMap[f]->l_con.dcpy[i].sPartner = j;	
	    }
    }
#ifdef USE_PERILLA_PTHREADS
    perilla::syncAllThreads();
#endif

    if(np == 1) return;

    //std::cout<< "local init done" <<std::endl;

#ifndef USE_PERILLA_PTHREADS
#pragma omp parallel shared(rg, mf, numfabs)
#endif
    {
	int tg = WorkerThread::perilla_wid();
	int fg, nsnd, nrcv;

	for(int f=0; f<numfabs; f++)
	{
	    //fg = f % (omp_get_num_threads()/perilla::NUM_THREADS_PER_TEAM);
	    //if((fg == tg) && ((tid%perilla::NUM_THREADS_PER_TEAM)==0))
	    if(WorkerThread::isMyRegion(tg,f) && WorkerThread::perilla_isMasterWorkerThread())
	    {
		//rg->lMap[f]->r_con.sndLock = new omp_lock_t;
		//rg->lMap[f]->r_con.rcvLock = new omp_lock_t;
		//omp_init_lock(rg->lMap[f]->r_con.sndLock);
		//omp_init_lock(rg->lMap[f]->r_con.rcvLock);
		//rg->lMap[f]->r_con.snd = new RemoteCommDescriptor[rg->lMap[f]->r_con.nsnd];
		//rg->lMap[f]->r_con.rcv = new RemoteCommDescriptor[rg->lMap[f]->r_con.nrcv];
		nrcv= -1;
		//for(int i=0; i<bxasc->r_con.nrcv; i++)
		for(int i=0; i<n_rcvs_mf; i++)
		{
		    const FabArrayBase::CopyComTagsContainer& cctc = *recv_cctc1[i];
		    int pr = recv_pr1[i];
		    for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
			    it != cctc.end(); ++it)
		    {	      
			//if(f == local_index(mf,bxasc->r_con.rcv[i].nd)) //LocalIndex


			if(mf.IndexArray()[f] == it->dstIndex)
			{
			    nrcv++;
			    rg->lMap[f]->r_con.rcv[nrcv].ns = it->srcIndex; //bxasc->r_con.rcv[i].ns;
			    //rg->lMap[f]->r_con.rcv[nrcv].lnd = ; //local_index(mf,bxasc->r_con.rcv[i].nd); // not used anywhere so deferred ---------????????
			    //rg->lMap[f]->r_con.rcv[nrcv].lns = -1; //undefined
			    rg->lMap[f]->r_con.rcv[nrcv].nd = it->dstIndex; //bxasc->r_con.rcv[i].nd;
			    rg->lMap[f]->r_con.rcv[nrcv].lnd = mf.localindex(it->dstIndex);
			    rg->lMap[f]->r_con.rcv[nrcv].lns = mf.localindex(it->srcIndex);
			    rg->lMap[f]->r_con.rcv[nrcv].sbx = it->sbox; //bxasc->r_con.rcv[i].sbx;
			    rg->lMap[f]->r_con.rcv[nrcv].dbx = it->dbox; //bxasc->r_con.rcv[i].dbx;
			    rg->lMap[f]->r_con.rcv[nrcv].pr = pr; //bxasc->r_con.rcv[i].pr;
			    rg->lMap[f]->r_con.rcv[nrcv].cnt = 0;		    
			    //!create queues for ghost cells
			    //call queue_init(mf%fbs(f)%r_con%rcv(nrcv)%pQueue)
			    //call queue_init(mf%fbs(f)%r_con%rcv(nrcv)%recycleQueue)
			    int psize = it->sbox.numPts() * mf.nComp(); //---------------------------------------------------------------????????????????
			    /*
			       p => dataptr(mf%fbs(f), mf%fbs(f)%r_con%rcv(nrcv)%dbx, 1, mf%nc)
			       s1= size(p,1)
			       s2= size(p,2)
			       s3= size(p,3)
			       s4= size(p,4)
			       s1*s2*s3*s4
			     */
			    rg->lMap[f]->r_con.rcv[nrcv].sz = psize;
			    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
			    {
				Package *tmpPkg = new Package(psize);
				for(int j=0; j<psize; j++)
				    tmpPkg->databuf[j] = 0;				  
				rg->lMap[f]->r_con.rcv[nrcv].pQueue.enqueue(tmpPkg);			
			    }
			    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
				rg->lMap[f]->r_con.rcv[nrcv].recycleQueue.enqueue(rg->lMap[f]->r_con.rcv[nrcv].pQueue.dequeue());
			}
		    }
		} // for(i<n_rcvs_mf)

		nsnd = -1;
		//for(int i=0; i<bxasc->r_con.nsnd; i++)
		for(int i=0; i<n_snds_mf; i++)
		{
		    const FabArrayBase::CopyComTagsContainer& cctc = *send_cctc1[i];
		    int pr = send_pr1[i];
		    for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
			    it != cctc.end(); ++it)
		    { 
			//if(f == local_index(mf,bxasc->r_con.snd[i].ns)) //LocalIndex
			if(mf.IndexArray()[f] == it->srcIndex )
			{
			    nsnd++;
			    rg->lMap[f]->r_con.snd[nsnd].ns = it->srcIndex; //bxasc->r_con.snd[i].ns;
			    rg->lMap[f]->r_con.snd[nsnd].nd = it->dstIndex; //bxasc->r_con.snd[i].nd;
			    //rg->lMap[f]->r_con.snd[nsnd].lns = ; //local_index(mf,bxasc->r_con.snd[i].ns); //not used anywhere so deferred ------?????????
			    //rg->lMap[f]->r_con.snd[nsnd].lnd = -1; //undefined
			    rg->lMap[f]->r_con.snd[nsnd].lns = mf.localindex(it->srcIndex);
			    rg->lMap[f]->r_con.snd[nsnd].lnd = mf.localindex(it->dstIndex);
			    rg->lMap[f]->r_con.snd[nsnd].sbx = it->sbox; //bxasc->r_con.snd[i].sbx;
			    rg->lMap[f]->r_con.snd[nsnd].dbx = it->dbox; //bxasc->r_con.snd[i].dbx;
			    rg->lMap[f]->r_con.snd[nsnd].pr = pr; //bxasc->r_con.snd[i].pr;
			    rg->lMap[f]->r_con.snd[nsnd].cnt = 0;
			    //!create queues for ghost cells		
			    //call queue_init(mf%fbs(f)%r_con%snd(nsnd)%pQueue)
			    //call queue_init(mf%fbs(f)%r_con%snd(nsnd)%recycleQueue)
			    int psize = it->sbox.numPts() * mf.nComp(); //---------------------------------------------------------------????????????????
			    /*
			       p => dataptr(mf%fbs(f), mf%fbs(f)%r_con%snd(nsnd)%sbx, 1, mf%nc)
			       s1= size(p,1)
			       s2= size(p,2)
			       s3= size(p,3)
			       s4= size(p,4)
			       s1*s2*s3*s4
			     */		    
			    rg->lMap[f]->r_con.snd[nsnd].sz = psize;
			    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
			    {
				Package *tmpPkg = new Package(psize);
				for(int j=0; j<psize; j++)
				    tmpPkg->databuf[j] = 0;				  
				rg->lMap[f]->r_con.snd[nsnd].pQueue.enqueue(tmpPkg);			
			    }
			    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
				rg->lMap[f]->r_con.snd[nsnd].recycleQueue.enqueue(rg->lMap[f]->r_con.snd[nsnd].pQueue.dequeue());

			    //std::cout<< "RQ f "<< f << " i "<< nsnd <<std::endl;
			}
		    }
		} // for(i<n_snds_mf)
		//std::cout<< "tid "<< tid << " f "<< f << " nfabs "<< numfabs <<std::endl;
	    }// if(fg==tg...)
	    //#pragma omp barrier
	}//for(f<numfabs)

	//std::cout<< "Barr 1 tid " << tid <<std::endl;

	//#pragma omp barrier      //----------------------------------- Barrier ------------------------------------------      
	perilla::syncAllThreads();

	//if(tid == 0)
	if(WorkerThread::perilla_isMasterWorkerThread() && tg==0)
	{
	    for(int f=0; f<numfabs; f++)	  
	    {
		for(int i=0; i<rg->lMap[f]->r_con.nsnd; i++)
		{
		    rg->sMap[f]->r_con.snd[i].ns = rg->lMap[f]->r_con.snd[i].ns;
		    rg->sMap[f]->r_con.snd[i].nd = rg->lMap[f]->r_con.snd[i].nd;
		    rg->sMap[f]->r_con.snd[i].lns = rg->lMap[f]->r_con.snd[i].lns;
		    rg->sMap[f]->r_con.snd[i].lnd = rg->lMap[f]->r_con.snd[i].lnd;
		    rg->sMap[f]->r_con.snd[i].r_gid = rg->graphID-1;
		    rg->sMap[f]->r_con.snd[i].r_grids = rg->numFabs;
		    rg->sMap[f]->r_con.snd[i].sbx = rg->lMap[f]->r_con.snd[i].sbx;
		    rg->sMap[f]->r_con.snd[i].dbx = rg->lMap[f]->r_con.snd[i].dbx;
		    rg->sMap[f]->r_con.snd[i].pr = rg->lMap[f]->r_con.snd[i].pr;
		    rg->sMap[f]->r_con.snd[i].sz = rg->lMap[f]->r_con.snd[i].sz;
		    rg->sMap[f]->r_con.snd[i].cnt = 0;
		    rg->lMap[f]->r_con.snd[i].cnt = 0;

		    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
		    {
			Package *tmpPkg = new Package(rg->lMap[f]->r_con.snd[i].sz);
			for(int j=0; j<rg->lMap[f]->r_con.snd[i].sz; j++)
			    tmpPkg->databuf[j] = 0;				  
			rg->sMap[f]->r_con.snd[i].pQueue.enqueue(tmpPkg);			
		    }
		    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
			rg->sMap[f]->r_con.snd[i].recycleQueue.enqueue(rg->sMap[f]->r_con.snd[i].pQueue.dequeue());
		}
		for(int i=0; i<rg->lMap[f]->r_con.nrcv; i++)
		{
		    rg->rMap[f]->r_con.rcv[i].ns = rg->lMap[f]->r_con.rcv[i].ns;
		    rg->rMap[f]->r_con.rcv[i].nd = rg->lMap[f]->r_con.rcv[i].nd;
		    rg->rMap[f]->r_con.rcv[i].lns = rg->lMap[f]->r_con.rcv[i].lns;
		    rg->rMap[f]->r_con.rcv[i].lnd = rg->lMap[f]->r_con.rcv[i].lnd;
		    rg->rMap[f]->r_con.rcv[i].r_gid = rg->graphID-1;
		    rg->rMap[f]->r_con.rcv[i].r_grids = rg->numFabs;
		    rg->rMap[f]->r_con.rcv[i].sbx = rg->lMap[f]->r_con.rcv[i].sbx;
		    rg->rMap[f]->r_con.rcv[i].dbx = rg->lMap[f]->r_con.rcv[i].dbx;
		    rg->rMap[f]->r_con.rcv[i].pr = rg->lMap[f]->r_con.rcv[i].pr;
		    rg->rMap[f]->r_con.rcv[i].sz = rg->lMap[f]->r_con.rcv[i].sz;
		    rg->rMap[f]->r_con.rcv[i].cnt = 0;
		    rg->lMap[f]->r_con.rcv[i].cnt = 0;

		    if(Perilla::genTags)
		    {
			try
			{
			    int rcv_pr = rg->rMap[f]->r_con.rcv[i].pr;
			    int dstIndex = rg->rMap[f]->r_con.rcv[i].nd;
			    int srcIndex = rg->rMap[f]->r_con.rcv[i].ns;
			    int psize = rg->rMap[f]->r_con.rcv[i].sz;
			    std::map<int,int>::iterator itr = tagMap[rcv_pr][rg->graphID-1][dstIndex][srcIndex].find(psize);
			    if( itr != tagMap[rcv_pr][rg->graphID-1][dstIndex][srcIndex].end())
			    {
				//rg->rCopyMapHead->map[f]->r_con.rcv[dcnt].lnd = itr->second;
			    }
			    else
			    {
				tagMap[rcv_pr][rg->graphID-1][dstIndex][srcIndex][psize] = Perilla::uTags++;
				//rg->rCopyMapHead->map[f]->r_con.rcv[dcnt].lnd = Perilla::uTags++;
				std::map<int,int>::iterator itr2 = pTagCnt[rcv_pr].find(rg->graphID-1);
				if(itr2 != pTagCnt[rcv_pr].end())
				    pTagCnt[rcv_pr][rg->graphID-1] = pTagCnt[rcv_pr][rg->graphID-1] + 1;
				else
				    pTagCnt[rcv_pr][rg->graphID-1] = 1;
			    }
			}
			catch(std::exception& e)
			{
			    std::cout <<"Inside tagGeneration gID "<< rg->graphID <<" "<< e.what() << '\n';
			}
		    }
		    //tagMap[rcv_pr][rg->graphID][it->dstIndex][it->srcIndex] = pTagCnt[rcv_pr];				  

		    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
		    {
			Package *tmpPkg = new Package(rg->lMap[f]->r_con.rcv[i].sz);
			for(int j=0; j<rg->lMap[f]->r_con.rcv[i].sz; j++)
			    tmpPkg->databuf[j] = 0;				  
			rg->rMap[f]->r_con.rcv[i].pQueue.enqueue(tmpPkg);			
		    }
		    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
			rg->rMap[f]->r_con.rcv[i].recycleQueue.enqueue(rg->rMap[f]->r_con.rcv[i].pQueue.dequeue());
		}
	    }
	}// if(tid==0)

    }// omp parallel    
}// multifabBuildFabCon
#endif

void Perilla::serviceLocalRequests(RegionGraph* rg, int tg)    
{
    int numfabs = rg->lMap.size();

    for(int f=0; f<numfabs; f++)
    {
	if(WorkerThread::isMyRegion(tg,f))
	{
	    //int lockSucceeded = pthread_mutex_trylock(&(rg->lMap[f]->l_con.sLock));
	    //if(lockSucceeded != 0) // 0-Fail, otherwise-Succeed
	    {		
		for(int i=0; i<rg->lMap[f]->l_con.nscpy; i++){
		    if(rg->lMap[f]->l_con.scpy[i].pQueue.queueSize()>0)
		    {
	    		pthread_mutex_lock(&(rg->lMap[f]->l_con.sLock));
			Package *sPackage = rg->lMap[f]->l_con.scpy[i].pQueue.dequeue();
			if(perilla::LAZY_PUSH)
			{
			    //  Implemetation deffered. Currently not required
			}
			pthread_mutex_lock(&(rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dLock));
			int dPartner = rg->lMap[f]->l_con.scpy[i].dPartner;
			if(dPartner == -1)
			    std::cout<< " Caution rQ size dPrtn "<< rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.ndcpy << " " << dPartner <<" graph ID " <<rg->graphID<<std::endl;
			Package *dPackage = rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dcpy[dPartner].recycleQueue.dequeue(true);

			//for(int j=0; j<sPackage->bufSize; j++)
			//dPackage->databuf[j] = sPackage->databuf[j];        //copy data------------------------------???????????????

			std::memcpy(dPackage->databuf.local(), sPackage->databuf.local(), dPackage->bufSize * sizeof(double));

			rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dcpy[dPartner].pQueue.enqueue(dPackage,true);
			if(rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dcpy[dPartner].pQueue.queueSize(true)==1)
			    rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.firingRuleCnt++;
			pthread_mutex_unlock(&(rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dLock));
			rg->lMap[f]->l_con.scpy[i].recycleQueue.enqueue(sPackage,true);
			pthread_mutex_unlock(&(rg->lMap[f]->l_con.sLock));
		    }
		}
		//pthread_mutex_unlock(&(rg->lMap[f]->l_con.sLock));
	    }// if(!lock succeedded)
	    if(perilla::LAZY_PUSH)
	    {
		//  Implemetation deffered. Currently not required
	    }	    
	}// if(tg==fg)
    }// for(f<numfabs)    
}//serviceLocalRequests

void Perilla::serviceRemoteRequests(RegionGraph* rg, int graphID, int nGraphs)
{
    bool nextsReq, nextrReq;
    int np = ParallelDescriptor::NProcs();
    int myProc = ParallelDescriptor::MyProc();
    int numfabs = rg->rMap.size();

    for(int f=0; f<numfabs; f++)
    {
	//int lockSucceeded = pthread_mutex_trylock(&(rg->rMap[f]->r_con.rcvLock));
	//if(lockSucceeded != 0)
	{
	    //if(pthread_mutex_trylock(&(rg->lMap[f]->r_con.rcvLock)) != 0)
	    {
		for(int i=0; i<rg->lMap[f]->r_con.nrcv; i++)
		{
		    if(rg->rMap[f]->r_con.rcv[i].pQueue.queueSize(true) == 0) //!no message has been received or all received messages have been claimed
			nextrReq = true;
		    else
		    {
			Package *rearPackage = rg->rMap[f]->r_con.rcv[i].pQueue.getRear(true);//!CHECK THIS POINT LATER
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

			rPackage->request = new future<>;
			rPackage->tag = tag;
			rg->lMap[f]->r_con.rcv[i].pQueue.enqueue(rPackage,true);   //!this is not done yet
			rg->rMap[f]->r_con.rcv[i].pQueue.enqueue(rMetaPackage,true);   //!this is not done yet
			pthread_mutex_lock(&(rMsgMap.lock));
			rMsgMap.map[rg->rMap[f]->r_con.rcv[i].pr][tag].push_back(rPackage);
                        rMsgMap.size++;
			pthread_mutex_unlock(&(rMsgMap.lock));
                        pthread_mutex_unlock(&(rg->lMap[f]->r_con.rcvLock));
                        pthread_mutex_unlock(&(rg->rMap[f]->r_con.rcvLock));
		    }
		}
		//pthread_mutex_unlock(&(rg->lMap[f]->r_con.rcvLock));
	    }// if(omp_test_lock)
	    //pthread_mutex_unlock(&(rg->rMap[f]->r_con.rcvLock));
	}// if(lockSucceeded)
    }// for(f<numfabs)

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
		    int ns = rg->sMap[f]->r_con.snd[i].ns;
		    int nd = rg->sMap[f]->r_con.snd[i].nd;
		    int r_gid = rg->sMap[f]->r_con.snd[i].r_gid;
		    int r_grids = rg->sMap[f]->r_con.snd[i].r_grids;
		    int tag = Perilla::myTagMap[r_gid][nd][ns][rg->sMap[f]->r_con.snd[i].sz];
		    int src= upcxx::rank_me();
		    //register send request so that the receiver can send back confirmation upon pull completion
                    sPackage->completed = false;
                    pthread_mutex_lock(&(sMsgMap.lock));
                    sMsgMap.map[rg->sMap[f]->r_con.snd[i].pr][tag].push_back(sPackage);
                    sMsgMap.size++;
                    pthread_mutex_unlock(&(sMsgMap.lock));
		    int size= rg->sMap[f]->r_con.snd[i].sz;
		    upcxx::global_ptr<double> sbuf= sPackage->databuf; //static_cast<upcxx::global_ptr<double> >((double*)sPackage->databuf);
		    int dst= rg->sMap[f]->r_con.snd[i].pr;
		    upcxx::rpc(dst, 
		       [=](){
			  //at destination rank, look up recv buffer and pull remote data and store data in the buffer
			  bool posted_recv=false;
			  double* localbuf= NULL;
			  pthread_mutex_lock(&(rMsgMap.lock));
			  if(rMsgMap.map.find(src) != rMsgMap.map.end()){
			      if(rMsgMap.map[src].find(tag) != rMsgMap.map[src].end())
				 if(rMsgMap.map[src][tag].size() >0){
				     posted_recv=true;
				     localbuf= (rMsgMap.map[src][tag].front()->databuf).local(); //(double*) static_cast<upcxx::global_ptr<void> > (rMsgMap.map[src][tag].front()->databuf).local(); 
			    	     *(rMsgMap.map[src][tag].front()->request)= upcxx::rget(sbuf, localbuf, size);
			             rMsgMap.map[src][tag].pop_front();
                                     rMsgMap.size--;
		                 }
			  } 
                          pthread_mutex_unlock(&(rMsgMap.lock));
			  //save pull request for later when recv buffer is posted 
			  if(posted_recv==false){
			      getReq_t *req= new getReq_t(src, tag, sbuf, size);
			      pendingGetList.add(req);
			  }
		       }
		    );
		}
	    }
	} // for(i<nsnd)
    } // for(f<numfabs)

    pendingGetList.process();

    //!now we test if send and receive requests have been serviced
    for(int f=0; f<numfabs; f++)
    {
	//!receive requests
	for(int i=0; i<rg->rMap[f]->r_con.nrcv; i++)
	{
	    if(rg->rMap[f]->r_con.rcv[i].pQueue.queueSize(true) > 0) //!all messages before rear have completed
	    {
		//if(pthread_mutex_trylock(&(rg->lMap[f]->r_con.rcvLock)) != 0) // 0-Fail, otherwise-Succeed
		{
		    Package *rearPackage =  rg->lMap[f]->r_con.rcv[i].pQueue.getRear(true);
		    if(!rearPackage->completed)
		    {
			bool flag = false;
			int ret_flag;
			if(rearPackage->request->ready())
			{
                            pthread_mutex_lock(&(rg->lMap[f]->r_con.rcvLock));
                              int ns =  rg->lMap[f]->r_con.rcv[i].ns;
                              int nd =  rg->lMap[f]->r_con.rcv[i].nd;
                              int lnd =  rg->lMap[f]->r_con.rcv[i].lnd;
                              int r_grids =  rg->lMap[f]->r_con.rcv[i].r_grids;
                              int tag = rearPackage->tag;
                              //int tag = tagMap[ rg->lMap[f]->r_con.rcv[i].pr][graphID][nd][ns][ rg->lMap[f]->r_con.rcv[i].sz];
                              int dst = upcxx::rank_me();
                              int src=  rg->lMap[f]->r_con.rcv[i].pr;
                                upcxx::rpc(src,
                                    [=](){
                                        pthread_mutex_lock(&(sMsgMap.lock));
                                        sMsgMap.map[dst][tag].front()->completed=true;
                                        sMsgMap.map[dst][tag].pop_front();
		                        sMsgMap.size--;
                                        pthread_mutex_unlock(&(sMsgMap.lock));
                                    }
                                );

			    delete rearPackage->request;
			    rearPackage->completeRequest();
			    rg->lMap[f]->r_con.rcv[i].pQueue.getRear()->completeRequest();
			    if(rg->rMap[f]->r_con.rcv[i].pQueue.queueSize(true) == 1)
				rg->lMap[f]->r_con.firingRuleCnt++;
                            pthread_mutex_unlock(&(rg->lMap[f]->r_con.rcvLock));
			}
		    }
		    //pthread_mutex_unlock(&(rg->lMap[f]->r_con.rcvLock));
		} // if(omp_test_lock)
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
		if(frontPackage->served) //!latest receive request has NOT been completed
		{
		    bool flag = false;
		    int ret_flag;
		    if(frontPackage->request==0)
		    {
			pthread_mutex_lock(&(rg->sMap[f]->r_con.sndLock));
			frontPackage = rg->sMap[f]->r_con.snd[i].pQueue.dequeue(true);
			frontPackage->completed = false;
			frontPackage->served = false;
			frontPackage->request = 0;
			frontPackage->tag = 0;
			frontPackage->notified = false;
			rg->sMap[f]->r_con.snd[i].recycleQueue.enqueue(frontPackage,true);
			pthread_mutex_unlock(&(rg->sMap[f]->r_con.sndLock));

			pthread_mutex_lock(&(rg->lMap[f]->r_con.sndLock));
			frontPackage = rg->lMap[f]->r_con.snd[i].pQueue.dequeue(true);
			frontPackage->completed = false;
			frontPackage->served = false;
			frontPackage->request = 0;
                        frontPackage->notified = false;
			frontPackage->tag = 0;
			rg->lMap[f]->r_con.snd[i].recycleQueue.enqueue(frontPackage,true);
			pthread_mutex_unlock(&(rg->lMap[f]->r_con.sndLock));			
		    }
		}
	    } // if(queueSize > 0)
	} // for(i<nsnd)
    }// for(f<numfabs)
}// serviceRemoteRequests

void Perilla::serviceRemoteRequests(RegionGraph* rg)
{
    serviceRemoteRequests(rg,0,1);
}

void Perilla::serviceSingleGraphComm(RegionGraph* graph, int tid)
{
    int np = ParallelDescriptor::NProcs();
    int tg = WorkerThread::perilla_wid();
    while(true)
    {
	if(graph->totalFinishes < perilla::NUM_THREAD_TEAMS)
	{
	    serviceLocalRequests(graph, tg);
	    if((np>1) & (tg==0))
		serviceRemoteRequests(graph);
	}
	else
	{
	    if(tg==0)
	    {
		while(graph->totalFinishes < perilla::NUM_THREAD_TEAMS)
		{
		}		
		//call parallel_barrier()  ---????????
		ParallelDescriptor::Barrier("serviceSingleGraph-1");
		graph->graphTeardown();
		graph->workerTeardown();
		//call parallel_barrier() ------?????????
		ParallelDescriptor::Barrier("serviceSingleGraph-2");
	    }
	    break;
	}
    } // while(true)
} //serviceSingleGraphComm

void Perilla::serviceMultipleGraphComm(RegionGraph graphArray[], int nGraphs, bool cpyAcross, int tid)
{
    int tg = WorkerThread::perilla_wid();
    int np = ParallelDescriptor::NProcs();    
    int graphFinishCnt = 0;
    while(true)
    {
	for(int g=0; g<nGraphs; g++)
	{
	    if(graphArray[g].totalFinishes < perilla::NUM_THREAD_TEAMS)
	    {
		serviceLocalRequests(&graphArray[g], tg);
		//if(cpyAcross)
		//serviceLocalGridCopyRequests(graphArray,g,tg);
		if(np > 1)
		    if(tg==0)
		    {
			serviceRemoteRequests(&graphArray[g],g,nGraphs);
			//if(cpyAcross)
			//serviceRemoteGridCopyRequests(graphArray,g,nGraphs,tg);
		    }
	    }
	}
	//!check if we have finished all the graph execution
	bool noMoreWork = true;
	for(int g=0; g<nGraphs; g++)
	    if(graphArray[g].totalFinishes < perilla::NUM_THREAD_TEAMS)
		noMoreWork = false;
	if(noMoreWork)
	    break;	
    } // while(true)
    if(tg==0)
	for(int g=0; g<nGraphs; g++)
	{
	    ParallelDescriptor::Barrier("serviceMultipleGraph-1");
	    graphArray[g].graphTeardown();
	    graphArray[g].workerTeardown();
	    ParallelDescriptor::Barrier("serviceMultipleGraph-2");
	}
} // serviceMultipleGraphComm

#if 1
void Perilla::serviceMultipleGraphCommDynamic(std::vector<RegionGraph*> graphArray, bool cpyAcross, int tid)
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

    //while(true)
    {	
	for(int g=0; g<graphArray.size(); g++)
	{
	    nGraphs = graphArray.size();
	    {
		serviceLocalRequests(graphArray[g], tg);
		if(cpyAcross)
		{
		    serviceLocalGridCopyRequests(graphArray,g,tg);
		}
                if(np > 1)//if(tg==0)
                {
		    serviceRemoteRequests(graphArray[g],g,nGraphs);
		    if(cpyAcross)
		        if(tg==0)serviceRemoteGridCopyRequests(graphArray,g,nGraphs,tg);
                }
	    }
	}
	//if( Perilla::numTeamsFinished == perilla::NUM_THREAD_TEAMS)
	//{
	//    if(doublechecked) // double check if there are still something to send
	//	break;
	//    else
	//	doublechecked = true;
	//}
	numloops++;
	avgltime += ltime;
    } // while(true)

    //nGraphs = graphArray.size();
    //for(int g=0; g<nGraphs; g++)
    //{
	//ParallelDescriptor::Barrier("serviceMultipleGraph-1");
	//graphArray[g]->graphTeardown(tg);
	//graphArray[g]->workerTeardown(tg);
	//ParallelDescriptor::Barrier("serviceMultipleGraph-2");
    //}
} // serviceMultipleGraphCommDynamic
#endif

void Perilla::serviceMultipleGraphComm(RegionGraph graphArray[], int nGraphs, int tid)
{
    serviceMultipleGraphComm(graphArray,nGraphs,false,tid);
} // serviceMultipleGraphComm

void Perilla::flattenGraphHierarchy(std::vector<std::vector<RegionGraph*> > graphArrayHierarchy, std::vector<RegionGraph*> &graphArray){
    graphArray.clear();
    int gCnt=0;
    for(int l=0; l<graphArrayHierarchy.size(); l++) gCnt+= graphArrayHierarchy[l].size();
    for(int l=0; l<graphArrayHierarchy.size(); l++)
        for(int g=0; g<graphArrayHierarchy[l].size(); g++)
            graphArray.push_back(graphArrayHierarchy[l][g]);
}


void Perilla::fillBoundaryPush(RegionGraph* graph, MultiFab* mf, int f)
{

    int nComp = mf->nComp();
    int tg= WorkerThread::perilla_wid();
    int ntid = WorkerThread::perilla_wtid();

    //if(graph->graphID == 1 && f == 1)
    //std::cout << "fillBPush for gID 1 f 1 ntid "<< ntid <<std::endl;

    if(perilla::LAZY_PUSH)
    { }
    else
    {
	if(ntid == 0)
	    pthread_mutex_lock(&(graph->lMap[f]->l_con.sLock));
	graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads

	if(perilla::PACKING_FINEGRAIN)
	{}
	else
	{
	    for(int i=0; i<graph->lMap[f]->l_con.nscpy; i++)
		if( (i%(perilla::NUM_THREADS_PER_TEAM-1)) == ntid)
		{

		    //if(graph->graphID == 1 && graph->lMap[f]->l_con.scpy[i].nd == 1)
		    //std::cout << "fillBPush for gID 1 nd 1 pQenQ f " << f << " i " << i <<std::endl;
		    Package *sPackage = graph->lMap[f]->l_con.scpy[i].recycleQueue.getFront(true);
		    mf->m_fabs_v[f]->copyToMem(graph->lMap[f]->l_con.scpy[i].sbx,0,nComp,sPackage->databuf.local());

		    //for(int d=0; d<sPackage->bufSize; d++)
			//if(sPackage->databuf[d] == 0)
			//{
			    //std::cout<< "in fbPush Sending 0 from f "<< f <<std::endl;
			    //BL_ASSERT(sPackage->databuf[d] != 0);
			//}
		    //if(graph->lMap[f]->l_con.scpy[i].sbx.smallEnd() == graph->lMap[f]->l_con.scpy[i].sbx.bigEnd())
		    //if(graph->lMap[f]->l_con.scpy[i].sbx.smallEnd(0)==7 && graph->lMap[f]->l_con.scpy[i].sbx.smallEnd(1)==7 && graph->lMap[f]->l_con.scpy[i].sbx.smallEnd(2)==4)
		    //  std::cout<< "Corner Push for f "<< f << " data0 " <<sPackage->databuf[0]<< " size " <<sPackage->bufSize << " se "<< graph->lMap[f]->l_con.scpy[i].sbx.smallEnd() <<std::endl;

		}	    
	} // if(PACKING_FINEGRAIN) - else

	graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads

	if(ntid==0)
	{
	    //if(graph->graphID == 1 && f == 1)
	    //std::cout << "fillBPush for gID 1 f 1 pQ enQ" <<std::endl;
	    for(int i=0; i<graph->lMap[f]->l_con.nscpy; i++)
	    {
		//if(graph->graphID == 1 && graph->lMap[f]->l_con.scpy[i].nd == 1)
		//std::cout << "fillBPush for gID 1 nd 1 pQ enQ from f "<< f <<std::endl;
		graph->lMap[f]->l_con.scpy[i].pQueue.enqueue( graph->lMap[f]->l_con.scpy[i].recycleQueue.dequeue(true),true );
	    }
	    pthread_mutex_unlock(&(graph->lMap[f]->l_con.sLock));
	}
	graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads
    } // if(LAZY_PUSH) - else

    int np = ParallelDescriptor::NProcs();
    if (np==1) return;

    if(ntid==0)
	pthread_mutex_lock(&(graph->lMap[f]->r_con.sndLock));
    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads

    for(int i=0; i<graph->lMap[f]->r_con.nsnd; i++)
	if((i%(perilla::NUM_THREADS_PER_TEAM-1))==ntid)
	{
	    //std::cout << "RQS " << graph->lMap[f]->r_con.snd[i].recycleQueue.queueSize() << std::endl;

	    Package *sndPackage = graph->lMap[f]->r_con.snd[i].recycleQueue.dequeue(true);	  
	    mf->m_fabs_v[f]->copyToMem(graph->lMap[f]->r_con.snd[i].sbx,0,nComp,sndPackage->databuf.local());
	    sndPackage->notified = false;
	    graph->lMap[f]->r_con.snd[i].pQueue.enqueue( sndPackage,true );
	    //!the local message handler will detect the change and notify the remote message handler =>read access
	    //!the remote message handler first modifies the front item of this queue, then it push this item back to the message pool
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

void Perilla::fillBoundaryPull(RegionGraph* graph, MultiFab* mf, int f, bool singleT)
{

    int nComp = mf->nComp();
    int tg= WorkerThread::perilla_wid();
    int ntid = WorkerThread::perilla_wtid();

    if(ntid==0)
	pthread_mutex_lock(&(graph->lMap[f]->l_con.dLock));
    if(!singleT)
        graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads    

    if(perilla::LAZY_PUSH)
    { }
    else
    {
	if(perilla::UNPACKING_FINEGRAIN)
	{}
	else
	{
	    for(int i=0; i<graph->lMap[f]->l_con.ndcpy; i++)
		if( (i%(perilla::NUM_THREADS_PER_TEAM-1)) == ntid)
		{
		    Package *dPackage = graph->lMap[f]->l_con.dcpy[i].pQueue.getFront(true);

		    //for(int d=0; d<dPackage->bufSize; d++)
			//if(dPackage->databuf[d] == 0)
			//{
			    //std::cout<< "in fbPull Reciving 0 for f "<< f <<std::endl;
			    //BL_ASSERT(dPackage->databuf[d] != 0);
			//}
		    /*
		       if(f==0)
		    //if(graph->lMap[f]->l_con.dcpy[i].dbx.smallEnd() == graph->lMap[f]->l_con.dcpy[i].dbx.bigEnd())
		    //if(graph->lMap[f]->l_con.dcpy[i].dbx.smallEnd(0)==-1 && graph->lMap[f]->l_con.dcpy[i].dbx.smallEnd(1)==-1 && graph->lMap[f]->l_con.dcpy[i].dbx.smallEnd(2)==4)
		    std::cout<< "Corner Pull for f "<< f << " data0 " <<dPackage->databuf[0]<< " size " <<dPackage->bufSize <<" se " <<graph->lMap[f]->l_con.dcpy[i].dbx.smallEnd()<<std::endl;
		     */
		    mf->m_fabs_v[f]->copyFromMem(graph->lMap[f]->l_con.dcpy[i].dbx,0,nComp,dPackage->databuf.local());		  
		}
	} // if(UNPACKING_FINEGRAIN) - else
    } // if(LAZY_PUSH) - else

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
	    mf->m_fabs_v[f]->copyFromMem(graph->lMap[f]->r_con.rcv[i].dbx,0,nComp,rcvPackage->databuf.local());
	    rcvPackage->completed = false;
	    rcvPackage->notified = false;
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

  void Perilla::fillBoundaryPull(amrex::RGIter& rgi, RegionGraph* rg, amrex::MultiFab& mf, bool singleT)
  {
    if(rgi.currentItr != 1)
      return;

    int f = rgi.currentRegion;
    fillBoundaryPull(rg, &mf, f, singleT);
  }

  void Perilla::fillBoundaryPull(amrex::RGIter& rgi, amrex::MultiFab& mf, bool singleT)
  {
    if(rgi.currentItr != 1)
      return;

    int f = rgi.currentRegion;
    fillBoundaryPull(rgi.itrGraph, &mf, f, singleT);
  }


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if 0
void Perilla::multifabExtractCopyAssoc(RegionGraph* gDst, RegionGraph* gSrc, const MultiFab& mfDst, const MultiFab& mfSrc, int nc, int ng, int ngSrc, const Periodicity& period)
{
    //    MultiFab* mfSrc = gSrc->assocMF;
    //    MultiFab* mfDst = gDst->assocMF;
    int myProc = ParallelDescriptor::MyProc();
    int np = ParallelDescriptor::NProcs();

    try{

	if(true)//if(!(*mfSrc == *mfDst))
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

	    if(gSrc->numTasks != mfSrc.IndexArray().size())
		std::cout<< "before " <<gSrc->numTasks << " now " <<mfSrc.size() << " at gID " << gSrc->graphID << std::endl;

	    gSrc->numFabs = mfSrc.size();
	    gDst->numFabs = mfDst.size();

	    gSrc->numTasks = mfSrc.IndexArray().size();
	    gDst->numTasks = mfDst.IndexArray().size();

	    int nfabsSrc = mfSrc.IndexArray().size();
	    int nfabsDst = mfDst.IndexArray().size();

	    const FabArrayBase::CPC& TheCPC = mfDst.getCPC(IntVect(ng), mfSrc, IntVect(ngSrc), period);

	    const int nloc_cpAsc = TheCPC.m_LocTags->size();
	    const int nsnds_cpAsc = TheCPC.m_SndTags->size();
	    const int nrcvs_cpAsc = TheCPC.m_RcvTags->size();

	    Vector<const FabArrayBase::CopyComTagsContainer*> send_cctc;
	    Vector<int> send_pr;
	    send_cctc.reserve(nsnds_cpAsc);

	    for (FabArrayBase::MapOfCopyComTagContainers::const_iterator m_it = TheCPC.m_SndTags->begin(),
		    m_End = TheCPC.m_SndTags->end();
		    m_it != m_End;
		    ++m_it)
	    {
		if(m_it->first != myProc)      // Not destined to me.
		{
		    send_pr.push_back(m_it->first);
		    send_cctc.push_back(&(m_it->second));
		}
	    }

	    //  std::cout<< "Loop 1" <<std::endl;

	    Vector<const FabArrayBase::CopyComTagsContainer*> recv_cctc;
	    Vector<int> recv_pr;
	    recv_cctc.reserve(nrcvs_cpAsc);

	    for (FabArrayBase::MapOfCopyComTagContainers::const_iterator m_it = TheCPC.m_RcvTags->begin(),
		    m_End = TheCPC.m_RcvTags->end();
		    m_it != m_End;
		    ++m_it)
	    {
		if(m_it->first != myProc)      // I am not the source for this receipt
		{
		    recv_pr.push_back(m_it->first);
		    recv_cctc.push_back(&(m_it->second));
		}
	    }

	    //std::cout<< "Before parallel at gID " << gDst->graphID << " numTask " << gDst->numTasks << " numFabs " << gDst->numFabs <<std::endl;      

//#pragma omp parallel shared(gSrc, gDst, mfSrc, mfDst, nfabsSrc, nfabsDst)
	    {
//		int tid = omp_get_thread_num();//perilla::tid();//omp_get_thread_num();   
//		int tg = tid/perilla::NUM_THREADS_PER_TEAM;//perilla::wid();//WorkerThread::perilla_wid();
//		int nt= tid%perilla::NUM_THREADS_PER_TEAM;
//		int fg;

		for(int f=0; f<nfabsSrc; f++)
		{
//		    if(nt==0)
//			if(WorkerThread::isMyRegion(tg,f))// && WorkerThread::perilla_isMasterWorkerThread())
			{
			    //if(gDst->graphID > 25)
			    //std::cout<< "Inside parallel Generating Send at tid " << tid << " f " << f << " gID " << gDst->graphID <<std::endl;   

			    FabCopyAssoc *cpSrc;
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
				const FabArrayBase::CopyComTag& tag = (*TheCPC.m_LocTags)[i];
				//if(f == tag.srcIndex)
				if(mfSrc.IndexArray()[f] == tag.srcIndex)
				    cpSrc->l_con.nscpy++;
			    }
			    cpSrc->l_con.scpy = new LocalCopyDescriptor[cpSrc->l_con.nscpy];
			    int scnt = 0;
			    //if(gDst->graphID == 4 && tag.dstIndex == 60 )
			    //std::cout<< "Inside parallel Generating Local Copy send at tid " << tid << " f " << f << " gID " << gDst->graphID <<std::endl;        

			    for(int i=0; i<nloc_cpAsc; i++)
			    {
				const FabArrayBase::CopyComTag& tag = (*TheCPC.m_LocTags)[i];
				//if(f == tag.srcIndex)
				if(mfSrc.IndexArray()[f] == tag.srcIndex)
				{

				    //if(gDst->graphID == 4 && (tag.dstIndex == 60 || tag.dstIndex == 59) )
				    //std::cout <<"myP " <<myProc<< " Added in S LocDep nd " << tag.dstIndex << " ns "<< tag.srcIndex << " f " << f << " i "<< scnt << " tg " <<tg << std::endl;

				    cpSrc->l_con.scpy[scnt].ns = mfSrc.localindex(tag.srcIndex);
				    cpSrc->l_con.scpy[scnt].nd = mfDst.localindex(tag.dstIndex);
				    cpSrc->l_con.scpy[scnt].sbx = tag.sbox;
				    cpSrc->l_con.scpy[scnt].dbx = tag.dbox;
				    int psize = tag.sbox.numPts() * mfSrc.nComp(); //---------------------------------------------------------------????????????????
				    //std::cout<< " gSrc ID "<< gSrc->graphID << " f "<<f<< " sndPkgsize " << psize <<std::endl;
				    for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
				    {
					Package *tmpPkg = new Package(psize);
			                void* local_ptr= tmpPkg->databuf.local(); //(static_cast<upcxx::global_ptr<void> >(tmpPkg->databuf)).local();
					assert(local_ptr!=0);
					for(int j=0; j<psize; j++){
					    //tmpPkg->databuf[j] = 0;
					    ((double*)local_ptr)[j]=0;
				        }
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

					    //if(gDst->graphID == 17 && (it->srcIndex == 1198 || it->srcIndex == 1198 || it->srcIndex == 978 || it->srcIndex == 978))
					    //std::cout <<"myP " <<myProc<< " Added in S Dep r " << it->dstIndex << " s "<< it->srcIndex << " f " << f << " i "<< scnt << " tg " <<tg <<" nsnd " << gSrc<< std::endl;

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
			    		        void* local_ptr= tmpPkg->databuf.local();//(static_cast<upcxx::global_ptr<void> >(tmpPkg->databuf)).local();
						for(int j=0; j<psize; j++){
						    //tmpPkg->databuf[j] = 0;
						    ((double*)local_ptr)[j]=0;
						}
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
//#pragma omp barrier
		    //        std::cout<< "Barrier 1" <<std::endl;            
		    if(np > 1)
		    {
			//if(WorkerThread::perilla_isMasterWorkerThread() && tg==0)
			//if(tid==0)
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
		//        std::cout<< "Barrier 2 " <<" tid "<<tid<<std::endl;     
//#pragma omp barrier               
		for(int f=0; f<nfabsDst; f++)
		{
		    //if(nt==0)
			//if(WorkerThread::isMyRegion(tg,f))// && perilla::isMasterWorkerThread())
			{
			    //        std::cout <<"tid: "<< tid << " f: "<< f << " is master "<<WorkerThread::isTeamMasterThread(tid) << " is my region "<<WorkerThread::isMyRegion(tg,f)<<std::endl;                 

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
				const FabArrayBase::CopyComTag& tag = (*TheCPC.m_LocTags)[i];
				//if(f == tag.dstIndex)
				if(mfDst.IndexArray()[f] == tag.dstIndex)
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
				const FabArrayBase::CopyComTag& tag = (*TheCPC.m_LocTags)[i];
				//if(f == tag.dstIndex)
				if(mfDst.IndexArray()[f] == tag.dstIndex)
				{

				    //if(gDst->graphID == 4 && (tag.dstIndex == 60 || tag.dstIndex == 59))
				    //std::cout<< "dcpy tid " << tid << " f " << f << " i " << i << " dcnt " << dcnt << " ns "<<tag.srcIndex << " nd "<<tag.dstIndex << " lo " << tag.dbox.smallEnd() << " hi " << tag.dbox.bigEnd() <<std::endl;   

				    cpDst->l_con.dcpy[dcnt].ns = mfSrc.localindex(tag.srcIndex);
				    cpDst->l_con.dcpy[dcnt].nd = mfDst.localindex(tag.dstIndex);
				    cpDst->l_con.dcpy[dcnt].sbx = tag.sbox;
				    cpDst->l_con.dcpy[dcnt].dbx = tag.dbox;

				    // if(gDst->graphID > 25 && f == 633)
				    //std::cout<< " Generating Package tid " << tid << " i " << i <<std::endl;      

				    int psize = tag.dbox.numPts() * mfSrc.nComp(); //---------------------------------------------------------------????????????????
				    cpDst->l_con.dcpy[dcnt].sz = psize;

				    if(!gDst->isDepGraph)
				    {
					for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
					{
					    Package *tmpPkg = new  Package(psize);
			    		    void* local_ptr= tmpPkg->databuf.local();//(static_cast<upcxx::global_ptr<void> >(tmpPkg->databuf)).local();

					    // if(tmpPkg == nullptr)
					    //std::cout<<"Found the culprit tid " << tid << " f " << f << " i " << i << std::endl;

					    for(int j=0; j<psize; j++){
						//tmpPkg->databuf[j] = 0;
						((double*)local_ptr)[j]=0;
					    }
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
			    		void* local_ptr= tmpPkg->databuf.local(); //(static_cast<upcxx::global_ptr<void> >(tmpPkg->databuf)).local();
					for(int j=0; j<psize; j++){
					    //tmpPkg->databuf[j] = 0;
					    ((double*)local_ptr)[j]=0;
				        }
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
			    			    void* local_ptr= tmpPkg->databuf.local(); //(static_cast<upcxx::global_ptr<void> >(tmpPkg->databuf)).local();
						    for(int j=0; j<psize; j++){
							//tmpPkg->databuf[j] = 0;
							((double*)local_ptr)[j]=0;
						    }
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
			    		    void* local_ptr= tmpPkg->databuf.local(); //(static_cast<upcxx::global_ptr<void> >(tmpPkg->databuf)).local();
					    for(int j=0; j<psize; j++){
						//tmpPkg->databuf[j] = 0;
						((double*)local_ptr)[j]=0;
					    }
					    cpdDst->r_con.rcv[i].pQueue.enqueue(tmpPkg);
					}
					for(int p=0; p<perilla::NUM_PREGENERATED_PACKAGES; p++)
					    cpdDst->r_con.rcv[i].recycleQueue.enqueue(cpdDst->r_con.rcv[i].pQueue.dequeue());
				    }
				}
			    } // if(np > 1)
			}// if(fg==tg)

//#pragma omp barrier               
		    if(np > 1)
		    {
			//if(WorkerThread::perilla_isMasterWorkerThread() && tg==0)
			//if(tid==0)
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
//#pragma omp barrier
	    for(int f=0; f<nfabsSrc; f++)
	    {
		//if(nt==0)
		    //if(WorkerThread::isMyRegion(tg,f))// && WorkerThread::perilla_isMasterWorkerThread())
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
//#pragma omp barrier
	    for(int f=0; f<nfabsDst; f++)
	    {
		//if(nt==0)
		    //if(WorkerThread::isMyRegion(tg,f))// && WorkerThread::perilla_isMasterWorkerThread())
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

} // multifabExtractCopyAssoc
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
		int tid = omp_get_thread_num();//perilla::tid();//omp_get_thread_num();	  
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


void Perilla::multifabExtractCopyAssoc(RegionGraph* gDst, RegionGraph* gSrc, const MultiFab& mfDst, const MultiFab& mfSrc, const Periodicity& period) 
{
    multifabExtractCopyAssoc(gDst, gSrc, mfDst, mfSrc, 1 /*component*/, 0/*ghost*/, 0/*src ghost*/, period);
}

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

	//if(srcGraph->graphID==18 && f ==316 && ntid == 0)
	//std::cout << "srgG chk see " << srcGraph << " " <<myProc <<std::endl;

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
		mfSrc->m_fabs_v[f]->copyToMem(cpSrc->l_con.scpy[i].sbx,srccomp,nc,sndPackage->databuf.local());
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
		    mfSrc->m_fabs_v[f]->copyToMem(cpSrc->l_con.scpy[i].sbx,srccomp,nc,sndPackage->databuf.local());
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
		mfSrc->m_fabs_v[f]->copyToMem(cpSrc->r_con.snd[i].sbx,srccomp,nc,sndPackage->databuf.local());
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
	    srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);

	    for(int i=0; i<cpSrc->r_con.nsnd; i++)
		if((i%(perilla::NUM_THREADS_PER_TEAM-1)) == ntid)
		{

		    // if(myProc==4 && srcGraph->graphID==2 && (f ==0 || f ==2))
		    //std::cout << " Pushing 2 316 164"<<std::endl;

		    Package* sndPackage = cpSrc->r_con.snd[i].recycleQueue.dequeue(true);
		    mfSrc->m_fabs_v[f]->copyToMem(cpSrc->r_con.snd[i].sbx,srccomp,nc,sndPackage->databuf.local());
		    sndPackage->notified = false;
		    sndPackage->notified = false;
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
		mfDst->m_fabs_v[f]->copyFromMem(cpDst->l_con.dcpy[i].dbx,dstcomp,nc,rcvPackage->databuf.local());
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
		    mfDst->m_fabs_v[f]->copyFromMem(cpDst->l_con.dcpy[i].dbx,dstcomp,nc,rcvPackage->databuf.local());
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
	    //pthread_mutex_lock(&(destGraph->rCopyMapHead->map[f]->r_con.rcvLock));
	    pthread_mutex_lock(&(cpDst->r_con.rcvLock));
	    for(int i=0; i<cpDst->r_con.nrcv; i++)
	    {
		///*
		//Package *rcvMetaPackage = destGraph->rCopyMapHead->map[f]->r_con.rcv[i].pQueue.dequeue(true);
                Package* rcvPackage = cpDst->r_con.rcv[i].pQueue.dequeue(true);
		mfDst->m_fabs_v[f]->copyFromMem(cpDst->r_con.rcv[i].dbx,dstcomp,nc,rcvPackage->databuf.local());	      
		rcvPackage->completed = false;
		rcvPackage->served = false;
		rcvPackage->request = 0;	 
                cpDst->r_con.rcv[i].recycleQueue.enqueue(rcvPackage, true); 
		//destGraph->rCopyMapHead->map[f]->r_con.rcv[i].recycleQueue.enqueue(rcvMetaPackage,true);

		//Package* rcvPackage = cpDst->r_con.rcv[i].pQueue.dequeue(true);                               // corrected from recycleQ to pQ
		//rcvPackage->notified = false;
		//rcvPackage->completed = false;
		//cpDst->r_con.rcv[i].recycleQueue.enqueue(rcvPackage,true);                         // corrected from pQ to recycleQ	      
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
	    //pthread_mutex_unlock(&(destGraph->rCopyMapHead->map[f]->r_con.rcvLock));

	}
	else
	{	
	    if(ntid==0)
	    {
		//pthread_mutex_lock(&(destGraph->rCopyMapHead->map[f]->r_con.rcvLock));
		pthread_mutex_lock(&(cpDst->r_con.rcvLock));
	    }
	    destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);

	    for(int i=0; i<cpDst->r_con.nrcv; i++)
		if((i%(perilla::NUM_THREADS_PER_TEAM-1)) == ntid)
		{
		    ///*
		    //Package *rcvMetaPackage = destGraph->rCopyMapHead->map[f]->r_con.rcv[i].pQueue.dequeue(true);
                    Package* rcvPackage = cpDst->r_con.rcv[i].pQueue.dequeue(true); 
		    mfDst->m_fabs_v[f]->copyFromMem(cpDst->r_con.rcv[i].dbx,dstcomp,nc,rcvPackage->databuf.local());	      
		    rcvPackage->completed = false;
		    rcvPackage->served = false;
		    rcvPackage->request = 0;	  
                    cpDst->r_con.rcv[i].recycleQueue.enqueue(rcvPackage, true);
		    //destGraph->rCopyMapHead->map[f]->r_con.rcv[i].recycleQueue.enqueue(rcvMetaPackage,true);

		    //Package* rcvPackage = cpDst->r_con.rcv[i].pQueue.dequeue(true);                               // corrected from recycleQ to pQ
		    //rcvPackage->notified = false;
		    //rcvPackage->completed = false;
		    //cpDst->r_con.rcv[i].recycleQueue.enqueue(rcvPackage,true);                         // corrected from pQ to recycleQ	      
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
		//pthread_mutex_unlock(&(destGraph->rCopyMapHead->map[f]->r_con.rcvLock));
	    }
	    destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
	}
    } // if(!(*mfDst == *mfSrc))

} // multifabCopyPull

void Perilla::multifabCopyPull(RegionGraph* destGraph, RegionGraph* srcGraph, MultiFab* mfDst, MultiFab* mfSrc, int f, bool singleT) 
{
    multifabCopyPull(destGraph, srcGraph, mfDst, mfSrc, f, 1, 1, 1, 0, 0,singleT);
}

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
		//int lockSucceeded = pthread_mutex_trylock(&(cpSrc->l_con.sLock));
		//if(lockSucceeded != 0)
		{
		    for(int i=0; i<cpSrc->l_con.nscpy; i++)
		    {
			if(cpSrc->l_con.scpy[i].pQueue.queueSize()>0)
			{
                            pthread_mutex_lock(&(cpSrc->l_con.sLock));
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
			    /*
			       for(int j=0; j<dPackage->bufSize; j++)
			       {
			       dPackage->databuf[j] = sPackage->databuf[j];
			       }
			     */
			    std::memcpy(dPackage->databuf.local(), sPackage->databuf.local(), dPackage->bufSize * sizeof(double));
			    //std::swap(dPackage->databuf, sPackage->databuf);


			    cpDst->l_con.dcpy[dPartner].pQueue.enqueue(dPackage,true);
			    if(cpDst->l_con.dcpy[dPartner].pQueue.queueSize(true) == 1)
				cpDst->l_con.firingRuleCnt++;
			    pthread_mutex_unlock(&(cpDst->l_con.dLock));
			    cpSrc->l_con.scpy[i].recycleQueue.enqueue(sPackage,true);
                            pthread_mutex_unlock(&(cpSrc->l_con.sLock));
			}
		    } // for
		    //pthread_mutex_unlock(&(cpSrc->l_con.sLock));
		} // if(lockSucceeded)
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

    for(int f=0; f<numfabs; f++)
    {	
	FabCopyAssoc* cpDst = graphArray[g]->task[f]->cpAsc_dstHead;
	while(cpDst != 0)
	{
	    //if(pthread_mutex_trylock(&(graphArray[g]->rCopyMapHead->map[f]->r_con.rcvLock)) != 0)
	    {
		//if(pthread_mutex_trylock(&(cpDst->r_con.rcvLock)) != 0)
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
			    // Also check the recycle queue because when rear is completed it may cause unlimited recv posts
                            if(rearPackage->completed && cpDst->r_con.rcv[i].pQueue.queueSize(true) == 1) //!latest receive request has been completed
			    {
				nextrReq = true;
			    }
			    else //!expected message is still on the way
				nextrReq = false;
			}
			if(nextrReq) //!take a message from recycle pool and post a receive
			{
			    //!create a package to keep track of receive requests
		            pthread_mutex_lock(&(cpDst->r_con.rcvLock));
			    //Package *rMetaPackage = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].recycleQueue.dequeue(true);
			    //!extract a package from the recycle pool at the destination NUMA node to buffer incoming data
			    int ns = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].ns;
			    int nd = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].nd;
			    int lnd = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].lnd;
			    int r_grids = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].r_grids;
			    Package *rPackage = cpDst->r_con.rcv[i].recycleQueue.dequeue(true);
			    int tag = tagMap[graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pr][g][nd][ns][graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].sz];

                            rPackage->request = new future<>;
                            rPackage->completed=false;
			    rPackage->tag = tag;
			    cpDst->r_con.rcv[i].pQueue.enqueue(rPackage,true);   //!this is not done yet
			    //graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pQueue.enqueue(rMetaPackage,true);   //!this is not done yet	 

                            pthread_mutex_lock(&(rMsgMap.lock));
                            rMsgMap.map[graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pr][tag].push_back(rPackage);
                            rMsgMap.size++;
                            pthread_mutex_unlock(&(rMsgMap.lock));
		            pthread_mutex_unlock(&(cpDst->r_con.rcvLock));
			}						
		    } // for (i<i<cpDst->r_con.nrcv)
		    //pthread_mutex_unlock(&(cpDst->r_con.rcvLock));
		} // if(ga locked)
		//pthread_mutex_unlock(&(graphArray[g]->rCopyMapHead->map[f]->r_con.rcvLock));
	    } // if(mf locked)
	    cpDst = cpDst->next;
	} // while(cpDst != 0)	
    } // for(f<nfabs)

    for(int f=0; f<numfabs; f++)
    {	
	FabCopyAssoc* cpSrc = graphArray[g]->task[f]->cpAsc_srcHead;
	while(cpSrc != 0)
	{
	    for(int i=0; i<cpSrc->r_con.nsnd; i++)
	    {
		//if(graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].pQueue.queueSize(true) == 0) //!no message has been received or all received messages have been claimed	       	
                if(cpSrc->r_con.snd[i].pQueue.queueSize(true) == 0)
		    nextrReq = false;
		else
		    nextrReq = true;

		if(nextrReq) //!take a message from recycle pool and post a receive
		{
		    //Package *sMetaPackage = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].pQueue.getFront(true);
                    Package *sPackage = cpSrc->r_con.snd[i].pQueue.getFront(true);
		    if(!sPackage->served)
		    {		    
			sPackage->completed = false;
			sPackage->served = true;
                        //sMetaPackage->request = new future<>;
			int ns = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].ns;
			int nd = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].nd;
			int r_gid = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].r_gid;
			int r_grids = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].r_grids;
			//int tag = tagGen(ns, nd, r_gid-1, np*r_grids, nGraphs);
			int tag = Perilla::myTagMap[r_gid][nd][ns][graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].sz];
                        int src= upcxx::rank_me();
			int dst= graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].pr;
		        int size= graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].sz;

                        //sPackage->request = new future<>;
                        pthread_mutex_lock(&(sMsgMap.lock));
                        sMsgMap.map[dst][tag].push_back(sPackage);
                        sMsgMap.size++;
                        pthread_mutex_unlock(&(sMsgMap.lock));
                        upcxx::global_ptr<double> sbuf= sPackage->databuf; //static_cast<upcxx::global_ptr<double> >((double*)sPackage->databuf);

                        upcxx::rpc(dst,
                            [=](){
                                //at destination rank, look up recv buffer and pull remote data and store data in the buffer
                                bool posted_recv=false;
                                double* localbuf= NULL;
                                pthread_mutex_lock(&(rMsgMap.lock));
                                if(rMsgMap.map.find(src) != rMsgMap.map.end()){
                                    if(rMsgMap.map[src].find(tag) != rMsgMap.map[src].end())
                                         if(rMsgMap.map[src][tag].size() >0){
                                             posted_recv=true;
                                             localbuf= (rMsgMap.map[src][tag].front()->databuf).local();//(double*) (static_cast<upcxx::global_ptr<void> > (rMsgMap.map[src][tag].front()->databuf).local());
   	                                     rMsgMap.map[src][tag].front()->tag= tag;
					     if(localbuf){
                                                 *(rMsgMap.map[src][tag].front()->request)= upcxx::rget(sbuf, localbuf, size);
                                                 rMsgMap.map[src][tag].pop_front();
			                         rMsgMap.size--;
				             }
                                         }
                                }
                                pthread_mutex_unlock(&(rMsgMap.lock));
                                //save pull request for later when recv buffer is posted 
                                if(posted_recv==false){
                                    getReq_t *req= new getReq_t(src, tag, sbuf, size);
                                    pendingGetList.add(req);
                                }
                                //store send request to notify sender later upon completion
                                //sFutureMap[fu]= sMetaPackage->request;
                            }
                        );
		    } //served
		}//nextReq		
	    } // for (i<i<cpSrc->r_con.nsnd)	    
	    cpSrc = cpSrc->next;
	} // while(cpSrc != 0)	
    } // for(f<nfabs)

    pendingGetList.process();

    for(int f=0; f<numfabs; f++)
    {	
	FabCopyAssoc* cpDst = graphArray[g]->task[f]->cpAsc_dstHead;
	while(cpDst != 0)
	{
	    for(int i=0; i<cpDst->r_con.nrcv; i++)
	    {
		//if(graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pQueue.queueSize(true) > 0) //!all messages before rear have completed
                if(cpDst->r_con.rcv[i].pQueue.queueSize(true) > 0)
		{		    
		    //if(pthread_mutex_trylock(&(cpDst->r_con.rcvLock)) != 0)
		    {		    
			Package *rearPackage =  cpDst->r_con.rcv[i].pQueue.getRear(true); //graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pQueue.getRear(true);
			if(!rearPackage->completed)
			{
			    if(rearPackage->request->ready())
			    {
			      pthread_mutex_lock(&(cpDst->r_con.rcvLock));
                              int ns = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].ns;
                              int nd = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].nd;
                              int lnd = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].lnd;
                              int r_grids = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].r_grids;
                              int tag = rearPackage->tag;
                              //int tag = tagMap[graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pr][g][nd][ns][graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].sz];
			      int dst = upcxx::rank_me();
                              int src= graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pr; 
		                upcxx::rpc(src,
                	            [=](){
                                        pthread_mutex_lock(&(sMsgMap.lock));
                                        //upcxx::future<> *ft= sMsgMap.map[dst][tag].front()->request;
				        //delete ft;//so that sender know
					sMsgMap.map[dst][tag].front()->completed = true;
                                        sMsgMap.map[dst][tag].pop_front();
		                        sMsgMap.size--;
                                        pthread_mutex_unlock(&(sMsgMap.lock));
                                    }
                                );

				delete rearPackage->request;
				rearPackage->completed=true;
				//cpDst->r_con.rcv[i].pQueue.getRear()->completeRequest();
                                //graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pQueue.getRear()->completeRequest();

				//if(graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pQueue.queueSize(true) == 1)
                                if(cpDst->r_con.rcv[i].pQueue.queueSize(true) == 1)
				{
				    cpDst->r_con.firingRuleCnt++;
				}
			        pthread_mutex_unlock(&(cpDst->r_con.rcvLock));
			    }
			}		   		    
			//pthread_mutex_unlock(&(cpDst->r_con.rcvLock));
		    } // if(ga locked)
		} // if(pQueue.queueSize(true) > 0)		    
	    } // for (i<i<cpDst->r_con.nrcv)
	    cpDst = cpDst->next;
	} // while(cpDst != 0)	
    } // for(f<nfabs)

    for(int f=0; f<numfabs; f++)
    {	
	FabCopyAssoc* cpSrc = graphArray[g]->task[f]->cpAsc_srcHead;
	while(cpSrc != 0)
	{
	    for(int i=0; i<cpSrc->r_con.nsnd; i++)
	    {		
		//if(graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].pQueue.queueSize(true) > 0)
                if(cpSrc->r_con.snd[i].pQueue.queueSize(true) >0)
		{
		    //Package *frontPackage = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].pQueue.getFront(true);
                    Package *frontPackage = cpSrc->r_con.snd[i].pQueue.getFront(true);
		    if(frontPackage->served /*&& !frontPackage->completed*/) //!latest receive request has NOT been completed
		    {
			bool flag = false;
			int ret_flag;
			//if(frontPackage->request==NULL)//data have been received by receiver
			if(frontPackage->completed)//data have been received by receiver
			{
/*
			    pthread_mutex_lock(&(graphArray[g]->sCopyMapHead->map[f]->r_con.sndLock));
			    frontPackage = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].pQueue.dequeue(true);
			    frontPackage->completed = false;
			    frontPackage->served = false;
			    frontPackage->request = 0;
			    frontPackage->tag = 0;
			    frontPackage->notified = false;
			    graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].recycleQueue.enqueue(frontPackage,true);
			    pthread_mutex_unlock(&(graphArray[g]->sCopyMapHead->map[f]->r_con.sndLock));
*/

			    pthread_mutex_lock(&(cpSrc->r_con.sndLock));
			    frontPackage = cpSrc->r_con.snd[i].pQueue.dequeue(true);
			    frontPackage->completed = false;
			    frontPackage->served = false;
			    frontPackage->request = 0;
                            frontPackage->notified = false;
			    frontPackage->tag = 0;
			    cpSrc->r_con.snd[i].recycleQueue.enqueue(frontPackage,true);
			    pthread_mutex_unlock(&(cpSrc->r_con.sndLock));			
			}
		    }
		} // if(queueSize > 0)				
	    } // for (i<i<cpSrc->r_con.nsnd)	    
	    cpSrc = cpSrc->next;
	} // while(cpSrc != 0)	
    } // for(f<nfabs)
    upcxx::progress();
} // serviceRemoteGridCopyRequests

#if 0
void Perilla::resetRemoteGridCopyRequests(std::vector<RegionGraph*> graphArray, int g, int nGraphs, int tg)
{
    int np = ParallelDescriptor::NProcs();
    int myProc = ParallelDescriptor::MyProc();
    int numfabs = graphArray[g]->numTasks;
    //MultiFab* mf = graphArray[g]->assocMF;
    int graphID = graphArray[g]->graphID;

    for(int f=0; f<numfabs; f++)
    {	
	if(WorkerThread::isMyRegion(tg,f)) //tg == fg
	{
	    FabCopyAssoc* cpSrc = graphArray[g]->task[f]->cpAsc_srcHead;
	    while(cpSrc != 0)
	    {
		if(cpSrc->r_con.remotePushReady)
		{
		    pthread_mutex_lock(&(graphArray[g]->sCopyMapHead->map[f]->r_con.sndLock));
		    for(int i=0; i<cpSrc->r_con.nsnd; i++)
		    {
			graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].pQueue.enqueue(graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].recycleQueue.dequeue(true),true);	    
		    }
		    pthread_mutex_unlock(&(graphArray[g]->sCopyMapHead->map[f]->r_con.sndLock));
		    cpSrc->r_con.remotePushReady = false;
		}// if remotepushready
		cpSrc = cpSrc->next;
	    }
	}// ismyRegion
    }//for f<numfabs

    for(int f=0; f<numfabs; f++)
    {	
	if(WorkerThread::isMyRegion(tg,f)) //tg == fg
	{
	    FabCopyAssoc* cpDst = graphArray[g]->task[f]->cpAsc_dstHead;
	    while(cpDst != 0)
	    {
		if(pthread_mutex_trylock(&(graphArray[g]->rCopyMapHead->map[f]->r_con.rcvLock)) != 0)
		{
		    if(pthread_mutex_trylock(&(cpDst->r_con.rcvLock)) != 0)
		    {
			//if(f==1 && g==26 && myProc == 54)
			//std::cout<<"Completing Push f " << f << " gID " << g+1 << " myP " << myProc << " PDone "<< cpDst->r_con.remotePullDone <<std::endl;
			if(cpDst->r_con.remotePullDone)
			{
			    for(int i=0; i<cpDst->r_con.nrcv; i++)
			    {

				Package *rcvMetaPackage = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pQueue.dequeue(true);
				rcvMetaPackage->completed = false;
				rcvMetaPackage->served = false;
				rcvMetaPackage->request = 0;	 
				graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].recycleQueue.enqueue(rcvMetaPackage,true);

				Package* rcvPackage = cpDst->r_con.rcv[i].pQueue.dequeue(true);                               // corrected from recycleQ to pQ
				rcvPackage->notified = false;
				rcvPackage->completed = false;
				cpDst->r_con.rcv[i].recycleQueue.enqueue(rcvPackage,true);                         // corrected from pQ to recycleQ

				//cpDst->r_con.firingRuleCnt = cpDst->r_con.firingRuleCnt - 1;

				if(cpDst->r_con.rcv[i].pQueue.queueSize(true) >= 1)
				    if(cpDst->r_con.rcv[i].pQueue.getFront(true)->checkRequest())
					cpDst->r_con.firingRuleCnt++;


			    } // for (i<i<cpDst->r_con.nrcv)

			    cpDst->r_con.remotePullDone = false;

			    //if(f==1 && g==26 && myProc == 54)
			    // std::cout<<"Completed Push f " << f << " gID " << g+1 << " myP " << myProc << " PDone "<< cpDst->r_con.remotePullDone <<std::endl;

			}
			pthread_mutex_unlock(&(cpDst->r_con.rcvLock));
		    } // if(ga locked)
		    pthread_mutex_unlock(&(graphArray[g]->rCopyMapHead->map[f]->r_con.rcvLock));
		} // if(mf locked)
		cpDst = cpDst->next;
	    } // while(cpDst != 0)	
	    /*  
		if(false)
		for(int id=0; id<graphArray[g]->task[f]->depTaskIDs.size(); id++)
		{	    
		int df = graphArray[g]->task[f]->depTaskIDs[id];
		if(WorkerThread::isMyRegion(0,df))
		{
		int lgID = graphArray[g]->srcLinkGraph->graphID-1;

	    //if(f==1 && g==26 && myProc == 54)
	    //std::cout<<"Completing Dep Push f " << df << " gID " << lgID+1 << " myP " << myProc  <<std::endl;

	    FabCopyAssoc* cpdDst = graphArray[lgID]->task[df]->cpAsc_dstHead;
	    while(cpdDst != 0)
	    {
	    if(omp_test_lock(graphArray[lgID]->rCopyMapHead->map[df]->r_con.rcvLock) != 0)
	    {
	    if(omp_test_lock(cpdDst->r_con.rcvLock) != 0)
	    {
	    //if(f==1 && g==26 && myProc == 54)
	    //std::cout<<"Completing Push f " << f << " gID " << g+1 << " myP " << myProc << " PDone "<< cpdDst->r_con.remotePullDone <<std::endl;
	    if(cpdDst->r_con.remotePullDone)
	    {
	    for(int i=0; i<cpdDst->r_con.nrcv; i++)
	    {

	    Package *rcvMetaPackage = graphArray[lgID]->rCopyMapHead->map[df]->r_con.rcv[i].pQueue.dequeue(true);
	    rcvMetaPackage->completed = false;
	    rcvMetaPackage->served = false;
	    rcvMetaPackage->request = MPI_REQUEST_NULL;	  
	    graphArray[lgID]->rCopyMapHead->map[df]->r_con.rcv[i].recycleQueue.enqueue(rcvMetaPackage,true);

	    Package* rcvPackage = cpdDst->r_con.rcv[i].pQueue.dequeue(true);                               // corrected from recycleQ to pQ
	    rcvPackage->notified = false;
	    rcvPackage->completed = false;
	    cpdDst->r_con.rcv[i].recycleQueue.enqueue(rcvPackage,true);                         // corrected from pQ to recycleQ

	    //cpdDst->r_con.firingRuleCnt = cpdDst->r_con.firingRuleCnt - 1;

	    if(cpdDst->r_con.rcv[i].pQueue.queueSize(true) >= 1)
	    if(cpdDst->r_con.rcv[i].pQueue.getFront(true)->checkRequest())
	    cpdDst->r_con.firingRuleCnt++;


	    } // for (i<i<cpdDst->r_con.nrcv)

	    cpdDst->r_con.remotePullDone = false;

	    //if(df==10 && lgID==24 && myProc == 54)
	    // std::cout<<"Completed Push f " << df << " gID " << lgID+1 << " myP " << myProc << " PDone "<< cpdDst->r_con.remotePullDone <<std::endl;
	    }
	    omp_unset_lock(cpdDst->r_con.rcvLock);
	    } // if(ga locked)
	    omp_unset_lock(graphArray[lgID]->rCopyMapHead->map[df]->r_con.rcvLock);
	    } // if(mf locked)
	    cpdDst = cpdDst->next;
	    } // while(cpdDst != 0)	


	    } // if tg==0 region


	    } // for all dependents
	     */



	}
    } // for(f<nfabs)

}
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

void Perilla::multifabCopyPushAsync(RegionGraph* destGraph, RegionGraph* srcGraph, MultiFab* mfDst, MultiFab* mfSrc, int f, bool singleT)
{
    multifabCopyPushAsync(destGraph, srcGraph, mfDst, mfSrc, f, 1, 1, 1, 0, 0, singleT);
}


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

    multifabCopyPush_1Team(destGraph,srcGraph,mfDst,mfSrc,f,dstcomp,srccomp,nc,ng,ngsrc,singleT);
    if(!singleT)
      srcGraph->worker[perilla::wid()]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);

    //double end_time_wtime = omp_get_wtime();
    //if(ntid==0)
      //Perilla::getPPPTimeSplit[2] += end_time_wtime - start_time_wtime;
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
                mfSrc->m_fabs_v[f]->copyToMem(cpSrc->l_con.scpy[i].sbx,srccomp,nc,sndPackage->databuf.local());
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
                  mfSrc->m_fabs_v[f]->copyToMem(cpSrc->l_con.scpy[i].sbx,srccomp,nc,sndPackage->databuf.local());
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
                mfSrc->m_fabs_v[f]->copyToMem(cpSrc->r_con.snd[i].sbx,srccomp,nc,sndPackage->databuf.local());
                sndPackage->notified = false;
                sndPackage->served = false;
                sndPackage->completed = false;
                cpSrc->r_con.snd[i].pQueue.enqueue(sndPackage,true);
              }
            cpSrc->r_con.remotePushReady = true;

            pthread_mutex_unlock(&(cpSrc->r_con.sndLock));

            /*
            pthread_mutex_lock(&(srcGraph->sCopyMapHead->map[f]->r_con.sndLock));
            for(int i=0; i<cpSrc->r_con.nsnd; i++)
              srcGraph->sCopyMapHead->map[f]->r_con.snd[i].pQueue.enqueue(srcGraph->sCopyMapHead->map[f]->r_con.snd[i].recycleQueue.dequeue(true),true);
            pthread_mutex_unlock(&(srcGraph->sCopyMapHead->map[f]->r_con.sndLock));
            */
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
                  mfSrc->m_fabs_v[f]->copyToMem(cpSrc->r_con.snd[i].sbx,srccomp,nc,sndPackage->databuf.local());
                  sndPackage->notified = false;
                  sndPackage->served = false;
                  sndPackage->completed = false;
                  cpSrc->r_con.snd[i].pQueue.enqueue(sndPackage,true);
                }

            //fout.close();         
            srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
            if(ntid==0)
              {
                cpSrc->r_con.remotePushReady = true;
                /*
                pthread_mutex_lock(&(srcGraph->sCopyMapHead->map[f]->r_con.sndLock));
                for(int i=0; i<cpSrc->r_con.nsnd; i++)
                  srcGraph->sCopyMapHead->map[f]->r_con.snd[i].pQueue.enqueue(srcGraph->sCopyMapHead->map[f]->r_con.snd[i].recycleQueue.dequeue(true),true);
                pthread_mutex_unlock(&(srcGraph->sCopyMapHead->map[f]->r_con.sndLock));
                */
                pthread_mutex_unlock(&(cpSrc->r_con.sndLock));
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
    int ntid = perilla::wtid();//-perilla::NUM_COMM_THREADS;

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
                  mf.m_fabs_v[f]->copyFromMem(graph->lMap[f]->l_con.dcpy[i].dbx,0,nComp,dPackage->databuf.local());
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
          rcvMetaPackage->request = 0;
          graph->rMap[f]->r_con.rcv[i].recycleQueue.enqueue(rcvMetaPackage,true);
          Package *rcvPackage = graph->lMap[f]->r_con.rcv[i].pQueue.dequeue(true);

          mf.m_fabs_v[f]->copyFromMem(graph->lMap[f]->r_con.rcv[i].dbx,0,nComp,rcvPackage->databuf.local());
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

