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

volatile int Perilla::numTeamsFinished = 0;
volatile int Perilla::updateMetadata_request = 0;
volatile int Perilla::updateMetadata_noticed = 0;
volatile int Perilla::updateMetadata_done = 0;
int Perilla::max_step=1;
std::map<int,std::map<int,int>> Perilla::pTagCnt;
int Perilla::uTags=0;
bool Perilla::genTags=true;
std::map<int, std::map<int, std::map<int, std::map<int, std::map<int,int> > > > > Perilla::tagMap;
std::map<int, std::map<int, std::map<int, std::map<int, int> > > > Perilla::myTagMap;

void Perilla::clearTagMap(){
    Perilla::tagMap.clear();
}

void Perilla::clearMyTagMap(){
    Perilla::myTagMap.clear();
}

void Perilla::communicateTags()
{
    int myProc = ParallelDescriptor::MyProc();
    int nPs = ParallelDescriptor::NProcs();
    typedef std::map<int, int> tags_t;
    typedef std::map<int, std::map<int,int>> stags_t;
    typedef std::map<int, std::map<int,std::map<int,int>>> dstags_t; 
    typedef std::map<int, std::map<int,std::map<int,std::map<int,int>>>> gdstags_t;
    typedef std::map<int, std::map<int,std::map<int,std::map<int,std::map<int,int>>>>> pgdstags_t;

    int** tags = new int*[nPs];
    int** rtags = new int*[nPs];
    int* rTagCnt = new int[nPs*2];
    int* sTagCnt = new int[nPs*2];

    MPI_Request *srrequest;
    srrequest = new MPI_Request[nPs];
    MPI_Request *ssrequest;
    ssrequest = new MPI_Request[nPs];
    MPI_Request *trrequest;
    trrequest = new MPI_Request[nPs];
    MPI_Request *tsrequest;
    tsrequest = new MPI_Request[nPs];

    std::vector<bool> proc_communicated;

    proc_communicated.resize(nPs);
    for(int p=0; p<nPs; p++)
	proc_communicated[p]=false;

    for(int p=0; p<nPs; p++)
    {
	if(p!=myProc)
	{
	    MPI_Irecv(&rTagCnt[p*2], 2, MPI_INT, p , 0, MPI_COMM_WORLD, &srrequest[p]);
	}
    }

    for(pgdstags_t::iterator it1 = Perilla::tagMap.begin(); it1  != Perilla::tagMap.end(); it1++)
    {
	int tac=0;
	int ng=0;
	for(gdstags_t::iterator it2 = it1->second.begin(); it2  != it1->second.end(); it2++)
	{
	    tac++;
	    tac++;
	    ng++;
	    for(dstags_t::iterator it3 = it2->second.begin(); it3  != it2->second.end(); it3++)
		for(stags_t::iterator it4 = it3->second.begin(); it4  != it3->second.end(); it4++)
		    for(tags_t::iterator it5 = it4->second.begin(); it5  != it4->second.end(); it5++)
		    {
			tac+=4;
		    }
	}
	sTagCnt[it1->first*2] = tac;
	sTagCnt[it1->first*2+1] = ng;
	tags[it1->first] = new int[sTagCnt[it1->first*2]];
	MPI_Isend(&sTagCnt[it1->first*2], 2, MPI_INT, it1->first, 0, MPI_COMM_WORLD, &ssrequest[it1->first]);
	proc_communicated[it1->first]=true;
    }

    for(int p=0; p<nPs; p++)
	if(p!=myProc)
	    if(!proc_communicated[p])
	    {
		sTagCnt[p*2] = 0;
		sTagCnt[p*2+1] = 0;
		MPI_Isend(&sTagCnt[p*2], 2, MPI_INT, p, 0, MPI_COMM_WORLD, &ssrequest[p]);
	    }


    for(pgdstags_t::iterator it1 = Perilla::tagMap.begin(); it1  != Perilla::tagMap.end(); it1++)
    {
	int tac=0;
	for(gdstags_t::iterator it2 = it1->second.begin(); it2  != it1->second.end(); it2++)
	{
	    tags[it1->first][tac++] = it2->first;
	    tags[it1->first][tac++] = pTagCnt[it1->first][it2->first];
	    int gtagc = 0;
	    for(dstags_t::iterator it3 = it2->second.begin(); it3  != it2->second.end(); it3++)
		for(stags_t::iterator it4 = it3->second.begin(); it4  != it3->second.end(); it4++)
		    for(tags_t::iterator it5 = it4->second.begin(); it5  != it4->second.end(); it5++)
		    {
			tags[it1->first][tac++] = it3->first;
			tags[it1->first][tac++] = it4->first;
			tags[it1->first][tac++] = it5->first;
			tags[it1->first][tac++] = it5->second;
			gtagc++;
		    }
	    BL_ASSERT(pTagCnt[it1->first][it2->first] == gtagc);
	}
	MPI_Isend(tags[it1->first], tac, MPI_INT, it1->first, 1, MPI_COMM_WORLD, &tsrequest[it1->first]);
    }


    MPI_Status status;
    for(int p=0; p<nPs; p++)
    {      
	if(p!=myProc)
	{
	    MPI_Wait( &srrequest[p], &status );
	    if(rTagCnt[p*2] > 0)
	    {
		rtags[p] = new int[rTagCnt[p*2]];
		MPI_Irecv(rtags[p], rTagCnt[p*2], MPI_INT, p , 1, MPI_COMM_WORLD, &trrequest[p]);
	    }
	}
    }

    //      //MPI_Irecv(size) Wait


    //MPI_recive tags arra
    for(int p=0; p<nPs; p++)
    {
	if(p!=myProc)
	{
	    if(rTagCnt[p*2] > 0)
	    {
		MPI_Wait( &trrequest[p], &status );
		int tCnt=0;
		for(int g=0; g<rTagCnt[p*2+1];g++)
		{
		    int gi = rtags[p][tCnt++];
		    int sCnt = rtags[p][tCnt++];
		    for(int j=0; j<sCnt; j++)
		    {
			Perilla::myTagMap[gi][rtags[p][tCnt]][rtags[p][tCnt+1]][rtags[p][tCnt+2]]=rtags[p][tCnt+3];
			//std::cout<< "at myP "<<myProc<<" g " << gi << " d " << rtags[p][tCnt]<< " s " << rtags[p][tCnt+1] << " t "<< rtags[p][tCnt+2]<<std::endl;
			tCnt += 4;    
		    }
		}
		//std::cout<< "at P "<< myProc <<" rCnt "<<rTagCnt[p*2]<<" " << tCnt << std::endl;
	    }
	}      	
    }

    for(int p=0; p<nPs; p++)
    {
	if(p!=myProc)
	    if(rTagCnt[p*2] > 0)
	    {
		delete[] rtags[p];
	    }
    }


    for(int p=0; p<nPs; p++)
    {
	if(p!=myProc)
	    if(proc_communicated[p])
	    {
		MPI_Wait( &tsrequest[p], &status );
		delete[] tags[p];
	    }
    }

    delete[] srrequest;
    delete[] ssrequest;
    delete[] trrequest;
    delete[] tsrequest;
    delete[] tags;
    delete[] rtags;
    delete[] rTagCnt;
    delete[] sTagCnt;

    Perilla::genTags=false;
}



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
#pragma omp parallel shared(rg, mf, numfabs, np, TheFB, recv_cctc, send_cctc)
    {
	//int tg = WorkerThread::perilla_wid();
	int fg;
	//if(WorkerThread::perilla_isCommunicationThread())	
#pragma omp single
	{	  
	    //bool cc = !mf->is_nodal(); //  cc = multifab_cell_centered_q(mf)
	    //mf->sMap.reserve(numfabs);
	    //mf->rMap.reserve(numfabs);
	    //std::cout<< "Allocating sMap and rMap" <<std::endl;
	    rg->alloc_lMap(mf);	  
	    rg->alloc_sMap(mf);
	    rg->alloc_rMap(mf);
	}
#pragma omp barrier      
	//if(tid==0)            	          
	{	  
	    //bool cc = !mf->is_nodal(); //  cc = multifab_cell_centered_q(mf)
	    //mf->sMap.reserve(numfabs);
	    //mf->rMap.reserve(numfabs);
#pragma omp for
	    for(int f=0; f<numfabs; f++) //	   !create local communication metadata for each fab
	    {
		//if(WorkerThread::isMyRegion(tg,f) && WorkerThread::perilla_isMasterWorkerThread())
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
#pragma omp barrier
	//now we know how many copying segments each fab owns as source and destination allocate memory for metadata   
#pragma omp for
	for(int f=0; f<numfabs; f++)
	{
	    //fg = f % (omp_get_num_threads()/perilla::NUM_THREADS_PER_TEAM);   /// need to check if computing correct ???????
	    //if((fg == tg) && ((tid%perilla::NUM_THREADS_PER_TEAM)==1))
	    //if(WorkerThread::isMyRegion(tg,f) && WorkerThread::perilla_isMasterWorkerThread())	    
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
#pragma omp barrier
	if(np > 1)
	{
#pragma omp for
	    for(int f=0; f<numfabs; f++)
	    {
		//if(WorkerThread::perilla_isMasterWorkerThread() && WorkerThread::isMyRegion(tg,f))      
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
	    //if(WorkerThread::perilla_isMasterWorkerThread() && tg==0)      
	    {
#pragma omp for
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
#pragma omp parallel shared(mf, numfabs, TheFB, recv_cctc, send_cctc)
    {
	int tg = WorkerThread::perilla_wid();

	//      std::cout<< "Barr 4- "<< tid <<" "<< tg << " " << WorkerThread::isTeamMasterThread(tid) << std::endl;

	//      std::cout<< "Barr 5" <<std::endl;
	int fg, scnt, dcnt;


#pragma omp for
	for(int f=0; f<numfabs; f++)
	{
	    //fg = f % (omp_get_num_threads()/perilla::NUM_THREADS_PER_TEAM);
	    //if((fg == tg) && ((tid%perilla::NUM_THREADS_PER_TEAM)==1))

	    //if((fg == tg) && ((tid%perilla::NUM_THREADS_PER_TEAM)==0))
	    //if(WorkerThread::isMyRegion(tg,f) && WorkerThread::perilla_isMasterWorkerThread())
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

#pragma omp barrier	  

//	if(WorkerThread::perilla_isMasterWorkerThread() && tg==0)    
#pragma omp for
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
	//int tg = WorkerThread::perilla_wid();
	int fg, nsnd, nrcv;

	for(int f=0; f<numfabs; f++)
	{
	    //fg = f % (omp_get_num_threads()/perilla::NUM_THREADS_PER_TEAM);
	    //if((fg == tg) && ((tid%perilla::NUM_THREADS_PER_TEAM)==0))
	    //if(WorkerThread::isMyRegion(tg,f) && WorkerThread::perilla_isMasterWorkerThread())
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

	//if(tid == 0)
	//if(WorkerThread::perilla_isMasterWorkerThread() && tg==0)
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

void Perilla::serviceLocalRequests(RegionGraph* rg, int tg)    
{
    int numfabs = rg->lMap.size();

    for(int f=0; f<numfabs; f++)
    {
	//int fg = f % perilla::NUM_THREAD_TEAMS;

	//	if(tg==0)
	//  std::cout<< "I am tg 0 :) processing fg " << fg <<std::endl;

	if(WorkerThread::isMyRegion(tg,f))
	    //if(tg == fg)
	{
	    //if(tg == 0)
	    //std::cout<<"I am tg " << tg << " starting to process " << f << " in Graph " << graph->graphID <<std::endl;

	    //int lockSucceeded = omp_test_lock(&(rg->lMap[f]->l_con.sLock));
	    //if(lockSucceeded != 0) // 0-Fail, otherwise-Succeed
	    {		
		//if(graph->graphID == 1)
		//if(tg == 0)
		//std::cout<<"I am tg " << tg << " processing " << f << " in Graph " << graph->graphID <<std::endl;
		/*if(graph->graphID == 1 && (f == 2 || f == 1) )
		  {
		  std::cout<< "serviceLR for gID 1  f " << f << " nscpy "<< rg->lMap[f]->l_con.nscpy << std::endl;
		  for(int i=0; i<rg->lMap[f]->l_con.nscpy; i++)
		  std::cout<< " " << rg->lMap[f]->l_con.scpy[i].nd << " " <<  rg->lMap[f]->l_con.scpy[i].dPartner << " " << rg->lMap[f]->l_con.scpy[i].pQueue.queueSize();
		  std::cout<< std::endl;
		  }*/
		for(int i=0; i<rg->lMap[f]->l_con.nscpy; i++){

		    //std::cout<< "serviceLR nscpy " << rg->lMap[f]->l_con.nscpy <<std::endl;

		    //if(graph->graphID == 1 && rg->lMap[f]->l_con.scpy[i].nd == 1)
		    //std::cout<< "Processing gID 1 nd 1 from f " << f << " i " << i << std::endl;

		    if(rg->lMap[f]->l_con.scpy[i].pQueue.queueSize()>0)
		    {
                        omp_set_lock(&(rg->lMap[f]->l_con.sLock));
			Package *sPackage = rg->lMap[f]->l_con.scpy[i].pQueue.dequeue();
			if(perilla::LAZY_PUSH)
			{
			    //  Implemetation deffered. Currently not required
			}
			//if(graph->graphID == 1 && rg->lMap[f]->l_con.scpy[i].nd == 1)
			//std::cout<< "Processing gID 1 nd 1 from f " << f << " i " << i << std::endl;
			omp_set_lock(&(rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dLock));
			int dPartner = rg->lMap[f]->l_con.scpy[i].dPartner;

			//if(rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dcpy[dPartner].recycleQueue.queueSize() == 0 )
			if(dPartner == -1)
			    std::cout<< " Caution rQ size dPrtn "<< rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.ndcpy << " " << dPartner <<" graph ID " <<rg->graphID<<std::endl;
			//std::cout<< " Caution rQ size "<< rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dcpy[dPartner].recycleQueue.queueSize() <<std::endl;

			Package *dPackage = rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dcpy[dPartner].recycleQueue.dequeue(true);

			//for(int j=0; j<sPackage->bufSize; j++)
			//dPackage->databuf[j] = sPackage->databuf[j];        //copy data------------------------------???????????????

			std::memcpy(dPackage->databuf, sPackage->databuf, dPackage->bufSize * sizeof(double));

			rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dcpy[dPartner].pQueue.enqueue(dPackage,true);
			if(rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dcpy[dPartner].pQueue.queueSize(true)==1)
			    rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.firingRuleCnt++;
			//if(graph->graphID == 1 && rg->lMap[f]->l_con.scpy[i].nd == 1)
			//std::cout << "gID 1 frc " << rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.firingRuleCnt << " df " << rg->lMap[f]->l_con.scpy[i].nd <<std::endl;
			omp_unset_lock(&(rg->lMap[rg->lMap[f]->l_con.scpy[i].nd]->l_con.dLock));

			//if(graph->graphID == 1)
			//std::cout<< "Processed gID 1  f " << rg->lMap[f]->l_con.scpy[i].nd << std::endl;

			rg->lMap[f]->l_con.scpy[i].recycleQueue.enqueue(sPackage,true);
		    }}
		omp_unset_lock(&(rg->lMap[f]->l_con.sLock));
#pragma omp flush
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

    // !we first post send and receive  
    for(int f=0; f<numfabs; f++)
    {
	//int lockSucceeded = omp_test_lock(rg->rMap[f]->r_con.rcvLock);
	//if(lockSucceeded != 0)
	{
	    //if(omp_test_lock(rg->lMap[f]->r_con.rcvLock) != 0)
	    {
		for(int i=0; i<rg->lMap[f]->r_con.nrcv; i++)
		{
		    if(rg->rMap[f]->r_con.rcv[i].pQueue.queueSize(true) == 0) //!no message has been received or all received messages have been claimed
			nextsReq = true;
		    else
		    {
			Package *rearPackage = rg->rMap[f]->r_con.rcv[i].pQueue.getRear(true);//!CHECK THIS POINT LATER
			if(rearPackage->completed && rg->rMap[f]->r_con.rcv[i].pQueue.queueSize(true) == 1) //!latest receive request has been completed
			    nextsReq = true;
			else //!expected message is still on the way
			    nextsReq = false;
		    }
		    if(nextsReq) //!take a message from recycle pool and post a receive
		    {
                        omp_set_lock((rg->rMap[f]->r_con.rcvLock));
                        omp_set_lock((rg->lMap[f]->r_con.rcvLock));
			int ns = rg->rMap[f]->r_con.rcv[i].ns;
			int nd = rg->rMap[f]->r_con.rcv[i].nd;
			int lnd = rg->rMap[f]->r_con.rcv[i].lnd;
			int r_grids = rg->rMap[f]->r_con.rcv[i].r_grids;
			//!create a package to keep track of receive requests
			Package *rMetaPackage = rg->rMap[f]->r_con.rcv[i].recycleQueue.dequeue(true);
			//!extract a package from the recycle pool at the destination NUMA node to buffer incoming data
			Package *rPackage = rg->lMap[f]->r_con.rcv[i].recycleQueue.dequeue(true);
			//tag = tagGen(mf%rMap(f)%r_con%rcv(i)%ns, mf%rMap(f)%r_con%rcv(i)%nd, gid, parallel_nprocs()*nfabs(mf), ngr)---------??????
			//int tag = tagGen(rg->rMap[f]->r_con.rcv[i].ns, rg->rMap[f]->r_con.rcv[i].nd, graphID-1, np*numfabs, nGraphs);			
			int tag = tagMap[rg->rMap[f]->r_con.rcv[i].pr][graphID][nd][ns][rg->rMap[f]->r_con.rcv[i].sz];

			rMetaPackage->request = MPI_REQUEST_NULL;
			rg->lMap[f]->r_con.rcv[i].pQueue.enqueue(rPackage,true);   //!this is not done yet
			rg->rMap[f]->r_con.rcv[i].pQueue.enqueue(rMetaPackage,true);   //!this is not done yet
			//rMetaPackage->request = parallel_irecv_dv(rpackage%ptr%dataBuf,mf%rMap(f)%r_con%rcv(i)%sz, mf%rMap(f)%r_con%rcv(i)%pr, tag) --------- ????
			rMetaPackage->request = ParallelDescriptor::Arecv(rPackage->databuf,
				rg->rMap[f]->r_con.rcv[i].sz,
				rg->rMap[f]->r_con.rcv[i].pr, tag).req(); // tag == SeqNum in c++ ver
                        omp_unset_lock((rg->lMap[f]->r_con.rcvLock));
                        omp_unset_lock((rg->rMap[f]->r_con.rcvLock));
		    }
		}
		//omp_unset_lock(rg->lMap[f]->r_con.rcvLock);
	    }// if(omp_test_lock)
	    //omp_unset_lock(rg->rMap[f]->r_con.rcvLock);
	}// if(lockSucceeded)
    }// for(f<numfabs)


    for(int f=0; f<numfabs; f++)
    {
	for(int i=0; i<rg->sMap[f]->r_con.nsnd; i++)
	{
	    if(rg->sMap[f]->r_con.snd[i].pQueue.queueSize(true) == 0) //then !no message has been issued or all send requests have been fulfilled
		nextrReq = false;
	    else
		nextrReq = true;

	    if(nextrReq)
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
		    //tag = tagGen(mf%sMap(f)%r_con%snd(i)%ns, mf%sMap(f)%r_con%snd(i)%nd, gid, parallel_nprocs()*nfabs(mf), ngr) -???????
		    //int tag = tagGen(rg->sMap[f]->r_con.snd[i].ns, rg->sMap[f]->r_con.snd[i].nd, graphID-1, np*numfabs, nGraphs);
		    int tag = Perilla::myTagMap[r_gid][nd][ns][rg->sMap[f]->r_con.snd[i].sz];
		    //int tag = myTagMap[graphID-1][rg->sMap[f]->r_con.snd[i].nd][rg->sMap[f]->r_con.snd[i].ns];
		    //sMetaPackage%ptr%request = parallel_isend_dv(spackage%ptr%dataBuf,mf%sMap(f)%r_con%snd(i)%sz, mf%sMap(f)%r_con%snd(i)%pr, tag) --?????
		    sMetaPackage->request = ParallelDescriptor::Asend(sPackage->databuf,
			    rg->sMap[f]->r_con.snd[i].sz,
			    rg->sMap[f]->r_con.snd[i].pr, tag).req();  // tag == SeqNum in c++ ver
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
	    if(rg->rMap[f]->r_con.rcv[i].pQueue.queueSize(true) > 0) //!all messages before rear have completed
	    {
		//if(omp_test_lock(rg->lMap[f]->r_con.rcvLock) != 0) // 0-Fail, otherwise-Succeed
		{
		    Package *rearPackage =  rg->rMap[f]->r_con.rcv[i].pQueue.getRear(true);
		    if(!rearPackage->completed)
		    {
			bool flag = false;
			int ret_flag;
			MPI_Status status;

			std::cout<< "myP "<< myProc << " f "<< f << " i "<< i<< " Req "<<rearPackage->request << std::endl;

			ParallelDescriptor::Test(rearPackage->request, ret_flag, status);
			flag = (ret_flag == 0) ? false : true;//parallel_test_one(rearPackage%ptr%request) -------???????
			if(flag)
			{
                            omp_set_lock((rg->lMap[f]->r_con.rcvLock));
			    rearPackage->completeRequest();
			    rg->lMap[f]->r_con.rcv[i].pQueue.getRear()->completeRequest();
			    if(rg->rMap[f]->r_con.rcv[i].pQueue.queueSize(true) == 1)
				rg->lMap[f]->r_con.firingRuleCnt++;
			    omp_unset_lock((rg->lMap[f]->r_con.rcvLock));
#pragma omp flush			    
			}
		    }
		    //omp_unset_lock(rg->lMap[f]->r_con.rcvLock);
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
		if(frontPackage->served && !frontPackage->completed) //!latest receive request has NOT been completed
		{
		    bool flag = false;
		    int ret_flag;
		    MPI_Status status;
		    ParallelDescriptor::Test(frontPackage->request, ret_flag, status);
		    flag = (ret_flag == 0) ? false : true;//parallel_test_one(frontPackage%ptr%request) -------???????		    
		    if(flag)
		    {
			omp_set_lock(rg->sMap[f]->r_con.sndLock);
			frontPackage = rg->sMap[f]->r_con.snd[i].pQueue.dequeue(true);
			frontPackage->completed = false;
			frontPackage->served = false;
			frontPackage->request = MPI_REQUEST_NULL;
			frontPackage->notified = false;
			rg->sMap[f]->r_con.snd[i].recycleQueue.enqueue(frontPackage,true);
			omp_unset_lock(rg->sMap[f]->r_con.sndLock);
#pragma omp flush
			omp_set_lock(rg->lMap[f]->r_con.sndLock);
			frontPackage = rg->lMap[f]->r_con.snd[i].pQueue.dequeue(true);
			frontPackage->completed = false;
			frontPackage->served = false;
			frontPackage->request = MPI_REQUEST_NULL;
			rg->lMap[f]->r_con.snd[i].recycleQueue.enqueue(frontPackage,true);
			omp_unset_lock(rg->lMap[f]->r_con.sndLock);			
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
#pragma omp flush (graph)
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
#pragma omp flush (graph)
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
#pragma omp flush (graphArray)
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

void Perilla::flattenGraphHierarchy(std::vector<std::vector<RegionGraph*> > graphArrayHierarchy, std::vector<RegionGraph*> &graphArray){
    int gCnt=0;
    for(int l=0; l<graphArrayHierarchy.size(); l++) gCnt+= graphArrayHierarchy[l].size();
    for(int l=0; l<graphArrayHierarchy.size(); l++)
        for(int g=0; g<graphArrayHierarchy[l].size(); g++)
            graphArray.push_back(graphArrayHierarchy[l][g]);
}

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
	//lstime = omp_get_wtime();
	for(int g=0; g<graphArray.size(); g++)
	{
	    nGraphs = graphArray.size();
	    //if(graphArray[g]->graphID==13)
	    //std::cout<<"Processing Local GridCopy Req Graph "<< g+1 << " tg " << tg <<std::endl;
#pragma omp flush (graphArray)
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
		if(np > 1)
		    //if(tg==0)
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
	/*
	//!check if we have finished all the graph execution
	bool noMoreWork = true;
	//std::cout<<"Graph Not Completed ";
	for(int g=0; g<nGraphs; g++)
	if(graphArray[g]->totalFinishes < perilla::NUM_THREAD_TEAMS)
	{
	noMoreWork = false;
	//if(tg==0)
	//std::cout<< g << " tfs " << graphArray[g]->totalFinishes << std::endl;
	}
	//else
	// std::cout<<"Graph Completed "<< g <<std::endl;
	//std::cout<<std::endl;
	if(noMoreWork)
	break;
	 */

	//std::cin.ignore( std::numeric_limits <std::streamsize> ::max(), '\n' );

	//for(int g=0; g<graphArray.size(); g++)
	//std::cout << g+1 << ":" << graphArray[g]->totalFinishes << " | ";


	//f( Perilla::numTeamsFinished == perilla::NUM_THREAD_TEAMS)
	//{
	//    if(doublechecked) // double check if there are still something to send
	//	break;
	//    else
	//	doublechecked = true;
	//}

	//std::cout<<"Teams Completed "<< Perilla::numTeamsFinished << " tid "<< tid << " myProc " << myProc <<std::endl;

	//letime = omp_get_wtime();
	numloops++;
	//ltime = letime - lstime;

	//avgltime += ltime;
	//if(ltime < minltime)
	 //   minltime = ltime;
	//if(ltime > maxltime)
	//    maxltime = ltime;

    } // while(true)

    //if(myProc==0)
    //std::cout<< std::endl << "COMM HANDLER TIMES tg" << tg << " avg " << avgltime/numloops << " min " << minltime << " max " << maxltime <<std::endl;

    //std::cout<< std::endl << "COMM HANDLER " << tg << " FINISHED EXECUTION" << " myProc " << myProc << " nTF " << Perilla::numTeamsFinished << " nTT " << perilla::NUM_THREAD_TEAMS<<std::endl;

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


void Perilla::serviceMultipleGraphComm(RegionGraph graphArray[], int nGraphs, int tid)
{
    serviceMultipleGraphComm(graphArray,nGraphs,false,tid);
} // serviceMultipleGraphComm

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
	    omp_set_lock(&(graph->lMap[f]->l_con.sLock));
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
		    mf->m_fabs_v[f]->copyToMem(graph->lMap[f]->l_con.scpy[i].sbx,0,nComp,sPackage->databuf);

		    for(int d=0; d<sPackage->bufSize; d++)
			if(sPackage->databuf[d] == 0)
			{
			    //std::cout<< "in fbPush Sending 0 from f "<< f <<std::endl;
			    //BL_ASSERT(sPackage->databuf[d] != 0);
			}
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
	    omp_unset_lock(&(graph->lMap[f]->l_con.sLock));
	}
	graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads
    } // if(LAZY_PUSH) - else

    int np = ParallelDescriptor::NProcs();
    if (np==1) return;

    if(ntid==0)
	omp_set_lock(graph->lMap[f]->r_con.sndLock);
    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads

    for(int i=0; i<graph->lMap[f]->r_con.nsnd; i++)
	if((i%(perilla::NUM_THREADS_PER_TEAM-1))==ntid)
	{
	    //std::cout << "RQS " << graph->lMap[f]->r_con.snd[i].recycleQueue.queueSize() << std::endl;

	    Package *sndPackage = graph->lMap[f]->r_con.snd[i].recycleQueue.dequeue(true);	  
	    mf->m_fabs_v[f]->copyToMem(graph->lMap[f]->r_con.snd[i].sbx,0,nComp,sndPackage->databuf);
	    sndPackage->notified = false;
	    graph->lMap[f]->r_con.snd[i].pQueue.enqueue( sndPackage,true );
	    //!the local message handler will detect the change and notify the remote message handler =>read access
	    //!the remote message handler first modifies the front item of this queue, then it push this item back to the message pool
	}
    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads
    if(ntid==0)
    {
	omp_unset_lock(graph->lMap[f]->r_con.sndLock);
	omp_set_lock(graph->sMap[f]->r_con.sndLock);
	for(int i=0; i<graph->lMap[f]->r_con.nsnd; i++)
	    graph->sMap[f]->r_con.snd[i].pQueue.enqueue( graph->sMap[f]->r_con.snd[i].recycleQueue.dequeue(true),true );
	omp_unset_lock(graph->sMap[f]->r_con.sndLock);
    }    																					      

} // fillBoundaryPush

void Perilla::fillBoundaryPull(RegionGraph* graph, MultiFab* mf, int f)
{

    int nComp = mf->nComp();
    int tg= WorkerThread::perilla_wid();
    int ntid = WorkerThread::perilla_wtid();

    if(ntid==0)
	omp_set_lock(&(graph->lMap[f]->l_con.dLock));
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

		    for(int d=0; d<dPackage->bufSize; d++)
			if(dPackage->databuf[d] == 0)
			{
			    //std::cout<< "in fbPull Reciving 0 for f "<< f <<std::endl;
			    //BL_ASSERT(dPackage->databuf[d] != 0);
			}
		    /*
		       if(f==0)
		    //if(graph->lMap[f]->l_con.dcpy[i].dbx.smallEnd() == graph->lMap[f]->l_con.dcpy[i].dbx.bigEnd())
		    //if(graph->lMap[f]->l_con.dcpy[i].dbx.smallEnd(0)==-1 && graph->lMap[f]->l_con.dcpy[i].dbx.smallEnd(1)==-1 && graph->lMap[f]->l_con.dcpy[i].dbx.smallEnd(2)==4)
		    std::cout<< "Corner Pull for f "<< f << " data0 " <<dPackage->databuf[0]<< " size " <<dPackage->bufSize <<" se " <<graph->lMap[f]->l_con.dcpy[i].dbx.smallEnd()<<std::endl;
		     */
		    mf->m_fabs_v[f]->copyFromMem(graph->lMap[f]->l_con.dcpy[i].dbx,0,nComp,dPackage->databuf);		  
		}
	} // if(UNPACKING_FINEGRAIN) - else
    } // if(LAZY_PUSH) - else

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
	omp_unset_lock(&(graph->lMap[f]->l_con.dLock));
    }
    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads

    int np = ParallelDescriptor::NProcs();
    if (np==1) return;

    if(ntid==0)
    {
	omp_set_lock(graph->rMap[f]->r_con.rcvLock);
	omp_set_lock(graph->lMap[f]->r_con.rcvLock);
    }
    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads

    for(int i=0; i<graph->lMap[f]->r_con.nrcv; i++)
	if( (i%(perilla::NUM_THREADS_PER_TEAM-1)) == ntid)
	{
	    Package *rcvMetaPackage = graph->rMap[f]->r_con.rcv[i].pQueue.dequeue(true);
	    rcvMetaPackage->completed = false;
	    rcvMetaPackage->served = false;
	    rcvMetaPackage->request = MPI_REQUEST_NULL;
	    graph->rMap[f]->r_con.rcv[i].recycleQueue.enqueue(rcvMetaPackage,true);
	    Package *rcvPackage = graph->lMap[f]->r_con.rcv[i].pQueue.dequeue(true);
	    mf->m_fabs_v[f]->copyFromMem(graph->lMap[f]->r_con.rcv[i].dbx,0,nComp,rcvPackage->databuf);
	    rcvPackage->completed = false;
	    rcvPackage->notified = false;
	    graph->lMap[f]->r_con.rcv[i].recycleQueue.enqueue(rcvPackage,true);
	}
    graph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads

    if(ntid==0)
    {
	graph->lMap[f]->r_con.firingRuleCnt = graph->lMap[f]->r_con.firingRuleCnt - graph->lMap[f]->r_con.nrcv;
	for(int i=0; i<graph->lMap[f]->r_con.nrcv; i++)
	    if(graph->lMap[f]->r_con.rcv[i].pQueue.queueSize(true) >= 1)
		if(graph->lMap[f]->r_con.rcv[i].pQueue.getFront(true)->checkRequest())
		    graph->lMap[f]->r_con.firingRuleCnt++;
	omp_unset_lock(graph->lMap[f]->r_con.rcvLock);
	omp_unset_lock(graph->rMap[f]->r_con.rcvLock);
    }

} // fillBoundaryPull

void Perilla::fillBoundaryPull(RegionGraph* graph, MultiFab* mf, int f, bool singleT)
{
exit(0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

	    //	std::cout<< "Loop 1" <<std::endl;

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

#pragma omp parallel shared(gSrc, gDst, mfSrc, mfDst, nfabsSrc, nfabsDst)
	    {
		int tg = WorkerThread::perilla_wid();
		int fg;

		for(int f=0; f<nfabsSrc; f++)
		{
		    if(WorkerThread::isMyRegion(tg,f) && WorkerThread::perilla_isMasterWorkerThread())
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
#pragma omp barrier
		    //	      std::cout<< "Barrier 1" <<std::endl;	      
		    if(np > 1)
		    {
			if(WorkerThread::perilla_isMasterWorkerThread() && tg==0)
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
#pragma omp barrier	  	  
		for(int f=0; f<nfabsDst; f++)
		{
		    if(WorkerThread::isMyRegion(tg,f) && perilla::isMasterWorkerThread())		
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

#pragma omp barrier		  
		    if(np > 1)
		    {
			//if(tid==0)
			if(WorkerThread::perilla_isMasterWorkerThread() && tg==0)			
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
#pragma omp barrier
	    for(int f=0; f<nfabsSrc; f++)
	    {
		if(WorkerThread::isMyRegion(tg,f) && WorkerThread::perilla_isMasterWorkerThread())	      
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
#pragma omp barrier
	    for(int f=0; f<nfabsDst; f++)
	    {
		if(WorkerThread::isMyRegion(tg,f) && WorkerThread::perilla_isMasterWorkerThread())
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

void Perilla::multifabExtractCopyAssoc(RegionGraph* gDst, RegionGraph* gSrc, const MultiFab& mfDst, const MultiFab& mfSrc, const Periodicity& period) 
{
    multifabExtractCopyAssoc(gDst, gSrc, mfDst, mfSrc, 1, 0, 0, period);
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
	    omp_set_lock(&(cpSrc->l_con.sLock));		    
	    for(int i=0; i<cpSrc->l_con.nscpy; i++)
	    {
		Package* sndPackage = cpSrc->l_con.scpy[i].recycleQueue.getFront(true);
		mfSrc->m_fabs_v[f]->copyToMem(cpSrc->l_con.scpy[i].sbx,srccomp,nc,sndPackage->databuf);
	    }	    
	    for(int i=0;i<cpSrc->l_con.nscpy; i++)
		cpSrc->l_con.scpy[i].pQueue.enqueue(cpSrc->l_con.scpy[i].recycleQueue.dequeue(true),true);
	    omp_unset_lock(&(cpSrc->l_con.sLock));   
	}
	else
	{
	    if(ntid == 0)
		omp_set_lock(&(cpSrc->l_con.sLock));	
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
		omp_unset_lock(&(cpSrc->l_con.sLock));   
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
	    omp_set_lock(cpSrc->r_con.sndLock);
	    for(int i=0; i<cpSrc->r_con.nsnd; i++)
	    {

		Package* sndPackage = cpSrc->r_con.snd[i].recycleQueue.dequeue(true);
		mfSrc->m_fabs_v[f]->copyToMem(cpSrc->r_con.snd[i].sbx,srccomp,nc,sndPackage->databuf);
		sndPackage->notified = false;
		sndPackage->notified = false;
		cpSrc->r_con.snd[i].pQueue.enqueue(sndPackage,true);
	    }

	    omp_unset_lock(cpSrc->r_con.sndLock); 

	    cpSrc->r_con.remotePushReady = true;
	    ///*
	    omp_set_lock(srcGraph->sCopyMapHead->map[f]->r_con.sndLock);
	    for(int i=0; i<cpSrc->r_con.nsnd; i++)
		srcGraph->sCopyMapHead->map[f]->r_con.snd[i].pQueue.enqueue(srcGraph->sCopyMapHead->map[f]->r_con.snd[i].recycleQueue.dequeue(true),true);	    
	    omp_unset_lock(srcGraph->sCopyMapHead->map[f]->r_con.sndLock);
	}
	else
	{
	    if(ntid == 0)
		omp_set_lock(cpSrc->r_con.sndLock);
	    srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);

	    for(int i=0; i<cpSrc->r_con.nsnd; i++)
		if((i%(perilla::NUM_THREADS_PER_TEAM-1)) == ntid)
		{

		    // if(myProc==4 && srcGraph->graphID==2 && (f ==0 || f ==2))
		    //std::cout << " Pushing 2 316 164"<<std::endl;

		    Package* sndPackage = cpSrc->r_con.snd[i].recycleQueue.dequeue(true);
		    mfSrc->m_fabs_v[f]->copyToMem(cpSrc->r_con.snd[i].sbx,srccomp,nc,sndPackage->databuf);
		    sndPackage->notified = false;
		    sndPackage->notified = false;
		    cpSrc->r_con.snd[i].pQueue.enqueue(sndPackage,true);

		}

	    srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
	    if(ntid==0)
	    {
		omp_unset_lock(cpSrc->r_con.sndLock); 

		cpSrc->r_con.remotePushReady = true;
		///*
		omp_set_lock(srcGraph->sCopyMapHead->map[f]->r_con.sndLock);
		for(int i=0; i<cpSrc->r_con.nsnd; i++)
		    srcGraph->sCopyMapHead->map[f]->r_con.snd[i].pQueue.enqueue(srcGraph->sCopyMapHead->map[f]->r_con.snd[i].recycleQueue.dequeue(true),true);	    
		omp_unset_lock(srcGraph->sCopyMapHead->map[f]->r_con.sndLock);
		//*/
	    }
	    srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
	}
    } // if(!(*mfDst == *mfSrc))      													      
} // multifabCopyPushAsync

void Perilla::multifabCopyPushAsync(RegionGraph* destGraph, RegionGraph* srcGraph, MultiFab* mfDst, MultiFab* mfSrc, int f, bool singleT) 
{
    multifabCopyPushAsync(destGraph, srcGraph, mfDst, mfSrc, f, 1, 1, 1, 0, 0, singleT);
} 

void Perilla::multifabCopyPush(RegionGraph* destGraph, RegionGraph* srcGraph, amrex::MultiFab* mfDst, amrex::MultiFab* mfSrc, int f, int dstcomp, int srccomp, int nc, int ng, int ngsrc, bool singleT)
{
    if(nc<1) cout <<"MULTIFAB_COPY_C: nc must be >= 1"<< endl;
    if(mfDst->nComp() < (dstcomp-1)) cout <<"MULTIFAB_COPY_C: nc too large for dst multifab"<< endl;
    if(mfSrc->nComp() < (srccomp-1)) cout <<"MULTIFAB_COPY_C: nc too large for src multifab"<< endl;

    multifabCopyPush_1Team(destGraph,srcGraph,mfDst,mfSrc,f,dstcomp,srccomp,nc,ng,ngsrc,singleT);
    if(!singleT)
      srcGraph->worker[perilla::wid()]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
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
            omp_set_lock(&(cpSrc->l_con.sLock));
            for(int i=0; i<cpSrc->l_con.nscpy; i++)
              {
                Package* sndPackage = cpSrc->l_con.scpy[i].recycleQueue.getFront(true);
                mfSrc->m_fabs_v[f]->copyToMem(cpSrc->l_con.scpy[i].sbx,srccomp,nc,sndPackage->databuf);
              }
            for(int i=0;i<cpSrc->l_con.nscpy; i++)
              cpSrc->l_con.scpy[i].pQueue.enqueue(cpSrc->l_con.scpy[i].recycleQueue.dequeue(true));
            omp_unset_lock(&(cpSrc->l_con.sLock));
        }
        else
          {
            if(ntid == 0)
              omp_set_lock(&(cpSrc->l_con.sLock));
            srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);

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
                  cpSrc->l_con.scpy[i].pQueue.enqueue(cpSrc->l_con.scpy[i].recycleQueue.dequeue(true));
                omp_unset_lock(&(cpSrc->l_con.sLock));
              }
            srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
          }

        int np = amrex::ParallelDescriptor::NProcs();
        if(np == 1)
          return;
        if(singleT)
        {
            omp_set_lock((cpSrc->r_con.sndLock));
            for(int i=0; i<cpSrc->r_con.nsnd; i++)
            {
                Package* sndPackage = cpSrc->r_con.snd[i].recycleQueue.dequeue(true);
                mfSrc->m_fabs_v[f]->copyToMem(cpSrc->r_con.snd[i].sbx,srccomp,nc,sndPackage->databuf);
                sndPackage->notified = false;
                sndPackage->served = false;
                sndPackage->completed = false;
                cpSrc->r_con.snd[i].pQueue.enqueue(sndPackage, true);
            }
            cpSrc->r_con.remotePushReady = true;
            omp_unset_lock((cpSrc->r_con.sndLock));
        }
        else
        {
            if(ntid == 0)
            {
                omp_set_lock((cpSrc->r_con.sndLock));
            }
            srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);

            for(int i=0; i<cpSrc->r_con.nsnd; i++)
              if((i%(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS)) == ntid)
                {
                  Package* sndPackage = cpSrc->r_con.snd[i].recycleQueue.dequeue(true);
                  mfSrc->m_fabs_v[f]->copyToMem(cpSrc->r_con.snd[i].sbx,srccomp,nc,sndPackage->databuf);
                  sndPackage->notified = false;
                  sndPackage->served = false;
                  sndPackage->completed = false;
                  cpSrc->r_con.snd[i].pQueue.enqueue(sndPackage, true);
                }

            //fout.close();         
            srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
            if(ntid==0)
              {
                cpSrc->r_con.remotePushReady = true;
/*
                for(int i=0; i<cpSrc->r_con.nsnd; i++){
                  Package* sndPackage = srcGraph->sCopyMapHead->map[f]->r_con.snd[i].recycleQueue.dequeue(true);
                  sndPackage->served = false;
                  sndPackage->completed = false;
                  srcGraph->sCopyMapHead->map[f]->r_con.snd[i].pQueue.enqueue(sndPackage, true);
                }
*/
                omp_unset_lock((cpSrc->r_con.sndLock));
            }
            srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
         }
      } // if(!(*mfDst == *mfSrc))                                                                                                                    
  } // multifabCopyPush


void Perilla::multifabCopyPush(RegionGraph* destGraph, RegionGraph* srcGraph, amrex::MultiFab* mfDst, amrex::MultiFab* mfSrc, int f, bool singleT)
  {
    multifabCopyPush(destGraph, srcGraph, mfDst, mfSrc, f, 1, 1, 1, 0, 0, singleT);
  }

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
	    omp_set_lock(&(cpDst->l_con.dLock));
	    for(int i=0; i<cpDst->l_con.ndcpy; i++)
	    {
		Package* rcvPackage = cpDst->l_con.dcpy[i].pQueue.getFront(true); // corrected from recycleQ to pQ
		mfDst->m_fabs_v[f]->copyFromMem(cpDst->l_con.dcpy[i].dbx,dstcomp,nc,rcvPackage->databuf);
	    }	
	    for(int i=0; i<cpDst->l_con.ndcpy; i++)
		cpDst->l_con.dcpy[i].recycleQueue.enqueue(cpDst->l_con.dcpy[i].pQueue.dequeue(true),true); // corrected from pQ to recycleQ and from recycleQ to pQ
	    cpDst->l_con.firingRuleCnt = cpDst->l_con.firingRuleCnt - cpDst->l_con.ndcpy;
	    omp_unset_lock(&(cpDst->l_con.dLock));
	}
	else
	{
	    if(ntid==0)
		omp_set_lock(&(cpDst->l_con.dLock));
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
		omp_unset_lock(&(cpDst->l_con.dLock));
	    }
	    destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
	}

	int np = ParallelDescriptor::NProcs();
	if(np == 1)
	    return;

	if(singleT)
	{
	    omp_set_lock(destGraph->rCopyMapHead->map[f]->r_con.rcvLock);
	    omp_set_lock(cpDst->r_con.rcvLock);
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
	    omp_unset_lock(cpDst->r_con.rcvLock);
	    omp_unset_lock(destGraph->rCopyMapHead->map[f]->r_con.rcvLock);

	}
	else
	{	
	    if(ntid==0)
	    {
		omp_set_lock(destGraph->rCopyMapHead->map[f]->r_con.rcvLock);
		omp_set_lock(cpDst->r_con.rcvLock);
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
		omp_unset_lock(cpDst->r_con.rcvLock);
		omp_unset_lock(destGraph->rCopyMapHead->map[f]->r_con.rcvLock);
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
		//std::cout<<" "<<cpSrc << " ";
		int lockSucceeded = omp_test_lock(&(cpSrc->l_con.sLock));
		if(lockSucceeded != 0)
		{
		    for(int i=0; i<cpSrc->l_con.nscpy; i++)
		    {
			if(cpSrc->l_con.scpy[i].pQueue.queueSize()>0)
			{
assert(doublechecked==false);
			    FabCopyAssoc* cpDst = cpSrc->graphPartner->task[cpSrc->l_con.scpy[i].nd]->cpAsc_dstHead;
			    while(cpDst != 0)
			    {
				if(cpDst->graphPartner == graphArray[g])
				    break;
				cpDst = cpDst->next;
			    }			    
			    Package* sPackage = cpSrc->l_con.scpy[i].pQueue.dequeue(true);
			    omp_set_lock(&(cpDst->l_con.dLock));
			    int dPartner = cpSrc->l_con.scpy[i].dPartner;
			    Package* dPackage = cpDst->l_con.dcpy[dPartner].recycleQueue.dequeue(true);
			    /*
			       for(int j=0; j<dPackage->bufSize; j++)
			       {
			       dPackage->databuf[j] = sPackage->databuf[j];
			       }
			     */
			    std::memcpy(dPackage->databuf, sPackage->databuf, dPackage->bufSize * sizeof(double));
			    //std::swap(dPackage->databuf, sPackage->databuf);

			    cpDst->l_con.dcpy[dPartner].pQueue.enqueue(dPackage,true);
			    if(cpDst->l_con.dcpy[dPartner].pQueue.queueSize(true) == 1)
				cpDst->l_con.firingRuleCnt++;
			    omp_unset_lock(&(cpDst->l_con.dLock));
			    cpSrc->l_con.scpy[i].recycleQueue.enqueue(sPackage,true);
			}
		    } // for
		    omp_unset_lock(&(cpSrc->l_con.sLock));
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
    //MultiFab* mf = graphArray[g]->assocMF;
    int graphID = graphArray[g]->graphID;

    for(int f=0; f<numfabs; f++)
    {	
	FabCopyAssoc* cpDst = graphArray[g]->task[f]->cpAsc_dstHead;
	while(cpDst != 0)
	{
	    if(omp_test_lock(graphArray[g]->rCopyMapHead->map[f]->r_con.rcvLock) != 0)
	    {
		if(omp_test_lock(cpDst->r_con.rcvLock) != 0)
		{
		    for(int i=0; i<cpDst->r_con.nrcv; i++)
		    {
			if(graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pQueue.queueSize(true) == 0) //!no message has been received or all received messages have been claimed
			{
			    nextsReq = true;
			}
			else
			{			    
			    Package *rearPackage = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pQueue.getRear(true);//!CHECK THIS POINT LATER
			    // Also check the recycle queue because when rear is completed it may cause unlimited recv posts
			    if(rearPackage->completed && graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].recycleQueue.queueSize() > 1) //!latest receive request has been completed
			    {
				nextsReq = true;
			    }
			    else //!expected message is still on the way
				nextsReq = false;
			}
			if(nextsReq) //!take a message from recycle pool and post a receive
			{
			    //!create a package to keep track of receive requests

			    Package *rMetaPackage = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].recycleQueue.dequeue(true);
			    //!extract a package from the recycle pool at the destination NUMA node to buffer incoming data
			    int ns = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].ns;
			    int nd = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].nd;
			    int lnd = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].lnd;
			    int r_grids = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].r_grids;
			    Package *rPackage = cpDst->r_con.rcv[i].recycleQueue.dequeue(true);
			    //int tag = tagGen(ns, nd, graphID-1, np*r_grids, nGraphs);
			    //int tag = Perilla::myTagMap[graphID-1][nd][ns];
			    //int tag = graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].lnd;
			    int tag = tagMap[graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pr][g][nd][ns][graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].sz];

			    //  if(graphArray[g]->graphID == 25 && lnd==10 && myProc==54)
			    //std::cout << "R Posted g " << g << " myP " << myProc << " lnd " << lnd <<" nd "<< nd << " ns "<<ns << " tag "<<tag << " pr " <<graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pr << std::endl;

			    rMetaPackage->request = MPI_REQUEST_NULL;
			    cpDst->r_con.rcv[i].pQueue.enqueue(rPackage,true);   //!this is not done yet
			    graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pQueue.enqueue(rMetaPackage,true);   //!this is not done yet	 
			    rMetaPackage->request = ParallelDescriptor::Arecv(rPackage->databuf,
				    graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].sz,
				    graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pr, tag).req(); // tag == SeqNum in c++ ver

			}						
		    } // for (i<i<cpDst->r_con.nrcv)
		    omp_unset_lock(cpDst->r_con.rcvLock);
		} // if(ga locked)
		omp_unset_lock(graphArray[g]->rCopyMapHead->map[f]->r_con.rcvLock);
	    } // if(mf locked)
	    cpDst = cpDst->next;
	} // while(cpDst != 0)	
    } // for(f<nfabs)

    for(int f=0; f<numfabs; f++)
    {	

	//	if(g == 17 && f == 316 )
	///std::cout << "Trying S Post " << std::endl;

	FabCopyAssoc* cpSrc = graphArray[g]->task[f]->cpAsc_srcHead;
	while(cpSrc != 0)
	{
	    for(int i=0; i<cpSrc->r_con.nsnd; i++)
	    {
		//if(g == 17 && f == 316 && i == 164)
		//std::cout << "Comm Thread nsnd "<< cpSrc->r_con.nsnd << " " << graphArray[g]<< std::endl;
		if(graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].pQueue.queueSize(true) == 0) //!no message has been received or all received messages have been claimed	       	
		    nextrReq = false;
		else
		    nextrReq = true;

		if(nextrReq) //!take a message from recycle pool and post a receive
		{

		    Package *sMetaPackage = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].pQueue.getFront(true);
		    if(!sMetaPackage->served)
		    {		    
			Package *sPackage = cpSrc->r_con.snd[i].pQueue.getFront(true);
			sMetaPackage->completed = false;
			sMetaPackage->served = true;
			sMetaPackage->request = MPI_REQUEST_NULL;
			int ns = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].ns;
			int nd = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].nd;
			int r_gid = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].r_gid;
			int r_grids = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].r_grids;
			//int tag = tagGen(ns, nd, r_gid-1, np*r_grids, nGraphs);
			int tag = Perilla::myTagMap[r_gid][nd][ns][graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].sz];
			sMetaPackage->request = ParallelDescriptor::Asend(sPackage->databuf,
				graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].sz,
				graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].pr, tag).req();  // tag == SeqNum in c++ ver
			//if(g == 31 && nd == 519 )
			//std::cout << "S Posted r_g " << r_gid << " atP " << myProc << " nd "<< nd << " ns "<<ns << " tag "<<tag << " pr " <<graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].pr << std::endl;

		    }
		}		
	    } // for (i<i<cpSrc->r_con.nsnd)	    
	    cpSrc = cpSrc->next;
	} // while(cpSrc != 0)	
    } // for(f<nfabs)

    for(int f=0; f<numfabs; f++)
    {	
	FabCopyAssoc* cpDst = graphArray[g]->task[f]->cpAsc_dstHead;
	while(cpDst != 0)
	{
	    for(int i=0; i<cpDst->r_con.nrcv; i++)
	    {
		if(graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pQueue.queueSize(true) > 0) //!all messages before rear have completed
		{		    
		    if(omp_test_lock(cpDst->r_con.rcvLock) != 0)
		    {		    
			Package *rearPackage =  graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pQueue.getRear(true);
			if(!rearPackage->completed)
			{
			    bool flag = false;
			    int ret_flag=0;
			    MPI_Status status;
			    ParallelDescriptor::Test(rearPackage->request, ret_flag, status);

			    flag = (ret_flag == 0) ? false : true;//parallel_test_one(rearPackage%ptr%request) -------???????
			    if(flag)
			    {
				rearPackage->completeRequest();				
				cpDst->r_con.rcv[i].pQueue.getRear()->completeRequest();

				if(graphArray[g]->rCopyMapHead->map[f]->r_con.rcv[i].pQueue.queueSize(true) == 1)
				{
				    //if(graphArray[g]->graphID == 25 && f==0 && myProc==1)
				    //std::cout<<"Recieved fc++ for f " << f << " fc " << cpDst->r_con.firingRuleCnt <<std::endl;
				    cpDst->r_con.firingRuleCnt++;
				}
#pragma omp flush			    
			    }
			}		   		    
			omp_unset_lock(cpDst->r_con.rcvLock);
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
		if(graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].pQueue.queueSize(true) > 0)
		{
		    Package *frontPackage = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].pQueue.getFront(true);
		    if(frontPackage->served && !frontPackage->completed) //!latest receive request has NOT been completed
		    {
			bool flag = false;
			int ret_flag;
			MPI_Status status;
			ParallelDescriptor::Test(frontPackage->request, ret_flag, status);
			flag = (ret_flag == 0) ? false : true;//parallel_test_one(frontPackage%ptr%request) -------???????		    
			if(flag)
			{
			    omp_set_lock(graphArray[g]->sCopyMapHead->map[f]->r_con.sndLock);
			    frontPackage = graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].pQueue.dequeue(true);
			    frontPackage->completed = false;
			    frontPackage->served = false;
			    frontPackage->request = MPI_REQUEST_NULL;
			    frontPackage->notified = false;
			    graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].recycleQueue.enqueue(frontPackage,true);
			    omp_unset_lock(graphArray[g]->sCopyMapHead->map[f]->r_con.sndLock);
#pragma omp flush
			    omp_set_lock(cpSrc->r_con.sndLock);
			    frontPackage = cpSrc->r_con.snd[i].pQueue.dequeue(true);
			    frontPackage->completed = false;
			    frontPackage->served = false;
			    frontPackage->request = MPI_REQUEST_NULL;
			    cpSrc->r_con.snd[i].recycleQueue.enqueue(frontPackage,true);
			    omp_unset_lock(cpSrc->r_con.sndLock);			
			}
		    }
		} // if(queueSize > 0)				
	    } // for (i<i<cpSrc->r_con.nsnd)	    
	    cpSrc = cpSrc->next;
	} // while(cpSrc != 0)	
    } // for(f<nfabs)
} // serviceRemoteGridCopyRequests

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
		    omp_set_lock(graphArray[g]->sCopyMapHead->map[f]->r_con.sndLock);
		    for(int i=0; i<cpSrc->r_con.nsnd; i++)
		    {
			graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].pQueue.enqueue(graphArray[g]->sCopyMapHead->map[f]->r_con.snd[i].recycleQueue.dequeue(true),true);	    
		    }
		    omp_unset_lock(graphArray[g]->sCopyMapHead->map[f]->r_con.sndLock);
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
		if(omp_test_lock(graphArray[g]->rCopyMapHead->map[f]->r_con.rcvLock) != 0)
		{
		    if(omp_test_lock(cpDst->r_con.rcvLock) != 0)
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
				rcvMetaPackage->request = MPI_REQUEST_NULL;	  
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
			omp_unset_lock(cpDst->r_con.rcvLock);
		    } // if(ga locked)
		    omp_unset_lock(graphArray[g]->rCopyMapHead->map[f]->r_con.rcvLock);
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
