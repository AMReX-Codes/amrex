#ifndef AMREX_ASYNC_MFITER
#define AMREX_ASYNC_MFITER
//Question? email tannguyen@lbl.gov
//Created 09-07-2017
//Last modification 09-07-2017

#include "AMReX_AbstractTask.H"
#include "AMReX_TaskGraph.H"
#include "RTS.H"
#include <AMReX_FabArrayBase.H>
#include <AMReX_FabArray.H>
#include <AMReX_IntVect.H>
#include <AMReX_FArrayBox.H>
#include "AMReX_Connections.H"
#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>

//#ifdef _OPENMP
#include <omp.h>
//#endif

namespace amrex {
    typedef MFIter LocalFabIdx; 
    class Action :public Task{
	protected:
	    LocalConnection l_con;
	    RemoteConnection r_con;
	    struct TileArray
	    {
		Vector<int> numLocalTiles;
		Vector<int> localTileIndexMap;
		Vector<Box> tileArray;
	    } ta;
	    bool _do_tiling;
	    int current_tile;

	    size_t tagGen(int src, int dest, int amrlevels, int fabs)
	    {
		return (src*fabs + dest%fabs)*amrlevels;
	    }
	    void FillBoundary_Push(){
		for(int i=0; i<l_con.nscpy; i++)
		{
		    Data* msg= new Data(TaskName(l_con.scpy[i].ns), TaskName(l_con.scpy[i].nd), l_con.scpy[i].sz);
		    msg->SetTag(i); //this is local copy so the recipient will have the matching tag
		    msg->SetDestRank(ParallelDescriptor::MyProc());
		    _fab->copyToMem(l_con.scpy[i].sbx, 0, _mf->nComp(), msg->GetBuffer()); //pack box to the message
		    _outputs.push(msg);//let the runtime know that the output is available
		}
		int np = ParallelDescriptor::NProcs();
		if (np==1) return;
		for(int i=0; i<r_con.nsnd; i++)
		{
		    Data* msg= new Data(TaskName(r_con.snd[i].ns), TaskName(r_con.snd[i].nd), r_con.snd[i].sz);
		    msg->SetTag(tagGen(r_con.snd[i].ns, r_con.snd[i].nd, 1, _mf->size()));
		    msg->SetDestRank(r_con.snd[i].pr);
		    _fab->copyToMem(r_con.snd[i].sbx, 0, _mf->nComp(), msg->GetBuffer());
		    _outputs.push(msg);
		}
	    }
	    void FillBoundary_Pull(){
		for(int i=0; i<l_con.ndcpy; i++)
		{
		    int tag=l_con.dcpy[i].sPartner;
		    Data* msg= _neighbors_in.pop_front(TaskName(l_con.dcpy[i].ns), tag);
		    _fab->copyFromMem(l_con.dcpy[i].dbx,0, _mf->nComp(), msg->GetBuffer());//unpack message to box
		    msg->Free();
		}
		int np = ParallelDescriptor::NProcs();
		if (np==1) return;
		for(int i=0; i<r_con.nrcv; i++)
		{
		    size_t tag= tagGen(r_con.rcv[i].ns, r_con.rcv[i].nd, 1, _mf->size());
		    Data* msg= _neighbors_in.pop_front(TaskName(r_con.rcv[i].ns), tag);
		    _fab->copyFromMem(r_con.rcv[i].dbx, 0, _mf->nComp(), msg->GetBuffer());
		    msg->Free();
		}
	    }
	    bool isSatisfied(){
		if(!DependSignal()) return false; //wait for dependent signal if ANY
		bool satisfied=true;
		for(int i=0; i<l_con.ndcpy; i++){
		    satisfied= Depend_on(TaskName(l_con.dcpy[i].ns), l_con.dcpy[i].sPartner);
		    if(!satisfied) return false;
		}
		//remote connections
		if(ParallelDescriptor::NProcs() == 1) return true;
		for(int i=0; i<r_con.nrcv; i++)
		{
		    satisfied= Depend_on(TaskName(r_con.rcv[i].ns), tagGen(r_con.rcv[i].ns, r_con.rcv[i].nd, 1, _mf->size()));
		    if(!satisfied) return false;
		}
		return true;
	    }
	    void extendIters(int extra){_nIters+= extra;}
	    int _iter, _nIters;
	    FabArray<FArrayBox> *_mf;
	    FArrayBox* _fab;
	    int _idx;
	    int _lIdx;
	    bool _communicateFirstTimeStep; //exchange ghost cells before starting the first time step
	    bool _communicateUponCompletion; //exchange ghost cells after computing the last time step
	public:
	    LocalConnection& LCon(){return l_con;}
	    RemoteConnection& RCon(){return r_con;}
	    Action(){
		_iter=-1;
		_nIters=1;
		_communicateFirstTimeStep=true; //the default is we exchange ghost cells before the first time step
		_communicateUponCompletion=false;
		_do_tiling=false;
	    }
	    ~Action(){
		free(l_con.scpy);
		free(l_con.dcpy);
		free(r_con.snd);
		free(r_con.rcv);
	    }
	    void enable_tiling(){_do_tiling=true;}
	    vector<Box>& getTileArray(){return ta.tileArray;}
	    void SetFab(FArrayBox* fab){_fab= fab;}
	    void SetSteps(int nIters){
		assert(nIters>=1);
		_nIters= nIters;
	    }
	    void SetMF(const FabArray<FArrayBox> &mf){
		_mf= (FabArray<FArrayBox>*)&mf;
	    }
	    void SetIdx(int idx){
		_idx=idx;
	    }
	    void SetLocalIdx(int lIdx){_lIdx= lIdx;}
	    virtual void Compute(Box)=0;
	    virtual void Init(){};
	    virtual bool DependSignal(){return true;}

	    void tileIndices(int &beginIndex, int &endIndex){
#ifdef _OPENMP
		int tid = omp_get_thread_num();
		int nthreads= omp_get_num_threads();
		int ntot= ta.tileArray.size();
		int nr   = ntot / nthreads;
		int nlft = ntot - nr * nthreads;
		if (tid < nlft) {  // get nr+1 items
		    beginIndex += tid * (nr + 1);
		    endIndex = beginIndex + nr + 1;
		} else {           // get nr items
		    beginIndex += tid * nr + nlft;
		    endIndex = beginIndex + nr;
		}
#else
		beginIndex=0;
		endIndex= ta.tileArray.size();
#endif
	    }
	    void Job(){
		if(_communicateFirstTimeStep){
		    if(_iter==-1) FillBoundary_Push();
		    else if (_iter==0){
			FillBoundary_Pull();
		    }
		}else if(_iter==-1) _iter++; //go directly to the first compute step
		if(_iter>=0 && _iter<_nIters){ //always compute from time step 0 to _nIters-1
		    if(_iter>0)FillBoundary_Pull(); //communication at step 0 is already governed by _communicateFirstTimeStep 
		    if(!_do_tiling) Compute(validbox());//execute task at Fab level
		    else{
#ifdef _OPENMP
#pragma omp parallel
			{
#endif
			    int beginIndex, endIndex;
			    tileIndices(beginIndex, endIndex);
			    for(int tile=beginIndex; tile<endIndex; tile++){
				Compute(ta.tileArray[tile]);
			    }
#ifdef _OPENMP
			}
#endif
		    }
		    if(_iter<_nIters-1) FillBoundary_Push();
		}
		_iter++;
		if(_communicateUponCompletion){
		    if(_iter<=_nIters) {
			if(_iter==_nIters)FillBoundary_Push(); 
		    }else if (_iter==(_nIters+1)){
			FillBoundary_Pull();
			SelfDestroy();
		    }
		}else{
		    if(_iter==_nIters)SelfDestroy();
		}
	    }
	    bool Dependency(){
		if(_iter==-1) return true;
		return isSatisfied();
	    }
	    void SetFirstTimeStepComm(bool input){
		_communicateFirstTimeStep= input;
	    }
	    void SetCompletionComm(bool input){
		_communicateUponCompletion= input;
	    }
	    void PostCompletion(){
		//nothing
	    }
	    Box validbox() const{
		return _mf->box(_idx);
	    }
	    //! Return the Fab associated with this task
	    FArrayBox& validFab(){
		return *(_mf->m_fabs_v[_lIdx]);
	    }
	    //! Return the local index of the Fab associated with this task. It can then be used to locate corresponding Fabs on other multifabs that have the same layout.
	    int localFabIdx(){
		return _lIdx;
	    }
	    //! Locate in multifab mf the Fab that has the same coordinate as this task' Fab
	    FArrayBox& validFab(const FabArray<FArrayBox> &mf){
		return *(mf.m_fabs_v[_lIdx]);
	    }
	    //! Locate in multifab mf the Fab that has local index lIdx
	    FArrayBox& validFab(const FabArray<FArrayBox> &mf, int lIdx){
		return *(mf.m_fabs_v[lIdx]);
	    }

	    void generateTileArray(const IntVect& tileSize){
		IntVect nt_in_fab, tsize, nleft;
		int ntiles = 1;
		const Box& bx = validbox();

		for (int d=0; d<BL_SPACEDIM; d++) {
		    int ncells = bx.length(d);
		    nt_in_fab[d] = std::max(ncells/tileSize[d], 1);
		    tsize    [d] = ncells/nt_in_fab[d];
		    nleft    [d] = ncells - nt_in_fab[d]*tsize[d];
		    ntiles *= nt_in_fab[d];
		}

		IntVect small, big, ijk;  // note that the initial values are all zero.
		ijk[0] = -1;
		for (int t = 0; t < ntiles; ++t) {
		    ta.localTileIndexMap.push_back(t);
		    ta.numLocalTiles.push_back(ntiles);

		    for (int d=0; d<BL_SPACEDIM; d++) {
			if (ijk[d]<nt_in_fab[d]-1) {
			    ijk[d]++;
			    break;
			} else {
			    ijk[d] = 0;
			}
		    }

		    for (int d=0; d<BL_SPACEDIM; d++) {
			if (ijk[d] < nleft[d]) {
			    small[d] = ijk[d]*(tsize[d]+1);
			    big[d] = small[d] + tsize[d];
			} else {
			    small[d] = ijk[d]*tsize[d] + nleft[d];
			    big[d] = small[d] + tsize[d] - 1;
			}
		    }

		    Box tbx(small, big, IndexType::TheCellType());
		    tbx.shift(bx.smallEnd());

		    ta.tileArray.push_back(tbx);
		}
	    }
	    Box growntilebox (Box bx, int ng){
		const Box& vbx = _mf->m_fabs_v[_lIdx]->box();
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
    };//end class Action

    template <class T>
	class MFGraph: public AbstractTaskGraph<T>{
	    protected:
		string _graphName;
		bool _do_tiling;
	    public:
		MFGraph(const FabArray<FArrayBox> &mf, int nSteps, int rank, int nProcs, Periodicity period, bool do_tiling){
		    AbstractTaskGraph<T>::_nProcs= nProcs;
		    AbstractTaskGraph<T>::_rank= rank;
		    for(int i=0; i< mf.local_size(); i++){
			TaskName name(mf.IndexArray()[i]);
			T *t= new T();
			t->SetSteps(nSteps);
			t->SetMF(mf);
			t->SetFab(mf.m_fabs_v[i]);
			t->SetIdx(mf.IndexArray()[i]);
			t->SetName(name);
			t->SetLocalIdx(i);
			if(do_tiling){  
			    IntVect ts= FabArrayBase::mfiter_tile_size;
			    t->generateTileArray(ts); //create tile array associated with this FAB
			}
			t->Init();
			if(do_tiling) t->enable_tiling();
			AbstractTaskGraph<T>::_initialTasks.push_back(t);
			AbstractTaskGraph<T>::_taskPool[name]= t;
		    }
		    AbstractTaskGraph<T>::_begin= *(AbstractTaskGraph<T>::_initialTasks.begin());
		    AbstractTaskGraph<T>::_end= *(AbstractTaskGraph<T>::_initialTasks.end());
		    AbstractTaskGraph<T>::_currIt= AbstractTaskGraph<T>::_initialTasks.begin();
		    AbstractTaskGraph<T>::_current= *(AbstractTaskGraph<T>::_currIt);
		    AbstractTaskGraph<T>::_mode= _Push;
		    SetupFabConnections(mf, period);
		    _do_tiling= do_tiling;
		}

		MFGraph(const Amr* amr, int max_step, Real stop_time, int rank, int nProcs, bool do_tiling){
		    AbstractTaskGraph<T>::_nProcs= nProcs;
		    AbstractTaskGraph<T>::_rank= rank;
		    //create an initial graph corresponding to the coarsest AMR level, this graph will evolve with time
		    int n_intervals=1;
		    for(int l=0; l<= amr->finest_level; l++){
			n_intervals*= amr->n_cycle[l];
			for(int s=0; s< amr->amr_level[l]->numStates(); s++){
			    MultiFab& mf = amr->amr_level[l]->get_new_data(s);
			    amr->amr_level[l]->get_state_data(s).allocOldData();
			    for(int i=0; i< mf.local_size(); i++){
				for(int it=0; it< n_intervals; it++){
				    TaskName name(l /*level*/, s/*state*/, mf.IndexArray()[i] /*fab*/, it);
				    T *t= new T();
				    t->SetName(name);
				    t->SetLocalIdx(i);
				    t->SetMF(mf);
				    t->SetFab(mf.m_fabs_v[i]);
				    t->SetIdx(mf.IndexArray()[i]);
				    t->InitAmrTask(amr, max_step, stop_time);
				    t->CreateLevelTask(l);
				    if(i==0) t->SetMaster();
				    if(do_tiling){
					IntVect ts= FabArrayBase::mfiter_tile_size;
					t->generateTileArray(ts); //create tile array associated with this FAB
				    }
				    if(do_tiling) t->enable_tiling();
				    AbstractTaskGraph<T>::_initialTasks.push_back(t);
				    AbstractTaskGraph<T>::_taskPool[name]= t;
		    		    SetupFabConnections(mf, amr->amr_level[l]->Geom().periodicity());
				}
			    }
			}
		    }
		    AbstractTaskGraph<T>::_begin= *(AbstractTaskGraph<T>::_initialTasks.begin());
		    AbstractTaskGraph<T>::_end= *(AbstractTaskGraph<T>::_initialTasks.end());
		    AbstractTaskGraph<T>::_currIt= AbstractTaskGraph<T>::_initialTasks.begin();
		    AbstractTaskGraph<T>::_current= *(AbstractTaskGraph<T>::_currIt);
		    AbstractTaskGraph<T>::_mode= _Push;
		    _do_tiling= do_tiling;
		}

		int FindProcessAssociation(TaskName name){
		    assert(false);
		}
		void SetupFabConnections(const FabArray<FArrayBox> &mf, Periodicity period){
		    int np = ParallelDescriptor::NProcs();
		    int myProc = ParallelDescriptor::MyProc();
		    int numfabs = mf.size();
		    bool cross = false;
		    const FabArrayBase::FB& TheFB = mf.getFB(period);
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
		    for(int f=0; f<numfabs; f++) //          !create local communication metadata for each fab
		    {
			LocalConnection &l_con= ((Action*)AbstractTaskGraph<T>::_initialTasks[f])->LCon();
			l_con.nscpy = 0;
			l_con.ndcpy = 0;
			l_con.firingRuleCnt = 0;
			for(int i=0; i<n_loc_mf; i++)
			{
			    const FabArrayBase::CopyComTag& tag = (*TheFB.m_LocTags)[i];
			    if(f == tag.srcIndex)
				l_con.nscpy++;
			    if(f == tag.dstIndex)
				l_con.ndcpy++;
			}
		    }
		    if(np > 1)
		    {
			for(int f=0; f<numfabs; f++)
			{
			    RemoteConnection &r_con= ((Action*)AbstractTaskGraph<T>::_initialTasks[f])->RCon();
			    r_con.nrcv = 0;
			    r_con.nsnd = 0;
			    r_con.firingRuleCnt = 0;
			    for(int i=0; i<n_snds_mf; i++)
			    {
				const FabArrayBase::CopyComTagsContainer& cctc = *send_cctc[i];
				for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
					it != cctc.end(); ++it)
				{
				    if(f == it->srcIndex)
					r_con.nsnd++;
				}
			    }
			    for(int i=0; i<n_rcvs_mf; i++)
			    {
				const FabArrayBase::CopyComTagsContainer& cctc = *recv_cctc[i];
				for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
					it != cctc.end(); ++it)
				{
				    if(f == it->dstIndex)
					r_con.nrcv++;
				}
			    }
			}
		    }
		    int scnt, dcnt;
		    for(int f=0; f<numfabs; f++)
		    {
			scnt = -1;
			dcnt = -1;
			LocalConnection &l_con= ((Action*)AbstractTaskGraph<T>::_initialTasks[f])->LCon();
			l_con.scpy = new LocalCopyDescriptor[l_con.nscpy];
			l_con.dcpy = new LocalCopyDescriptor[l_con.ndcpy];
			for(int i=0; i<n_loc_mf; i++)
			{
			    const FabArrayBase::CopyComTag& tag = (*TheFB.m_LocTags)[i];
			    if(f == tag.srcIndex)
			    {
				scnt++;
				l_con.scpy[scnt].ns = tag.srcIndex;
				l_con.scpy[scnt].nd = tag.dstIndex;
				l_con.scpy[scnt].sbx = tag.sbox;
				l_con.scpy[scnt].dbx = tag.dbox;
				l_con.scpy[scnt].sz = tag.sbox.numPts() * mf.nComp() * sizeof(double);
			    }
			    if(f == tag.dstIndex)
			    {
				dcnt++;
				l_con.dcpy[dcnt].ns = tag.srcIndex;
				l_con.dcpy[dcnt].nd = tag.dstIndex;
				l_con.dcpy[dcnt].sbx = tag.sbox;
				l_con.dcpy[dcnt].dbx = tag.dbox; 
				l_con.dcpy[dcnt].sz = tag.dbox.numPts() * mf.nComp() * sizeof(double);
			    }
			} 
		    }
		    for(int f=0; f<numfabs; f++)
		    {
			LocalConnection &l_con= ((Action*)AbstractTaskGraph<T>::_initialTasks[f])->LCon();
			for(int i=0; i<l_con.nscpy; i++){
			    LocalConnection &l_con1= ((Action*)AbstractTaskGraph<T>::_initialTasks[l_con.scpy[i].nd])->LCon();
			    for(int j=0; j<l_con1.ndcpy; j++){
				if(l_con.scpy[i].dbx == l_con1.dcpy[j].dbx)
				    l_con.scpy[i].dPartner = j;
			    }
			}

			for(int i=0; i<l_con.ndcpy; i++){
			    LocalConnection &l_con1= ((Action*)AbstractTaskGraph<T>::_initialTasks[l_con.dcpy[i].ns])->LCon();
			    for(int j=0; j<l_con1.nscpy; j++)
				if(l_con.dcpy[i].dbx == l_con1.scpy[j].dbx)
				    l_con.dcpy[i].sPartner = j;
			}
		    }


		    if(np == 1) return;
		    int nsnd, nrcv;

		    for(int f=0; f<numfabs; f++)
		    {
			RemoteConnection &r_con= ((Action*)AbstractTaskGraph<T>::_initialTasks[f])->RCon();
			r_con.snd = new RemoteCommDescriptor[r_con.nsnd];
			r_con.rcv = new RemoteCommDescriptor[r_con.nrcv];
			nrcv= -1;
			for(int i=0; i<n_rcvs_mf; i++)
			{
			    const FabArrayBase::CopyComTagsContainer& cctc = *recv_cctc[i];
			    int pr = recv_pr[i];
			    for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
				    it != cctc.end(); ++it)
			    {
				if(f == it->dstIndex)
				{
				    nrcv++;
				    r_con.rcv[nrcv].ns = it->srcIndex;
				    r_con.rcv[nrcv].lns = -1;
				    r_con.rcv[nrcv].nd = it->dstIndex;
				    r_con.rcv[nrcv].lnd = mf.localindex(it->dstIndex);
				    r_con.rcv[nrcv].sbx = it->sbox;
				    r_con.rcv[nrcv].dbx = it->dbox;
				    r_con.rcv[nrcv].pr = pr;
				    r_con.rcv[nrcv].cnt = 0;
				    r_con.rcv[nrcv].sz = it->sbox.numPts() * mf.nComp() * sizeof(double);
				}
			    }
			}
			nsnd = -1;
			for(int i=0; i<n_snds_mf; i++)
			{
			    const FabArrayBase::CopyComTagsContainer& cctc = *send_cctc[i];
			    int pr = send_pr[i];
			    for (FabArrayBase::CopyComTagsContainer::const_iterator it = cctc.begin();
				    it != cctc.end(); ++it)
			    {
				if(f == it->dstIndex )
				{
				    nsnd++;
				    r_con.snd[nsnd].ns = it->srcIndex;
				    r_con.snd[nsnd].lns = mf.localindex(it->srcIndex);
				    r_con.snd[nsnd].nd = it->dstIndex;
				    r_con.snd[nsnd].lnd = -1;
				    r_con.snd[nsnd].sbx = it->sbox;
				    r_con.snd[nsnd].dbx = it->dbox;
				    r_con.snd[nsnd].pr = pr;
				    r_con.snd[nsnd].cnt = 0;
				    r_con.snd[nsnd].sz = it->dbox.numPts() * mf.nComp() * sizeof(double);
				}
			    }
			} // for(i<n_snds_mf)
		    }
		}
	};

    template <class A>
#define max_life 2^32-1
	class AMFIter{
	    private:
		RTS rts;
		MFGraph<A> *graph;
		bool _do_tiling;
	    public:
		AMFIter(const FabArray<FArrayBox> &mf, int nSteps, Periodicity period, bool do_tiling=false){
		    _do_tiling= do_tiling;
		    graph= new MFGraph<A>(mf, nSteps, ParallelDescriptor::MyProc(), ParallelDescriptor::NProcs(), period, do_tiling);
		}

		AMFIter(const FabArray<FArrayBox> &mf, Periodicity period, bool do_tiling=false){
		    _do_tiling= do_tiling;
		    graph= new MFGraph<A>(mf, max_life, ParallelDescriptor::MyProc(), ParallelDescriptor::NProcs(), period, do_tiling);
		}

		AMFIter(const Amr* amr, int max_step, Real stop_time, bool do_tiling=false){
		    _do_tiling= do_tiling;
		    graph= new MFGraph<A>(amr, max_step, stop_time, ParallelDescriptor::MyProc(), ParallelDescriptor::NProcs(), do_tiling);
		}

		~AMFIter(){
		    delete graph;
		}

		void Iterate(){
		    rts.Init(ParallelDescriptor::MyProc(), ParallelDescriptor::NProcs());
		    rts.Iterate(graph);
		    rts.Finalize();
		}
	};
#undef max_life
} //end namespace
#endif
