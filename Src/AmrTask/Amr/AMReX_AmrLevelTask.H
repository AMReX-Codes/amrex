#ifndef AMREX_AmrLevelAsync_H_
#define AMREX_AmrLevelAsync_H_

#include "AMReX_AmrLevel.H"

namespace amrex {

    class AmrLevelAsync: public AmrLevel
    {
	public:
	    friend class AsyncFillPatchIterator;
	    virtual void initPerilla (Real time) = 0;
	    virtual void finalizePerilla (Real time) = 0;
    };

    class AsyncFillPatchIterator :public MFIter
    {
	public:

	    friend class AmrLevelAsync;
	    friend class RGIter;

	    AsyncFillPatchIterator (AmrLevel& amrlevel,
		    MultiFab& leveldata,
		    int       boxGrow,
		    Real      time,
		    int       state_indx,
		    int       scomp,
		    int       ncomp,
		    int       f,
		    int iter,
		    int       tid);

	    void initFillPatch(int boxGrow,
		    int time,
		    int index,
		    int scomp,
		    int ncomp,
		    int iter);
	    static void  initialSend(amrex::Vector<amrex::AsyncFillPatchIterator*> afpi,
		    amrex::Vector<amrex::AsyncFillPatchIterator*> upper_afpi,
		    int  boxGrow,
		    Real time,
		    int  state_indx,
		    int  scomp,
		    int  ncomp,
		    int  iter,
		    int tid);

	    void PushOnly (int  boxGrow,
		    Real time,
		    int  state_indx,
		    int  scomp,
		    int  ncomp,
		    int f,
		    int tid,
		    unsigned char pushLevel,
		    int tf,
		    bool singleT=false);

	    void SendIntraLevel (RGIter& rgi,
		    int  boxGrow,
		    Real time,
		    int  state_indx,
		    int  scomp,
		    int  ncomp,
		    int  iter,
		    int f,
		    int tid,
		    bool singleT=false);

	    void SendIntraLevel (RGIter* rgi,
		    int  boxGrow,
		    Real time,
		    int  state_indx,
		    int  scomp,
		    int  ncomp,
		    int  iter,
		    int f,
		    int tid,
		    bool singleT=false);

	    void SendInterLevel (RGIter& rgi,
		    int  boxGrow,
		    Real time,
		    int  state_indx,
		    int  scomp,
		    int  ncomp,
		    int  iter,
		    int f,
		    int tid,
		    bool singleT=false);

	    void SendInterLevel (RGIter* rgi,
		    int  boxGrow,
		    Real time,
		    int  state_indx,
		    int  scomp,
		    int  ncomp,
		    int  iter,
		    int f,
		    int tid,
		    bool singleT=false);

	    void Receive (RGIter& rgi,
		    int  boxGrow,
		    Real time,
		    int  state_indx,
		    int  scomp,
		    int  ncomp,
		    int f,
		    int tid,
		    bool singleT=false);

	    void Receive (RGIter* rgi,
		    int  boxGrow,
		    Real time,
		    int  state_indx,
		    int  scomp,
		    int  ncomp,
		    int f,
		    int tid,
		    bool singleT=false);

	    void Receive (RGIter& rgi,
		    MultiFab& dest,
		    int  boxGrow,
		    Real time,
		    int  state_indx,
		    int  scomp,
		    int  ncomp,
		    int f,
		    int tid,
		    bool singleT=false);

	    void Receive (RGIter* rgi,
		    MultiFab& dest,
		    int  boxGrow,
		    Real time,
		    int  state_indx,
		    int  scomp,
		    int  ncomp,
		    int f,
		    int tid,
		    bool singleT=false);

	    void PullOnly (int  boxGrow,
		    Real time,
		    int  state_indx,
		    int  scomp,
		    int  ncomp,
		    int f,
		    int tid,
		    bool singleT=false);

	    void PullOnly (MultiFab& dest,
		    int  boxGrow,
		    Real time,
		    int  state_indx,
		    int  scomp,
		    int  ncomp,
		    int f,
		    int tid,
		    bool singleT=false);
	    void FillFromTwoLevelsPush (Real time,
		    int index,
		    int scomp,
		    int dcomp,
		    int ncomp,
		    int f,
		    int tid,
		    unsigned char pushLevel,
		    int tf,
		    bool singleT);
	    void FillFromTwoLevelsPull (Real time,
		    int index,
		    int scomp,
		    int dcomp,
		    int ncomp,
		    int f,
		    int tid,
		    bool singleT);

	    ~AsyncFillPatchIterator ();

	    FArrayBox& operator() () { return m_fabs[MFIter::index()]; }

	    Box UngrownBox () const { return MFIter::validbox(); }

	    MultiFab& get_mf() { return m_fabs; }
	    AsyncFillPatchIterator ();
	    AsyncFillPatchIterator (const AsyncFillPatchIterator& rhs);
	    AsyncFillPatchIterator& operator= (const AsyncFillPatchIterator& rhs);

	    AmrLevel&                         m_amrlevel;
	    MultiFab&                         m_leveldata;
	    std::vector< std::pair<int,int> > m_range;
	    MultiFab                          m_fabs;
	    int                               m_ncomp;

	public:
	    bool                              isProperlyNested;

	    Vector<MultiFab*>                  smf;
	    Vector<Real>                       stime;
	    StateDataPhysBCFunct*             physbcf;
	    Geometry*                         geom;


	    Vector<MultiFab*>                  smf_crse;
	    Vector<Real>                       stime_crse;
	    StateDataPhysBCFunct*             physbcf_crse;
	    Geometry*                         geom_crse;

	    Vector<MultiFab*>                  smf_fine;
	    Vector<Real>                       stime_fine;
	    StateDataPhysBCFunct*             physbcf_fine;
	    Geometry*                         geom_fine;


	    RegionGraph*                      destGraph;
	    RegionGraph*                      csrcGraph;
	    RegionGraph*                      fsrcGraph;

	    MultiFab*                         m_mf_crse_patch;
	    RegionGraph*                      m_rg_crse_patch;
	    const FabArrayBase::FPinfo*       m_fpc;

	    //PArray<MultiFab>                  raii;
	    MultiFab*                         dmf;
	    MultiFab*                         dmff;

	    void completeRegionGraphs(int tg)
	    {
		//std::cout << "Completing RGs ";

		if(destGraph != NULL)
		{
		    //std::cout << destGraph->graphID << " ";
		    destGraph->completeRegionGraph(tg);
		}
		if(csrcGraph != NULL)
		{
		    //std::cout << csrcGraph->graphID << " ";
		    csrcGraph->completeRegionGraph(tg);
		}
		if(fsrcGraph != NULL)
		{
		    //std::cout << fsrcGraph->graphID << " ";
		    fsrcGraph->completeRegionGraph(tg);
		}
		if(m_rg_crse_patch != NULL)
		{
		    //std::cout << m_rg_crse_patch->graphID << " ";
		    m_rg_crse_patch->completeRegionGraph(tg);
		}

		//std::cout <<" by tg " << tg << std::endl;
	    }
	    void Reset(int tg)
	    {
		//std::cout << "Resetting RGs ";
		if(destGraph != NULL)
		{
		    //std::cout << destGraph->graphID << " ";
		    destGraph->Reset(tg);
		}
		if(csrcGraph != NULL)
		{
		    //std::cout << csrcGraph->graphID << " ";
		    csrcGraph->Reset(tg);
		}
		if(fsrcGraph != NULL)
		{
		    //std::cout << fsrcGraph->graphID << " ";
		    fsrcGraph->Reset(tg);
		}
		if(m_rg_crse_patch != NULL)
		{
		    //std::cout << m_rg_crse_patch->graphID << " ";
		    m_rg_crse_patch->Reset(tg);
		}
		//std::cout <<" by tg " << tg << std::endl;
	    }

	    // Variables for optimization calls of two level push/pulll

    };


#endif
