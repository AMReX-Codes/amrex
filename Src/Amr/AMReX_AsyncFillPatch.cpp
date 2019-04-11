#include <sstream>
#include <unistd.h>
#include <memory>
#include <limits>

#include <AMReX_AmrLevel.H>
#include <AMReX_Derive.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>

#ifdef USE_PERILLA
#include <WorkerThread.H>

namespace amrex {
    using namespace perilla;

    AsyncFillPatchIterator::AsyncFillPatchIterator (AmrLevel& amrlevel,
	    MultiFab& leveldata,
	    int       boxGrow,
	    Real      time,
	    int       index,
	    int       scomp,
	    int       ncomp,
	    int       iter)
	:
	    MFIter(leveldata),
	    m_amrlevel(amrlevel),
	    m_leveldata(leveldata),
	    m_ncomp(ncomp),
	    physbcf(0),
	    physbcf_crse(0),
	    physbcf_fine(0),
	    destGraph(0),
	    fsrcGraph(0),
	    csrcGraph(0),
	    m_rg_crse_patch(0),
	    dmf(0),
	    dmff(0)
    {
	BL_ASSERT(scomp >= 0);
	BL_ASSERT(ncomp >= 1);
	BL_ASSERT(AmrLevel::desc_lst[index].inRange(scomp,ncomp));
	BL_ASSERT(0 <= index && index < AmrLevel::desc_lst.size());

	initFillPatch(boxGrow, time, index, scomp, ncomp, iter);

#ifdef BL_USE_TEAM
	ParallelDescriptor::MyTeam().MemoryBarrier();
#endif
    }


AsyncFillPatchIterator::~AsyncFillPatchIterator () {
#ifdef USE_PERILLA
        while(regionList.size()){
          RegionGraph* tmp= regionList.front();
          delete tmp;
          regionList.pop_front();
        }

        while(mfList.size()){
          MultiFab *tmp= mfList.front();
          delete tmp;
          mfList.pop_front();
        }

        while(stateDataList.size()){
          StateDataPhysBCFunct *tmp= stateDataList.front();
          delete tmp;
          stateDataList.pop_front();
        }
#endif
    }


    void AsyncFillPatchIterator::FillFromTwoLevelsPush (Real time,
	    int index,
	    int scomp,
	    int dcomp,
	    int ncomp,
	    int f,
	    unsigned char pushLevel,
	    bool singleT)
    {

	int ilev_fine = m_amrlevel.level;
	int ilev_crse = ilev_fine-1;

	BL_ASSERT(ilev_crse >= 0);

	AmrLevel& fine_level = m_amrlevel;
	AmrLevel& crse_level = m_amrlevel.parent->getLevel(ilev_crse);

	Geometry* tgeom_fine = &fine_level.geom;
	Geometry* tgeom_crse = &crse_level.geom;

	Vector<MultiFab*> tsmf_crse;
	Vector<MultiFab*> tsmf_fine;
	Vector<Real> tstime_crse;
	Vector<Real> tstime_fine;
	StateData& statedata_crse = crse_level.state[index];
	statedata_crse.getData(tsmf_crse,tstime_crse,time);
	StateDataPhysBCFunct* tphysbcf_crse = new StateDataPhysBCFunct(statedata_crse,scomp,*geom_crse);

	StateData& statedata_fine = fine_level.state[index];
	statedata_fine.getData(tsmf_fine,tstime_fine,time);
	StateDataPhysBCFunct* tphysbcf_fine = new StateDataPhysBCFunct(statedata_fine,scomp,*geom_fine);

	const StateDescriptor& desc = AmrLevel::desc_lst[index];

	FillPatchTwoLevelsPush(*(m_amrlevel.parent), m_fabs, time,
		tsmf_crse, tstime_crse,
		tsmf_fine, tstime_fine,
		destGraph, csrcGraph, fsrcGraph, f,
		this,
		dmf,
		dmff,
		scomp, dcomp, ncomp,
		*tgeom_crse, *tgeom_fine,
		*tphysbcf_crse, *tphysbcf_fine,
		crse_level.fineRatio(),
		desc.interp(scomp), desc.getBCs(), pushLevel,singleT);
    }

    void AsyncFillPatchIterator::FillFromTwoLevelsPull (Real time,
	    int index,
	    int scomp,
	    int dcomp,
	    int ncomp,
	    int f,
	    bool singleT)
    {

	int ilev_fine = m_amrlevel.level;
	int ilev_crse = ilev_fine-1;

	BL_ASSERT(ilev_crse >= 0);

	AmrLevel& fine_level = m_amrlevel;
	AmrLevel& crse_level = m_amrlevel.parent->getLevel(ilev_crse);

	Geometry* tgeom_fine = &fine_level.geom;
	Geometry* tgeom_crse = &crse_level.geom;

	Vector<MultiFab*> tsmf_crse;
	Vector<MultiFab*> tsmf_fine;
	Vector<Real> tstime_crse;
	Vector<Real> tstime_fine;
	StateData& statedata_crse = crse_level.state[index];
	statedata_crse.getData(tsmf_crse,tstime_crse,time);
	StateDataPhysBCFunct* tphysbcf_crse = new StateDataPhysBCFunct(statedata_crse,scomp,*geom_crse);

	StateData& statedata_fine = fine_level.state[index];
	statedata_fine.getData(tsmf_fine,tstime_fine,time);
	StateDataPhysBCFunct* tphysbcf_fine = new StateDataPhysBCFunct(statedata_fine,scomp,*geom_fine);


	const StateDescriptor& desc = AmrLevel::desc_lst[index];

	FillPatchTwoLevelsPull(m_fabs, time,
		tsmf_crse, tstime_crse,
		tsmf_fine, tstime_fine,
		destGraph, csrcGraph, fsrcGraph, f,
		this,
		scomp, dcomp, ncomp,
		*tgeom_crse, *tgeom_fine,
		*tphysbcf_crse, *tphysbcf_fine,
		crse_level.fineRatio(),
		desc.interp(scomp), desc.getBCs(), singleT);
    }

    void AsyncFillPatchIterator::FillFromTwoLevelsPull (MultiFab& dest,
	    Real time,
	    int index,
	    int scomp,
	    int dcomp,
	    int ncomp,
	    int f,
	    bool singleT)
    {
	int ilev_fine = m_amrlevel.level;
	int ilev_crse = ilev_fine-1;

	BL_ASSERT(ilev_crse >= 0);

	AmrLevel& fine_level = m_amrlevel;
	AmrLevel& crse_level = m_amrlevel.parent->getLevel(ilev_crse);

	Geometry* tgeom_fine = &fine_level.geom;
	Geometry* tgeom_crse = &crse_level.geom;

	Vector<MultiFab*> tsmf_crse;
	Vector<MultiFab*> tsmf_fine;
	Vector<Real> tstime_crse;
	Vector<Real> tstime_fine;
	StateData& statedata_crse = crse_level.state[index];
	statedata_crse.getData(tsmf_crse,tstime_crse,time);
	StateDataPhysBCFunct* tphysbcf_crse = new StateDataPhysBCFunct(statedata_crse,scomp,*geom_crse);

	StateData& statedata_fine = fine_level.state[index];
	statedata_fine.getData(tsmf_fine,tstime_fine,time);
	StateDataPhysBCFunct* tphysbcf_fine = new StateDataPhysBCFunct(statedata_fine,scomp,*geom_fine);


	const StateDescriptor& desc = AmrLevel::desc_lst[index];

	FillPatchTwoLevelsPull(dest, time,
		tsmf_crse, tstime_crse,
		tsmf_fine, tstime_fine,
		destGraph, csrcGraph, fsrcGraph, f,
		this,
		scomp, dcomp, ncomp,
		*tgeom_crse, *tgeom_fine,
		*tphysbcf_crse, *tphysbcf_fine,
		crse_level.fineRatio(),
		desc.interp(scomp), desc.getBCs(), singleT);
    }

    void AsyncFillPatchIterator::initFillPatch(int boxGrow,
	    Real time,
	    int index,
	    int scomp,
	    int ncomp,
	    int iter)
    {
	BL_ASSERT(scomp >= 0);
	BL_ASSERT(ncomp >= 1);
	BL_ASSERT(0 <= index && index < AmrLevel::desc_lst.size());

	int myProc = amrex::ParallelDescriptor::MyProc();

	const StateDescriptor& desc = AmrLevel::desc_lst[index];

	m_ncomp = ncomp;
	m_range = desc.sameInterps(scomp,ncomp);

	m_fabs.define(m_leveldata.boxArray(),m_leveldata.DistributionMap(),
		m_ncomp,boxGrow);

	BL_ASSERT(m_leveldata.DistributionMap() == m_fabs.DistributionMap());

	const IndexType& boxType = m_leveldata.boxArray().ixType();
	const int level = m_amrlevel.level;

	for (int i = 0, DComp = 0; i < m_range.size(); i++)
	{
	    const int SComp = m_range[i].first;
	    const int NComp = m_range[i].second;
	    int dcomp = DComp;

	    if (level == 0)
	    {
		BL_ASSERT(m_amrlevel.level == 0);
		StateData& statedata = m_amrlevel.state[index];
		statedata.getData(smf,stime,time);
		geom = &m_amrlevel.geom;
		physbcf = new StateDataPhysBCFunct(statedata,scomp,*geom);
                stateDataList.push_back(physbcf);
		BL_ASSERT(scomp+ncomp <= smf[0]->nComp());
		BL_ASSERT(dcomp+ncomp <= m_fabs.nComp());
		BL_ASSERT(smf.size() == stime.size());
		BL_ASSERT(smf.size() != 0);

		if (smf.size() == 1)
		{
		    dmf = new MultiFab(smf[0]->boxArray(), smf[0]->DistributionMap(), ncomp, 0);
		    destGraph = new RegionGraph(m_fabs.IndexArray().size());
		    fsrcGraph = new RegionGraph(smf[0]->IndexArray().size());
                    regionList.push_back(destGraph);
                    regionList.push_back(fsrcGraph);
		    Perilla::multifabExtractCopyAssoc( destGraph, fsrcGraph, m_fabs, *smf[0], ncomp, m_fabs.nGrow(), 0, geom->periodicity());
		    m_amrlevel.parent->graphArray[level].push_back(destGraph);
		    m_amrlevel.parent->graphArray[level].push_back(fsrcGraph);
		}
		else if (smf.size() == 2)
		{
		    BL_ASSERT(smf[0]->boxArray() == smf[1]->boxArray());
		    if (m_fabs.boxArray() == smf[0]->boxArray())
		    {
			dmf = &m_fabs;
			destGraph = new RegionGraph(m_fabs.IndexArray().size());
                        regionList.push_back(destGraph);
			Perilla::multifabBuildFabCon(destGraph, m_fabs, geom->periodicity());
			m_amrlevel.parent->graphArray[level].push_back(destGraph);
		    }
		    else
		    {
			dmf = new MultiFab(smf[0]->boxArray(), smf[0]->DistributionMap(), ncomp, 0);
			destGraph = new RegionGraph(m_fabs.IndexArray().size());
			fsrcGraph = new RegionGraph(dmf->IndexArray().size());
			fsrcGraph->buildTileArray(*dmf);
                        regionList.push_back(destGraph);
                        regionList.push_back(fsrcGraph);
                        mfList.push_back(dmf);

			Perilla::multifabExtractCopyAssoc( destGraph, fsrcGraph, m_fabs, *dmf, ncomp, m_fabs.nGrow(), 0, geom->periodicity());
			m_amrlevel.parent->graphArray[level].push_back(destGraph);
			m_amrlevel.parent->graphArray[level].push_back(fsrcGraph);
		    }
		}
		else
		{
		    amrex::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
		}
		//-------------------------------------------------- FillFromLevel0 initialization completed
	    }
	    else
	    {
		isProperlyNested = amrex::ProperlyNested(m_amrlevel.crse_ratio,
			m_amrlevel.parent->blockingFactor(m_amrlevel.level),
			boxGrow, boxType, desc.interp(SComp));
		if (level == 1 || isProperlyNested)
		{
		    int ilev_fine = m_amrlevel.level;
		    int ilev_crse = ilev_fine-1;
		    BL_ASSERT(ilev_crse >= 0);
		    AmrLevel& fine_level = m_amrlevel;
		    AmrLevel& crse_level = m_amrlevel.parent->getLevel(ilev_crse);
		    geom_fine = &fine_level.geom;
		    geom_crse = &crse_level.geom;
		    StateData& statedata_crse = crse_level.state[index];
		    statedata_crse.getData(smf_crse,stime_crse,time);
		    physbcf_crse = new StateDataPhysBCFunct(statedata_crse,scomp,*geom_crse);
		    StateData& statedata_fine = fine_level.state[index];
		    statedata_fine.getData(smf_fine,stime_fine,time);
		    physbcf_fine = new StateDataPhysBCFunct(statedata_fine,scomp,*geom_fine);

                    stateDataList.push_back(physbcf_crse);
                    stateDataList.push_back(physbcf_fine);

		    const StateDescriptor& desc = AmrLevel::desc_lst[index];
		    int ngrow = m_fabs.nGrow();
		    if (ngrow > 0 || m_fabs.getBDKey() != smf_fine[0]->getBDKey())
		    {
			InterpolaterBoxCoarsener coarsener = desc.interp(scomp)->BoxCoarsener(crse_level.fineRatio());
			Box fdomain = geom_fine->Domain();
			fdomain.convert(m_fabs.boxArray().ixType());
			Box fdomain_g(fdomain);
			for (int i = 0; i < BL_SPACEDIM; ++i) {
			    if (geom_fine->isPeriodic(i)) {
				fdomain_g.grow(i,ngrow);
			    }
			}
			// dummytostopcraycompilererror
			std::cout << "";

			Box c_dom= amrex::coarsen(geom_fine->Domain(), m_amrlevel.crse_ratio);

			m_fpc = &FabArrayBase::TheFPinfo(*smf_fine[0], m_fabs, fdomain_g, IntVect(ngrow), coarsener, c_dom);

			if (!m_fpc->ba_crse_patch.empty())
			{
			    m_mf_crse_patch = new MultiFab(m_fpc->ba_crse_patch, m_fpc->dm_crse_patch, ncomp, 0);
                            mfList.push_back(m_mf_crse_patch);
			    BL_ASSERT(scomp+ncomp <= smf_crse[0]->nComp());
			    BL_ASSERT(dcomp+ncomp <= m_mf_crse_patch->nComp());
			    BL_ASSERT(smf_crse.size() == stime_crse.size());
			    BL_ASSERT(smf_crse.size() != 0);

			    //if (smf_crse.size() == 1)
			    if (iter == 1)
			    {
				m_rg_crse_patch = new RegionGraph(m_mf_crse_patch->IndexArray().size());

				m_rg_crse_patch->isDepGraph = true;

				csrcGraph = new RegionGraph(smf_crse[0]->IndexArray().size());

                                regionList.push_back(m_rg_crse_patch);
                                regionList.push_back(csrcGraph);

				Perilla::multifabExtractCopyAssoc( m_rg_crse_patch, csrcGraph, *m_mf_crse_patch, *smf_crse[0], ncomp, m_mf_crse_patch->nGrow(), 0,geom_crse->periodicity());
#if 0
				MultiFab  temp_4_tile(m_fpc->ba_dst_boxes, m_fpc->dm_crse_patch, ncomp, 0);
				m_rg_crse_patch->buildTileArray(temp_4_tile);
#endif
				///m_rg_crse_patch->buildTileArray(*m_mf_crse_patch);
				m_amrlevel.parent->graphArray[level].push_back(m_rg_crse_patch);
				m_amrlevel.parent->graphArray[level].push_back(csrcGraph);
			    }

			    else if (iter > 1)
			    {
				if (m_mf_crse_patch->boxArray() == smf_crse[0]->boxArray())
				{
				    //dmf = m_mf_crse_patch;
				    m_rg_crse_patch = new RegionGraph(m_mf_crse_patch->IndexArray().size());

				    //std::cout<< " level " << level  << " rg_crs_ptch ID " << m_rg_crse_patch->graphID << std::endl;

				    Perilla::multifabBuildFabCon(m_rg_crse_patch, *m_mf_crse_patch,geom_crse->periodicity());
				    m_amrlevel.parent->graphArray[level].push_back(m_rg_crse_patch);
                                    regionList.push_back(m_rg_crse_patch);
				}
				else
				{
				    dmf = new MultiFab(smf_crse[0]->boxArray(), smf_crse[0]->DistributionMap(), ncomp, 0);
				    m_rg_crse_patch = new RegionGraph(m_mf_crse_patch->IndexArray().size());
				    m_rg_crse_patch->isDepGraph = true;
				    csrcGraph = new RegionGraph(dmf->IndexArray().size());
				    csrcGraph->buildTileArray(*dmf);

#if 0
				    MultiFab      temp_4_tile(m_fpc->ba_dst_boxes, m_fpc->dm_crse_patch, ncomp, 0);
				    m_rg_crse_patch->buildTileArray(temp_4_tile);
#endif

                                      regionList.push_back(m_rg_crse_patch);
                                      regionList.push_back(csrcGraph);
                                      mfList.push_back(dmf);


				    Perilla::multifabExtractCopyAssoc( m_rg_crse_patch, csrcGraph, *m_mf_crse_patch, *dmf, ncomp, m_mf_crse_patch->nGrow(), 0, geom_crse->periodicity());
				    m_amrlevel.parent->graphArray[level].push_back(m_rg_crse_patch);
				    m_amrlevel.parent->graphArray[level].push_back(csrcGraph);
				}
			    }
			    else
			    {
				amrex::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
			    }

			}
		    }

		    BL_ASSERT(scomp+ncomp <= smf_fine[0]->nComp());
		    BL_ASSERT(dcomp+ncomp <= m_fabs.nComp());
		    BL_ASSERT(smf_fine.size() == stime_fine.size());
		    BL_ASSERT(smf_fine.size() != 0);

		    //if (smf_fine.size() == 1)       // probabily it should aways be this because same level
		    //if (iter == 1)
		    if(true) // it will always be the case because same level comm and time will be available
		    {
			destGraph = new RegionGraph(m_fabs.IndexArray().size());
			fsrcGraph = new RegionGraph(smf_fine[0]->IndexArray().size());

                          regionList.push_back(destGraph);
                          regionList.push_back(fsrcGraph);


			if(m_rg_crse_patch != 0)
			{
			    destGraph->srcLinkGraph = m_rg_crse_patch;
			    //for(int lfi=0; lfi < destGraph->numTasks; lfi++ )

			    //std::cout << " m_mf_crse_patch->IndexArray().size() " << m_mf_crse_patch->IndexArray().size() << " size " << m_mf_crse_patch->size() << " myP " << myProc<< std::endl;

			    {
				for (MFIter mfi(*(m_mf_crse_patch),false); mfi.isValid(); ++mfi)
				{
				    int li = mfi.LocalIndex();
				    int gi = m_fpc->dst_idxs[li];
				    //if(gi == m_mf_crse_patch->IndexArray()[li])
				    {
					int lfi = m_fabs.localindex(gi);
					destGraph->task[lfi]->depTasksCompleted = false;
					destGraph->task[lfi]->depTaskIDs.push_back(li);

				    }
				}
			    }
			}
			//if(level == 2)
			//std::cout<< "Sending In "<<destGraph<<" "<< fsrcGraph << " myP " << ParallelDescriptor::MyProc()<<std::endl;                                      
			Perilla::multifabExtractCopyAssoc( destGraph, fsrcGraph, m_fabs, *smf_fine[0], ncomp, m_fabs.nGrow(), 0, geom_fine->periodicity());

			m_amrlevel.parent->graphArray[level].push_back(destGraph);
			m_amrlevel.parent->graphArray[level].push_back(fsrcGraph);
		    }
		    else if (smf_fine.size() == 2)
			//else if (iter > 1) 
		    {

			//BL_ASSERT(smf_fine[0]->boxArray() == smf_fine[1]->boxArray());
			//PArray<MultiFab> raii(PArrayManage);
			//MultiFab * dmf;

			if (m_fabs.boxArray() == smf_fine[0]->boxArray())
			{
			    //dmf = &m_fabs;
			    destGraph = new RegionGraph(m_fabs.IndexArray().size());

			    Perilla::multifabBuildFabCon(destGraph, m_fabs, geom_fine->periodicity());
			    m_amrlevel.parent->graphArray[level].push_back(destGraph);
                            regionList.push_back(destGraph);
			}
			else
			{
			    //dmf = raii.push_back(new MultiFab(smf_fine[0]->boxArray(), m_amrlevel.dmap, ncomp, 0));
			    dmff = new MultiFab(smf_fine[0]->boxArray(), smf_fine[0]->DistributionMap(), ncomp, 0);
			    //dmff->initVal(); // for Perilla NUMA
			    destGraph = new RegionGraph(m_fabs.IndexArray().size());
			    fsrcGraph = new RegionGraph(dmff->IndexArray().size());
			    fsrcGraph->buildTileArray(*dmff);
                              regionList.push_back(destGraph);
                              regionList.push_back(fsrcGraph);
                              mfList.push_back(dmff);

			    Perilla::multifabExtractCopyAssoc( destGraph, fsrcGraph, m_fabs, *dmff, ncomp, m_fabs.nGrow(), 0, geom_fine->periodicity());
			    m_amrlevel.parent->graphArray[level].push_back(destGraph);
			    m_amrlevel.parent->graphArray[level].push_back(fsrcGraph);
			}
		    }
		    else
		    {
			amrex::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
		    }
		    //-------------------- FillFromTwoLevels initialization completed

		} // if(level==1 OR ProperlyNested)
		else
		{
		    amrex::Abort("initFillPatch: level is not properly nested");
		}
	    }
	    DComp += NComp;
	}

	destGraph->buildTileArray(m_fabs);
	destGraph->buildTileArray_gtbx(m_leveldata,boxGrow);
	//MemOpt 
	//m_fabs.clear();
    }

    void AsyncFillPatchIterator::SendIntraLevel (RGIter& rgi,
	    int  boxGrow,
	    Real time,
	    int  index,
	    int  scomp,
	    int  ncomp,
	    int iteration,
	    int f,
	    bool singleT)
    {
	if(rgi.currentItr != rgi.totalItr)
	    return;

	const int level = m_amrlevel.level;



	int ncycle = m_amrlevel.parent->nCycle(level);
	unsigned char pushLevel = 0x02;
	PushOnly(boxGrow, time, index, scomp, ncomp, f, pushLevel, singleT);
    }
    void AsyncFillPatchIterator::SendInterLevel (RGIter* rgi,
	    int  boxGrow,
	    Real time,
	    int  index,
	    int  scomp,
	    int  ncomp,
	    int iteration,
	    int f,
	    bool singleT)
    {
	if(rgi->currentItr != rgi->totalItr)
	    return;

	if(m_amrlevel.level-1 < m_amrlevel.parent->finestLevel())
	{
	    unsigned char tuc = 0x01;
	    PushOnly(boxGrow, time+((iteration-1)*m_amrlevel.parent->dtLevel(m_amrlevel.level)), index, scomp, ncomp, f, tuc, singleT);
	}
    }

    void AsyncFillPatchIterator::SendInterLevel (RGIter& rgi,
            int  boxGrow,
            Real time,
            int  index,
            int  scomp,
            int  ncomp,
            int iteration,
            int f,
            bool singleT)
    {
	SendInterLevel(&rgi, boxGrow, time, index, scomp, ncomp, iteration, f, singleT);
    }


    void AsyncFillPatchIterator::PushOnly (int  boxGrow,
	    Real time,
	    int  index,
	    int  scomp,
	    int  ncomp,
	    int f,
	    unsigned char pushLevel,
	    bool singleT)
    {
	BL_PROFILE("FillPatchIterator::InitializePush");
	BL_ASSERT(scomp >= 0);
	BL_ASSERT(ncomp >= 1);
	BL_ASSERT(0 <= index && index < AmrLevel::desc_lst.size());

	//const IndexType& boxType = m_leveldata.boxArray().ixType();
	const int level = m_amrlevel.level;

	int myProc = amrex::ParallelDescriptor::MyProc();
	for (int i = 0, DComp = 0; i < m_range.size(); i++)
	{
	    if(i>0)
		amrex::Abort("**** Error in FillPatchIterator::Initialize:  non contigeous components not implemented");

	    const int SComp = m_range[i].first;
	    const int NComp = m_range[i].second;

	    if (level == 0)
	    {
		Vector<MultiFab*>                  tsmf;
		Vector<Real>                       tstime;
		StateData& statedata = m_amrlevel.state[index];
		statedata.getData(tsmf,tstime,time);
		FillPatchSingleLevelPush (*(m_amrlevel.parent), m_fabs, time, tsmf, tstime, destGraph, fsrcGraph, f, dmf, SComp, DComp, NComp, *geom, *physbcf, singleT);
	    }else{
		if (level == 1 || isProperlyNested)
		{
		    FillFromTwoLevelsPush(time, index, SComp, DComp, NComp, f, pushLevel, singleT);
		}else {
		    amrex::Abort("**** Error in FillPatchIterator::Initialize:  !ProperlyNested not implemented");
		}
	    }
	    DComp += NComp;
	}
    }

    void AsyncFillPatchIterator::PullOnly (int  boxGrow,
	    Real time,
	    int  index,
	    int  scomp,
	    int  ncomp,
	    int f,
	    bool singleT)
    {
	BL_PROFILE("FillPatchIterator::InitializePull");
	BL_ASSERT(scomp >= 0);
	BL_ASSERT(ncomp >= 1);
	BL_ASSERT(0 <= index && index < AmrLevel::desc_lst.size());

	//const IndexType& boxType = m_leveldata.boxArray().ixType();
	const int level = m_amrlevel.level;

	for (int i = 0, DComp = 0; i < m_range.size(); i++)
	{
	    if(i>0)
		amrex::Abort("**** Error in FillPatchIterator::Initialize:  non contigeous components not implemented");

	    const int SComp = m_range[i].first;
	    const int NComp = m_range[i].second;

	    if (level == 0)
	    {
		FillPatchSingleLevelPull (m_fabs, time, smf, stime, destGraph, fsrcGraph, f, SComp, DComp, NComp, *geom, *physbcf, singleT);
	    }
	    else
	    {
		if (level == 1 || isProperlyNested)
		{
		    FillFromTwoLevelsPull(time, index, SComp, DComp, NComp, f, singleT);
		} else {
		    amrex::Abort("**** Error in FillPatchIterator::Initialize:  !ProperlyNested not implemented");
		}
	    }
	    //if(WorkerThread::isTeamMasterThread(tid))
	    {
		const MultiFab& mf_fillpatched = m_fabs;
		if(singleT)
		{
		    for(int t=0; t<destGraph->fabTiles_gtbx[f]->numTiles; t++)
		    {
			const Box& bx = *(destGraph->fabTiles_gtbx[f]->tileBx[t]);
			MultiFab::Copy(m_leveldata, mf_fillpatched, f, 0, DComp, ncomp, bx);
		    }
		}
		else
		{
		    int nt = perilla::wtid();
		    int totalCompThreads= perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS;
		    for(int t=nt; t<destGraph->fabTiles_gtbx[f]->numTiles; t+= totalCompThreads)
		    {
			const Box& bx = *(destGraph->fabTiles_gtbx[f]->tileBx[t]);
			MultiFab::Copy(m_leveldata, mf_fillpatched, f, 0, DComp, ncomp, bx);
		    }
		}
	    }
	    DComp += NComp;
	}
    }


void
AsyncFillPatchIterator::PullOnly (MultiFab& dest,
                                  int  boxGrow,
                                  Real time,
                                  int  index,
                                  int  scomp,
                                  int  ncomp,
                                  int f,
                                  bool singleT)
{
    BL_PROFILE("FillPatchIterator::InitializePull");
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(0 <= index && index < AmrLevel::desc_lst.size());

    int myProc = amrex::ParallelDescriptor::MyProc();

    //const IndexType& boxType = m_leveldata.boxArray().ixType();
    const int level = m_amrlevel.level;

    for (int i = 0, DComp = 0; i < m_range.size(); i++)
      {
        if(i>0)
          amrex::Abort("**** Error in FillPatchIterator::Initialize:  non contigeous components not implemented");

        const int SComp = m_range[i].first;
        const int NComp = m_range[i].second;

        if (level == 0)
          {

            //double start_time_wtime = omp_get_wtime();            
            try{
              //MemOpt
              FillPatchSingleLevelPull (dest, time, smf, stime, destGraph, fsrcGraph, f, SComp, DComp, NComp, *geom, *physbcf, singleT);
              //amrex::FillPatchSingleLevelPull (m_fabs, time, smf, stime, destGraph, fsrcGraph, f, tid, SComp, DComp, NComp, *geom, *physbcf, singleT);
            }
            catch(std::exception& e){
              std::cout<< "AFPI_Receive_FPSLPull: Proc " <<myProc << " tid " << tid <<" exception: " << e.what() <<std::endl;
            }
          }
        else
          {
            if (level == 1 || isProperlyNested)
              {
                try{
                  //MemOpt
                  FillFromTwoLevelsPull(dest, time, index, SComp, DComp, NComp, f,singleT);
                  //FillFromTwoLevelsPull(time, index, SComp, DComp, NComp, f, tid, singleT);
                }
                catch(std::exception& e){
                  std::cout<< "AFPI_Receive_FPTLPull: Proc " <<myProc << " tid " << tid <<" exception: " << e.what() <<std::endl;
                }

              } else {
              amrex::Abort("**** Error in FillPatchIterator::Initialize:  !ProperlyNested not implemented");
            }
          }
    }
}


    void AsyncFillPatchIterator::FillPatchSingleLevelPush (Amr& amr, MultiFab& mf, Real time,
	    Vector<MultiFab*>& smf, const Vector<Real>& stime,
	    RegionGraph* destGraph, RegionGraph* srcGraph, int f,
	    MultiFab *dmf,
	    int scomp, int dcomp, int ncomp,
	    const Geometry& geom, StateDataPhysBCFunct& physbcf, bool singleT)
    {
	BL_PROFILE("FillPatchSingleLevel");
	BL_ASSERT(scomp+ncomp <= smf[0]->nComp());
	BL_ASSERT(dcomp+ncomp <= mf.nComp());
	BL_ASSERT(smf.size() == stime.size());
	BL_ASSERT(smf.size() != 0);

	int tg = perilla::wid();
	int nt = perilla::wtid();

	if (smf.size() == 1)
	{
	    //mf.copy(smf[0], scomp, dcomp, ncomp, 0, mf.nGrow(), geom.periodicity());        
	    Perilla::multifabCopyPush(destGraph, srcGraph, &mf, smf[0],  f, dcomp, scomp, ncomp, mf.nGrow(), 0, singleT);
	}
	else if (smf.size() == 2)
	{
	    BL_ASSERT(smf[0]->boxArray() == smf[1]->boxArray());
	    //PArray<MultiFab> raii(PArrayManage);
	    //MultiFab * dmf;
	    int destcomp;
	    bool sameba;
	    //if (false && mf.boxArray() == smf[0]->boxArray())
	    if (mf.boxArray() == smf[0]->boxArray())
	    {
		std::cout << "FillPatchUtil SLPush Nt Handled" << std::endl;

		//dmf = &mf;
		destcomp = dcomp;
		sameba = true;

		int fis = smf[0]->IndexArray()[f];
		int fid = mf.IndexArray()[f];

		const Box& bx = mf[fid].box();
		mf[fid].linInterp(smf[0]->get(fis),
			scomp,
			smf[1]->get(fis),
			scomp,
			stime[0],
			stime[1],
			time,
			bx,
			destcomp,
			ncomp);
		Perilla::fillBoundaryPush(destGraph, &mf, f);
	    }
	    else
	    {
		destcomp = 0;
		sameba = false;

		int fis = smf[0]->IndexArray()[f];
		int fid = dmf->IndexArray()[f];

		for(int t=0; t<srcGraph->fabTiles[f]->numTiles; t++)
		    if( singleT || t % (perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS) == nt)
		    {
			const Box& bx = *(srcGraph->fabTiles[f]->tileBx[t]);
			if(bx.ok())
			    (*dmf)[fid].linInterp(smf[0]->get(fis),
				    scomp,
				    smf[1]->get(fis),
				    scomp,
				    stime[0],
				    stime[1],
				    time,
				    bx,
				    destcomp,
				    ncomp);
		    }
		if(!singleT)
		    srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS); // Barrier to synchronize team threads
		int src_ngrow = 0;
		int dst_ngrow = mf.nGrow();
		Perilla::multifabCopyPush( destGraph, srcGraph, &mf, dmf, f, dcomp, 0, ncomp, mf.nGrow(), 0 ,singleT);

	    }
	}
	else
	{
	    amrex::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
	}
    }

    void AsyncFillPatchIterator::FillPatchSingleLevelPull (MultiFab& mf, Real time,
	    Vector<MultiFab*>& smf, const Vector<Real>& stime,
	    RegionGraph* destGraph, RegionGraph* srcGraph, int f,
	    int scomp, int dcomp, int ncomp,
	    const Geometry& geom, StateDataPhysBCFunct& physbcf, bool singleT)
    {
	BL_PROFILE("FillPatchSingleLevel");

	BL_ASSERT(scomp+ncomp <= smf[0]->nComp());
	BL_ASSERT(dcomp+ncomp <= mf.nComp());
	BL_ASSERT(smf.size() == stime.size());
	BL_ASSERT(smf.size() != 0);

	int tg = perilla::wid(); 

	if (smf.size() == 1)
	{
	    //mf.copy(smf[0], scomp, dcomp, ncomp, 0, mf.nGrow(), geom.periodicity());      
	    Perilla::multifabCopyPull( destGraph, srcGraph, &mf, smf[0], f, dcomp, scomp, ncomp, mf.nGrow(), 0, singleT);
	}
	else if (smf.size() == 2)
	{
	    BL_ASSERT(smf[0]->boxArray() == smf[1]->boxArray());
	    Vector<MultiFab*> raii;
	    MultiFab * dmf;
	    int destcomp;
	    bool sameba;
	    //if (false && mf.boxArray() == smf[0]->boxArray()) {
	    if (mf.boxArray() == smf[0]->boxArray()) {
		//dmf = &mf;
		destcomp = dcomp;
		sameba = true;
	    } else {

		//dmf = srcGraph->assocMF;              
		destcomp = 0;
		sameba = false;
	    }
	    if (sameba)
	    {
		// Note that when sameba is true mf's BoxArray is nonoverlapping.
		// So FillBoundary is safe.
		//mf.FillBoundary(dcomp,ncomp,geom.periodicity());

		Perilla::fillBoundaryPull(destGraph, &mf, f, singleT);

		//std::cout << "After sameba fBPull" << std::endl;
	    }
	    else
	    {
		int src_ngrow = 0;
		int dst_ngrow = mf.nGrow();
		MultiFab* dummyMF;

		//mf.copy(*dmf, 0, dcomp, ncomp, src_ngrow, dst_ngrow, geom.periodicity());

		Perilla::multifabCopyPull( destGraph, srcGraph, &mf, dummyMF, f, dcomp, 0, ncomp, mf.nGrow(), 0, singleT);

	    }
	}
	    else {
		amrex::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
	    }

	    if(!singleT)
		destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
	}


	void AsyncFillPatchIterator::FillPatchTwoLevelsPush (Amr& amr, MultiFab& mf, Real time,
		Vector<MultiFab*>& cmf, Vector<Real>& ct,
		Vector<MultiFab*>& fmf, Vector<Real>& ft,
		RegionGraph* destGraph, RegionGraph* csrcGraph, RegionGraph* fsrcGraph, int f,
		AsyncFillPatchIterator* fpIter,
		MultiFab *dmf,
		MultiFab *dmff,
		int scomp, int dcomp, int ncomp,
		const Geometry& cgeom, const Geometry& fgeom,
		StateDataPhysBCFunct& cbc, StateDataPhysBCFunct& fbc,
		const IntVect& ratio,
		Interpolater* mapper, const Vector<BCRec>& bcs, unsigned char pushLevel, bool singleT)
	{   
	    BL_PROFILE("FillPatchTwoLevels");

	    int ngrow = mf.nGrow();

	    if(f>=0){//fill only this fab
		if(pushLevel & 0x01 )
		{   
		    if (ngrow > 0 || mf.getBDKey() != fmf[0]->getBDKey())
		    {

			if (!fpIter->m_fpc->ba_crse_patch.empty())
			{   
			    FillPatchSingleLevelPush(amr, *(fpIter->m_mf_crse_patch), time, cmf, ct, fpIter->m_rg_crse_patch, csrcGraph, f, dmf, scomp, 0, ncomp, cgeom, cbc, singleT);
			}
		    }
		}
		if((pushLevel & 0x02) && (pushLevel != 0x03))
		{    
		    FillPatchSingleLevelPush(amr, mf, time, fmf, ft, destGraph, fsrcGraph, f, dmff, scomp, dcomp, ncomp, fgeom, fbc, singleT);
		}
	    }else{ //fill the whole multifab
		if(pushLevel & 0x01 && pushLevel & 0x02)
		{   
		    int tg = perilla::wid();
		    for(int fi=0; fi < fmf[0]->IndexArray().size(); fi++)
		    {   
			if(WorkerThread::isMyRegion(tg,fi))
			{ 
			    FillPatchSingleLevelPush(amr, mf, time, fmf, ft, destGraph, fsrcGraph, fi, dmff, scomp, dcomp, ncomp, fgeom, fbc, singleT);
			}
		    }
		}
		if(pushLevel & 0x04)
		{
		    int tg = perilla::wid();
		    for(int fi=0; fi < fmf[0]->IndexArray().size(); fi++)
		    {
			if(WorkerThread::isMyRegion(tg,fi))
			{
			    FillPatchSingleLevelPush(amr, mf, time, fmf, ft, destGraph, fsrcGraph, fi, dmff, scomp, dcomp, ncomp, fgeom, fbc, singleT);
			}
		    }
		}
	    }

#if 0
	    BL_PROFILE("FillPatchTwoLevels");

	    int ngrow = mf.nGrow();

	    if(pushLevel & 0x01 )
	    { 
		if (ngrow > 0 || mf.getBDKey() != fmf[0]->getBDKey())
		{ 
		    if ( ! fpIter->m_fpc->ba_crse_patch.empty())
		    {                 
			FillPatchSingleLevelPush(amr, *(fpIter->m_mf_crse_patch), time, cmf, ct, fpIter->m_rg_crse_patch, csrcGraph, f, tid, dmf, scomp, 0, ncomp, cgeom, cbc, singleT);

		    }
		}


		if(tf == 0 && (pushLevel & 0x02) )
		{ 
		    int tg = WorkerThread::groupID(tid);
		    for(int fi=0; fi < fmf[0]->IndexArray().size(); fi++)
		    { 
			if(WorkerThread::isMyRegion(tg,fi))
			{ 
			    FillPatchSingleLevelPush(amr, mf, time, fmf, ft, destGraph, fsrcGraph, fi, tid, dmff, scomp, dcomp, ncomp, fgeom, fbc, singleT);
			}
		    }
		}
	    }
	    if(tf == 0 && (pushLevel & 0x04) )
	    {
		int tg = WorkerThread::groupID(tid);
		for(int fi=0; fi < fmf[0]->IndexArray().size(); fi++)
		{
		    if(WorkerThread::isMyRegion(tg,fi))
		    {
			FillPatchSingleLevelPush(amr, mf, time, fmf, ft, destGraph, fsrcGraph, fi, tid, dmff, scomp, dcomp, ncomp, fgeom, fbc, singleT);
		    } 
		}   
	    }   

	    if((pushLevel & 0x02) && (pushLevel != 0x03))
	    {
		FillPatchSingleLevelPush(amr, mf, time, fmf, ft, destGraph, fsrcGraph, f, tid, dmff, scomp, dcomp, ncomp, fgeom, fbc, singleT);
	    }
#endif
	}     


	void AsyncFillPatchIterator::FillPatchTwoLevelsPull (MultiFab& mf, Real time,
		Vector<MultiFab*>& cmf, Vector<Real>& ct,
		Vector<MultiFab*>& fmf, Vector<Real>& ft,
		RegionGraph* destGraph, RegionGraph* csrcGraph, RegionGraph* fsrcGraph, int f,
		AsyncFillPatchIterator* fpIter,
		int scomp, int dcomp, int ncomp,
		const Geometry& cgeom, const Geometry& fgeom,
		StateDataPhysBCFunct& cbc, StateDataPhysBCFunct& fbc,
		const IntVect& ratio,
		Interpolater* mapper, const Vector<BCRec>& bcs, bool singleT)
	{
	    BL_PROFILE("FillPatchTwoLevels");
	    int ngrow = mf.nGrow();
	    int tg = WorkerThread::perilla_wid();
	    int nt = WorkerThread::perilla_wtid();

	    if (ngrow > 0 || mf.getBDKey() != fmf[0]->getBDKey())
	    {

		if ( ! fpIter->m_fpc->ba_crse_patch.empty())
		{

		    int idummy1=0, idummy2=0;
		    bool cc = fpIter->m_fpc->ba_crse_patch.ixType().cellCentered();
		    {
			int gi = mf.IndexArray()[f];
			for(int i=0; i<destGraph->task[f]->depTaskIDs.size();i++)
			{ 
			    int li = destGraph->task[f]->depTaskIDs[i];
			    int mfi = fpIter->m_mf_crse_patch[0].IndexArray()[li];
			    FillPatchSingleLevelPull(*(fpIter->m_mf_crse_patch), time, cmf, ct, fpIter->m_rg_crse_patch, csrcGraph, li, scomp, 0, ncomp, cgeom, cbc, singleT);
			}
			if(!singleT)
			    destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);

			int nt = WorkerThread::perilla_wtid();
			Box fdomain = fgeom.Domain();
			for(int i=0; i<destGraph->task[f]->depTaskIDs.size();i++)
			{
			    int li = destGraph->task[f]->depTaskIDs[i];
			    int mfi = fpIter->m_mf_crse_patch[0].IndexArray()[li];
			    if(singleT)
			    {
				const Box& dbx = fpIter->m_fpc->dst_boxes[li];
				//Array<BCRec> bcr(ncomp);
				Vector<BCRec> bcr(ncomp);
				amrex::setBC(dbx,fdomain,scomp,0,ncomp,bcs,bcr);

				mapper->interp(fpIter->m_mf_crse_patch[0][mfi],
					0,
					mf[gi],
					dcomp,
					ncomp,
					dbx,
					ratio,
					cgeom,
					fgeom,
					bcr,
					idummy1, idummy2);
			    }
			    else
			    {
				if(!cc)
				{
				    if(WorkerThread::perilla_isMasterWorkerThread())
				    {
					const Box& dbx = fpIter->m_fpc->dst_boxes[li];
					//Box fdomain = fgeom.Domain();

					Vector<BCRec> bcr(ncomp);
					amrex::setBC(dbx,fdomain,scomp,0,ncomp,bcs,bcr);

					mapper->interp(fpIter->m_mf_crse_patch[0][mfi],
						0,
						mf[gi],
						dcomp,
						ncomp,
						dbx,
						ratio,
						cgeom,
						fgeom,
						bcr,
						idummy1, idummy2);

				    }
				}
				else
				{

				    for(int j=0; j < fpIter->m_rg_crse_patch->fabTiles[li]->numTiles; j++)
					if(j % (perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS) == nt)
					{
					    const Box& dbx = *(fpIter->m_rg_crse_patch->fabTiles[li]->tileBx[j]);
					    if(dbx.ok())
					    {
						Vector<BCRec> bcr(ncomp);
						amrex::setBC(dbx,fdomain,scomp,0,ncomp,bcs,bcr);
						mapper->interp(fpIter->m_mf_crse_patch[0][mfi],
							0,
							mf[gi],
							dcomp,
							ncomp,
							dbx,
							ratio,
							cgeom,
							fgeom,
							bcr,
							idummy1, idummy2);
					    }
					}

				}
			    }
			    //destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
			}
			if(!singleT)
			    destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
		    }
		}
	    }

	    FillPatchSingleLevelPull(mf, time, fmf, ft, destGraph, fsrcGraph, f, scomp, dcomp, ncomp, fgeom, fbc, singleT);
	}

#if 0
	void AsyncFillPatchIterator::FillPatchTwoLevelsPull (MultiFab& mf, Real time,
		Vector<MultiFab*>& cmf, const Vector<Real>& ct,
		Vector<MultiFab*>& fmf, const Vector<Real>& ft,
		RegionGraph* destGraph, RegionGraph* csrcGraph, RegionGraph* fsrcGraph, int f, int tid,
		AsyncFillPatchIterator* fpIter,
		int scomp, int dcomp, int ncomp,
		const Geometry& cgeom, const Geometry& fgeom,
		PhysBCFunctBase& cbc, PhysBCFunctBase& fbc,
		const IntVect& ratio,
		Interpolater* mapper, const Vector<BCRec>& bcs, bool singleT)
	{
	    BL_PROFILE("FillPatchTwoLevels");

	    int ngrow = mf.nGrow();

	    int tg = perilla:wid();//WorkerThread::groupID(tid);
	    int nt = perilla::wtid();//WorkerThread::numaTID(tid);


	    int myProc = ParallelDescriptor::MyProc();
	    //std::ofstream fout;
	    //fout.open(std::to_string(myProc)+ "_" + std::to_string(tid) + ".txt", std::fstream::app);


	    if (ngrow > 0 || mf.getBDKey() != fmf[0]->getBDKey())
	    {

		//fout << "FPTL fpIter->m_fpc->ba_crse_patch.empty() " << fpIter->m_fpc->ba_crse_patch.empty() << std::endl;

		if ( ! fpIter->m_fpc->ba_crse_patch.empty())
		{

		    int idummy1=0, idummy2=0;
		    bool cc = fpIter->m_fpc->ba_crse_patch.ixType().cellCentered();

		    //std::cout << "Check CC : " << cc << std::endl;
		    //#ifdef _OPENMP
		    //#pragma omp parallel if (cc)
		    //#endif


		    //for (MFIter mfi(*(fpIter->m_mf_crse_patch),false,false); mfi.isValid(); ++mfi)
		    {
			//int li = mfi.LocalIndex();
			//int gi = fpIter->m_fpc->dst_idxs[li];
			//if(gi == mf.IndexArray()[f])
			int gi = mf.IndexArray()[f];
			//if(gi == f)

			// fout << "FPTL gi " << gi << " f " << f << " destGraph->task[f]->depTaskIDs.size() " << destGraph->task[f]->depTaskIDs.size() << std::endl;

			//double start_time_wtime = omp_get_wtime();

			for(int i=0; i<destGraph->task[f]->depTaskIDs.size();i++)
			{
			    int li = destGraph->task[f]->depTaskIDs[i];
			    int mfi = fpIter->m_mf_crse_patch[0].IndexArray()[li];

			    //fout << "Calling FPSL for dependent gi "<< gi << " li " << li << std::endl;
			    FillPatchSingleLevelPull(*(fpIter->m_mf_crse_patch), time, cmf, ct, fpIter->m_rg_crse_patch, csrcGraph, li, tid, scomp, 0, ncomp, cgeom, cbc, singleT);
			}
			//if(!singleT)
			//destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);

			//double end_time_wtime = omp_get_wtime();
			/*if(singleT)
			  Perilla::getPPPTimeSplit[8] += end_time_wtime - start_time_wtime;
			  else
			  if(WorkerThread::isTeamMasterThread(tid))
			  Perilla::getPPPTimeSplit[8] += end_time_wtime - start_time_wtime;
			 */
			//if(myProc == 0 && nt == 4)
			//std::cout << "CC " << cc << " DepTasks " << destGraph->task[f]->depTaskIDs.size() << std::endl;


			//start_time_wtime = omp_get_wtime();
			int nt = perilla::wtid();//WorkerThread::numaTID(tid);
			Box fdomain = fgeom.Domain();
			for(int i=0; i<destGraph->task[f]->depTaskIDs.size();i++)
			{
			    int li = destGraph->task[f]->depTaskIDs[i];
			    int mfi = fpIter->m_mf_crse_patch[0].IndexArray()[li];
			    /*              
					    const Box& dbx1 = fpIter->m_fpc->dst_boxes[li];

					    std::ofstream fout;
					    fout.open(std::to_string(myProc)+ "_" + std::to_string(tid) + ".txt", std::fstream::app);
					    fout << "i "<< i << " depsiz " << destGraph->task[f]->depTaskIDs.size() <<" li "<< li << " f "<< f << " mfi "<< mfi << " gi " << gi << std::endl; 
					    fout <<" numFabs " << fpIter->m_rg_crse_patch->fabTiles.size() << " numTls " << fpIter->m_rg_crse_patch->fabTiles[li]->numTiles << " ndbs " << fpIter->m_fpc->dst_boxes.size() << std::endl;
			    //fout <<"dbx " << dbx << std::endl;
			    fout <<"dst_bxs " << dbx1 << std::endl;
			    fout <<"fine bx " << mf[gi].box() << std::endl;
			    fout <<"crse bx " << fpIter->m_mf_crse_patch[0][mfi].box() << std::endl;
			    fout.close();
			     */

			    if(singleT)
			    {
				const Box& dbx = fpIter->m_fpc->dst_boxes[li];
				//Box fdomain = fgeom.Domain();

				Vector<BCRec> bcr(ncomp);
				amrex::setBC(dbx,fdomain,scomp,0,ncomp,bcs,bcr);

				mapper->interp(fpIter->m_mf_crse_patch[0][mfi],
					0,
					mf[gi],
					dcomp,
					ncomp,
					dbx,
					ratio,
					cgeom,
					fgeom,
					bcr,
					idummy1, idummy2);
			    }
			    else
			    {
				if(!cc)
				{
				    if(WorkerThread::isTeamMasterThread(tid))
				    {
					const Box& dbx = fpIter->m_fpc->dst_boxes[li];
					//Box fdomain = fgeom.Domain();

					Array<BCRec> bcr(ncomp);
					amrex::setBC(dbx,fdomain,scomp,0,ncomp,bcs,bcr);

					mapper->interp(fpIter->m_mf_crse_patch[0][mfi],
						0,
						mf[gi],
						dcomp,
						ncomp,
						dbx,
						ratio,
						cgeom,
						fgeom,
						bcr,
						idummy1, idummy2);

				    }
				}
				else
				{
				    //std::cout << "myP " << myProc << " nt "<< nt << " li "<< li << " mfi " << mfi << " ntiles " << fpIter->m_rg_crse_patch->fabTiles.size()<< std::endl;
				    ///for(int j=0; j < fpIter->m_rg_crse_patch->fabTiles[f]->numTiles; j++)
				    for(int j=0; j < fpIter->m_rg_crse_patch->fabTiles[li]->numTiles; j++)
					if(j % (perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS) == nt-perilla::NUM_COMM_THREADS)
					    ///if(i % (perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS) == nt-perilla::NUM_COMM_THREADS)
					{

					    //if(myProc == 0 && nt == 4)
					    //std::cout << "CC " << cc << " DepTasks " << destGraph->task[f]->depTaskIDs.size() << " i " << i  << std::endl;
					    ///const Box& dbx = fpIter->m_fpc->dst_boxes[li];
					    ///const Box& dbx = *(fpIter->m_rg_crse_patch->fabTiles[f]->tileBx[j]);
					    const Box& dbx = *(fpIter->m_rg_crse_patch->fabTiles[li]->tileBx[j]);
					    //Box fdomain = fgeom.Domain();

					    //if(myProc == 0 && nt == 4)
					    //std::cout << "CC " << cc << " DepTasks " << destGraph->task[f]->depTaskIDs.size() << " i " << i  << " dbx " << dbx << std::endl;


					    //fout << "FPTL interping gi "<< gi << " li " << li<< " dbx " << dbx << std::endl;
					    if(dbx.ok())
					    {
						Vector<BCRec> bcr(ncomp);
						amrex::setBC(dbx,fdomain,scomp,0,ncomp,bcs,bcr);
						mapper->interp(fpIter->m_mf_crse_patch[0][mfi],
							0,
							mf[gi],
							dcomp,
							ncomp,
							dbx,
							ratio,
							cgeom,
							fgeom,
							bcr,
							idummy1, idummy2);
					    }
					}
				}
			    }
			    //destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
			}
			if(!singleT)
			    destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-perilla::NUM_COMM_THREADS);
			/*
			   end_time_wtime = omp_get_wtime();
			   if(singleT)
			   Perilla::getPPPTimeSplit[9] += end_time_wtime - start_time_wtime;
			   else
			   if(WorkerThread::isTeamMasterThread(tid))
			   Perilla::getPPPTimeSplit[9] += end_time_wtime - start_time_wtime;
			 */
		    }
		}
	    }
	    /*
	       int mfi = mf.IndexArray()[f];
	       const Box& bx = mf[mfi].box();
	    //if(mfi == 0 )
	    {
	    fout << "Before second FPSL at FPSL mfi " << mfi << " bx " << bx.smallEnd() << bx.bigEnd() <<  std::endl;
	    for(int i=bx.smallEnd(0); i<=bx.smallEnd(0); i++)
	    {
	    for(int j=bx.smallEnd(1); j<=bx.bigEnd(1); j++)
	    {
	    for(int k=bx.smallEnd(2); k<=bx.bigEnd(2); k++)
	    {
	    fout << mf[mfi](IntVect(i,j,k))  << " ";
	    }
	    fout << std::endl;
	    }
	    fout << std::endl;
	    }
	    }
	     */
	    //fout.close();
	    //double start_time_wtime = omp_get_wtime();

	    FillPatchSingleLevelPull(mf, time, fmf, ft, destGraph, fsrcGraph, f, tid, scomp, dcomp, ncomp, fgeom, fbc, singleT);
	    /*      
		    double end_time_wtime = omp_get_wtime();
		    if(singleT)
		    Perilla::getPPPTimeSplit[10] += end_time_wtime - start_time_wtime;
		    else
		    if(WorkerThread::isTeamMasterThread(tid))
		    Perilla::getPPPTimeSplit[10] += end_time_wtime - start_time_wtime;
	     */
	}


#endif



	void AsyncFillPatchIterator::initialSend(Vector<amrex::AsyncFillPatchIterator*> afpi,
		Vector<amrex::AsyncFillPatchIterator*> upper_afpi,
		int  boxGrow,
		Real time,
		int  state_indx,
		int  scomp,
		int  ncomp,
		int  iteration)
	{
	    int myProc = amrex::ParallelDescriptor::MyProc();
	    int level = afpi[iteration-1]->m_amrlevel.level;
	    if(level == 0 && iteration == 1)
	    {
		int tg = perilla::wid();
		for(int f=0; f < afpi[iteration-1]->m_fabs.IndexArray().size(); f++)
		{
		    if(WorkerThread::isMyRegion(tg, f))
		    {
			for(int i=0; i < afpi[iteration-1]->m_amrlevel.parent->nCycle(level); i++){
			    //fill neighbor fabs of the same AMR level
			    afpi[i]->PushOnly( boxGrow, time+(i*afpi[iteration-1]->m_amrlevel.parent->dtLevel(level)), state_indx, scomp, ncomp, f, 0xFF, false);
			}
		    }
		}
	    }

	    if(level < afpi[iteration-1]->m_amrlevel.parent->finestLevel())
	    {
		int i = 0;
		unsigned char tuc = 0x04;
		//init Fill Patch at the next finer AMR level
		upper_afpi[i]->PushOnly(boxGrow, time+(i*afpi[iteration-1]->m_amrlevel.parent->dtLevel(level+1)), state_indx, scomp, ncomp, -1/* all FABs*/, tuc, false);
	    }
	}

	void AsyncFillPatchIterator::Receive (RGIter& rgi,
		    int  boxGrow,
		    Real time,
		    int  index,
		    int  scomp,
		    int  ncomp,
		    int f,
		    bool singleT)
	    {
		if(rgi.currentItr != 1)
		    return;

		PullOnly(boxGrow, time, index, scomp, ncomp, f, singleT);
	    }

	void AsyncFillPatchIterator::Receive (RGIter* rgi,
		    int  boxGrow,
		    Real time,
		    int  index,
		    int  scomp,
		    int  ncomp,
		    int f,
		    bool singleT)
	    {
		if(rgi->currentItr != 1)
		    return;

		PullOnly(boxGrow, time, index, scomp, ncomp, f, singleT);
	    }


	void AsyncFillPatchIterator::Receive (RGIter& rgi,
		    MultiFab& dest,
		    int  boxGrow,
		    Real time,
		    int  index,
		    int  scomp,
		    int  ncomp,
		    int f,
		    bool singleT)
	    {
		if(rgi.currentItr != 1)
		    return;

		PullOnly(dest, boxGrow, time, index, scomp, ncomp, f, singleT);
	    }


	void AsyncFillPatchIterator::Receive (RGIter* rgi,
		    MultiFab& dest,
		    int  boxGrow,
		    Real time,
		    int  index,
		    int  scomp,
		    int  ncomp,
		    int f,
		    bool singleT)
	    {
		if(rgi->currentItr != 1)
		    return;

		PullOnly(dest, boxGrow, time, index, scomp, ncomp, f, singleT);
	    }
    }//end amrex namespace

#endif
    //end USE_PERILLA
