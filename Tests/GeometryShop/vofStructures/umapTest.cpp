#include <cmath>
#include <type_traits>

#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_GeometryShop.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_ParmParse.H"
#include "AMReX_RealVect.H"
#include "AMReX_SphereIF.H"
#include "AMReX_FlatPlateGeom.H"
#include "AMReX_EBGraph.H"
#include "AMReX_EBISBox.H"
#include "AMReX_MultiFab.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_BLProfiler.H"
#include "AMReX_BLFort.H"
#include "AMReX_VisMF.H"
#include "AMReX_PlotFileUtil.H"
#include "AMReX_BaseUmap.H"

#include "AMReX_MeshRefine.H"
#include "AMReX_EBStruct.H"

using namespace amrex;

static void
get_EBLGs(EBLevelGrid& eblg_fine,
          EBLevelGrid& eblg_crse,
          int&         ratio)
{
    GridParameters paramsCrse, paramsFine;
    getGridParameters(paramsCrse, 0, true);

    paramsFine = paramsCrse;
    ratio = 2;
    paramsFine.refine(ratio);

    ParmParse pp;

    EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
    int which_geom = paramsCrse.whichGeom;
    bool verbosity=true;
    pp.query("verbosity",verbosity);
    if (which_geom == 0)
    {
        RealVect origin = RealVect::Zero;
        Real radius = 0.5;
        pp.get("sphere_radius", radius);
        std::vector<Real> centervec(SpaceDim);
        pp.getarr("sphere_center", centervec, 0, SpaceDim);
        RealVect center;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
            center[idir] = centervec[idir];
        }

        //bool insideRegular = false;
        bool insideRegular = true;
        SphereIF sphere(radius, center, insideRegular);

        GeometryShop gshop(sphere, verbosity);
        ebisPtr->define(paramsFine.coarsestDomain, origin, paramsFine.coarsestDx, gshop);
    }
    else if (which_geom == 1) 
    {
        RealVect origin = RealVect::Zero;
        std::vector<Real>  platelovec(SpaceDim);
        std::vector<Real>  platehivec(SpaceDim);

        int normalDir;
        Real plateLoc;
    
        pp.getarr("plate_lo", platelovec, 0, SpaceDim);
        pp.getarr("plate_hi", platehivec, 0, SpaceDim);
        pp.get("plate_location", plateLoc);
        pp.get("plate_normal", normalDir);

        RealVect plateLo, plateHi;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
            plateLo[idir] = platelovec[idir];
            plateHi[idir] = platehivec[idir];
        }
        FlatPlateGeom flat_plate(normalDir, plateLoc, plateLo, plateHi);
        ebisPtr->define(paramsFine.coarsestDomain, origin, paramsFine.coarsestDx, flat_plate);
    }
    else
    {
        amrex::Abort("Unknown geom_type");
    }

    std::vector<EBLevelGrid> eblg_tmpf;
    getAllIrregEBLG(eblg_tmpf, paramsFine);
    eblg_fine = eblg_tmpf[0];

    std::vector<EBLevelGrid> eblg_tmpc;
    getAllIrregEBLG(eblg_tmpc, paramsCrse);
    eblg_crse = eblg_tmpc[0];

}

struct tface
{
    /*
      A simple class to simplify uniquifying faces
     */
    tface() : mL(-1), mR(-1) {}
    tface(int L, int R, const SideIterator& s) : mL(L), mR(R) {if (s()==Side::Lo) flip();}
    bool operator< (const tface& rhs) const {
	if (mL == rhs.mL)
	{
	    return mR < rhs.mR;
	}
	return mL < rhs.mL;
    }
    void flip ()
	{
	    int t=mL;
	    mL = mR;
	    mR = t;
	}
    int mL, mR;
};

static void
BuildFortranGraph(FabArray<BaseUmap<CutCell> >& cmap,
                  std::array<FabArray<BaseUmap<CutFace> >, BL_SPACEDIM>& fmap,
                  const EBLevelGrid& eblg)
{
    for (MFIter mfi(cmap); mfi.isValid(); ++mfi)
    {
        EBISBox ebis_box = eblg.getEBISL()[mfi];
        const Box gbox = cmap[mfi].box() & ebis_box.getDomain();

	//int ebCellID = 0; // Reset to zero for each box
        if (!ebis_box.isAllRegular())
        {
            BaseUmap<CutCell>& cmap_fab = cmap[mfi];
            const EBGraph& ebgr = ebis_box.getEBGraph();
	    std::array<std::map<IntVect, std::set<tface> >, BL_SPACEDIM> tf;

            for (BoxIterator bit(gbox); bit.ok(); ++bit)
            {
                const IntVect& iv = bit();
                if (ebis_box.isIrregular(iv))
                {
		    // Set up cnode data structure
                    std::vector<VolIndex> gbox_vofs = ebis_box.getVoFs(iv);
                    if (gbox_vofs.size() > 0)
                    {
                        int nCells = gbox_vofs.size();
                        for (int icc=0; icc<nCells; ++icc)
                        {
                            CutCell cc;
                            //cc.ebCellID = ebCellID++;

                            const VolIndex& vof = gbox_vofs[icc];
                            for (int idir = 0; idir < SpaceDim; idir++)
                            {
                                for (SideIterator sit; sit.ok(); ++sit)
                                {
                                    int iside = sit() == Side::Lo ? 0 : 1;
                                    std::vector<FaceIndex> faces = ebgr.getFaces(vof,idir,sit());

				    IntVect iv_nbr = iv + sign(sit())*BASISV(idir);
				    IntVect iv_face = iv + iside*BASISV(idir);

                                    int nValid = 0;
                                    for (int iface=0; iface<faces.size(); ++iface)
                                    {
                                        if (!ebis_box.isCovered(iv_nbr))
                                        {
                                            if (ebis_box.isRegular(iv_nbr))
                                            {
                                                cc.nbr[idir][iside][nValid] = REGULAR_CELL;
                                            }
                                            else 
                                            {
                                                cc.nbr[idir][iside][nValid] = faces[iface].cellIndex(sit());
                                                tf[idir][iv_face].insert(tface(icc,cc.nbr[idir][iside][nValid],sit));
                                            }
                                            nValid++;
					    BL_ASSERT(nValid<NCELLMAX);
                                        }
                                    }
                                    cc.Nnbr[idir][iside] = nValid;
                                }
                            }
                            cmap_fab.setVal(cc,iv,0,icc);
                        }
                    }
                }
            }

            const Box& vbox = mfi.validbox();
            for (int idir = 0; idir < SpaceDim; idir++)
            {
                //int ebFaceID = 0;
                const Box fbox = amrex::surroundingNodes(vbox,idir);
                for (BoxIterator bit(fbox); bit.ok(); ++bit)
                {
                    const IntVect& iv_face = bit();
                    std::map<IntVect, std::set<tface> >::const_iterator it=tf[idir].find(iv_face);
                    if (it!=tf[idir].end()) {
                        CutFace cf;
                        const std::set<tface>& tfiv=it->second;
                        int cnt = 0;

                        for (std::set<tface>::iterator fit=tfiv.begin(); fit!=tfiv.end(); ++fit)
                        {
                            cf.cellLo = fit->mL;
                            cf.cellHi = fit->mR;
                            bool is_face_between_reg_and_cut_cells
                                = cf.cellLo==REGULAR_CELL
                                || cf.cellHi==REGULAR_CELL;
                                
                            if ( !is_face_between_reg_and_cut_cells )
                            {
                                /*
                                  Number this face, add it to the umap of faces, and set cmap face ids consistently
                                 */
                                //cf.ebFaceID = ebFaceID++;
                                fmap[idir][mfi].setVal(cf,iv_face,0,cnt);

                                bool found_lo=false;
                                IntVect iv_lo = iv_face - BASISV(idir);
                                CutCell& cc_lo = cmap_fab.getVal(iv_lo,0,cf.cellLo);
                                for(int inbr=0; inbr<cc_lo.Nnbr[idir][1]; ++inbr) {
                                    if (cc_lo.nbr[idir][1][inbr] == cf.cellHi) {
                                        cc_lo.faceID[idir][1][inbr] = cnt;
                                        found_lo = true;
                                    }
                                }
                                if (!found_lo) amrex::Abort();

                                bool found_hi = false;
                                IntVect iv_hi = iv_face;
                                CutCell& cc_hi = cmap_fab.getVal(iv_hi,0,cf.cellHi);
                                for(int inbr=0; inbr<cc_hi.Nnbr[idir][0]; ++inbr) {
                                    if (cc_hi.nbr[idir][0][inbr] == cf.cellLo) {
                                        cc_hi.faceID[idir][0][inbr] = cnt;
                                        found_hi = true;
                                    }
                                }
                                if (!found_hi) amrex::Abort();
                                        
                                cnt++;
                            }
                            else
                            {
                                //cf.ebFaceID = REGULAR_FACE;
                            }
                        }
                    }
                }
            }

        }
    }

    // Dump graphs
    std::cout << "cmap: " << std::endl;
    for (MFIter mfi(cmap); mfi.isValid(); ++mfi)
    {
        const BaseUmap<CutCell>& cmap_fab = cmap[mfi];
        for (BaseUmap<CutCell>::const_iterator umi = cmap_fab.begin(); umi<cmap_fab.end(); ++umi )
        {
            const CutCell& cc = *umi;
            //std::cout << "ebCellID: " << cc.ebCellID;
            const BaseUmap<CutCell>::Tuple& tuple = umi.tuple();
            std::cout << " " << tuple.pos << " L = " << tuple.l << " ncomp= " << tuple.ncomp;
            for (int idir = 0; idir < SpaceDim; idir++)
            {
                for (SideIterator sit; sit.ok(); ++sit)
                {
                    int iside = sit() == Side::Lo ? 0 : 1;
                    std::cout << " (" << idir << ", " << iside << "): " << cc.Nnbr[idir][iside];
                }
            } 
            std::cout << std::endl;
        }
    }

    for (int idir = 0; idir < SpaceDim; idir++)
    {
        std::cout << "fmap[" << idir << "]: " << std::endl;
        for (MFIter mfi(cmap); mfi.isValid(); ++mfi)
        {
            const BaseUmap<CutFace>& fmap_fab = fmap[idir][mfi];
            for (BaseUmap<CutFace>::const_iterator umi = fmap_fab.begin(); umi<fmap_fab.end(); ++umi )
            {
                const CutFace& cf = *umi;
                //std::cout << "ebFaceID: " << cf.ebFaceID;
                const BaseUmap<CutFace>::Tuple& tuple = umi.tuple();
                std::cout << " " << tuple.pos << " L = " << tuple.l << " ncomp= " << tuple.ncomp;
                std::cout << " (" << idir << "): " << cf.cellLo << ", " << cf.cellHi << std::endl;
            }
        }
    }
}

template <class T>
class BaseUmapFactory
    : public FabFactory<BaseUmap<T> >
{
public:
    BaseUmapFactory(int lmax) : m_lmax(lmax) {}

    virtual BaseUmap<T>* create (const Box& box, int nComp, const FabInfo& info, int box_index) const
        {
            return new BaseUmap<T>(box, nComp, m_lmax);
        }

protected:
    int m_lmax;
};

int myTest()
{
    std::cout << "CutCell is POD: " << std::is_pod<CutCell>::value << std::endl;
    std::cout << "CutFace is POD: " << std::is_pod<CutFace>::value << std::endl;
    std::cout << "Size of CutCell: " << sizeof(CutCell)/sizeof(int) << " ints" << std::endl;
    std::cout << "Size of CutFace: " << sizeof(CutFace)/sizeof(int) << " ints" << std::endl;

    EBLevelGrid eblg_fine, eblg_crse;
    int ratio;
    get_EBLGs(eblg_fine,eblg_crse,ratio);
    const BoxArray& ba_fine = eblg_fine.getBoxArray();
    const DistributionMapping& dm_fine = eblg_fine.getDM();

    int nComp =1;
    int nGrow =1;
    BaseUmapFactory<CutCell> cmap_factory(NCELLMAX);
    FabArray<BaseUmap<CutCell> > cmap_fine(ba_fine, dm_fine, nComp, nGrow, MFInfo(), cmap_factory);

    BaseUmapFactory<CutFace> fmap_factory(NFACEMAX);
    std::array<BoxArray,BL_SPACEDIM> fba_fine;
    std::array<FabArray<BaseUmap<CutFace> >, BL_SPACEDIM> fmap_fine;
    for (int idir=0; idir<BL_SPACEDIM; ++idir)
    {
        fba_fine[idir] = ba_fine;
        fba_fine[idir].surroundingNodes(idir);
        fmap_fine[idir].define(fba_fine[idir], dm_fine, nComp, nGrow, MFInfo(), fmap_factory);
    }

    BuildFortranGraph(cmap_fine,fmap_fine,eblg_fine);

#if 1
    MultiFab vfrac(ba_fine, dm_fine,  nComp, 0);

    for(MFIter  mfi(vfrac); mfi.isValid(); ++mfi)
    {
        const EBISBox& ebis_box = eblg_fine.getEBISL()[mfi];
        const Box& vbox = mfi.validbox();
        FArrayBox& vfrac_fab = vfrac[mfi];
        vfrac[mfi].setVal(1);
        if ( !(ebis_box.isAllRegular()) )
        {
            if (ebis_box.isAllCovered()) {
                vfrac_fab.setVal(0);
            }
            else
            {
                for (BoxIterator bit(vbox); bit.ok(); ++bit)
                {
                    const std::vector<VolIndex>& vofs = ebis_box.getVoFs(bit());
                    Real vfrac_tot = 0;
                    for (int ivof=0; ivof<vofs.size(); ++ivof)
                    {
                        vfrac_tot += ebis_box.volFrac(vofs[ivof]);
                    }
                    vfrac_fab(bit(),0) = vfrac_tot;
                }
            }
        }
    }

    Array<std::string> name(1,"vfrac");
    RealBox rb(AMREX_D_DECL(0,0,0),
               AMREX_D_DECL(1,1,1));
    Geometry geom_fine(eblg_fine.getDomain(),&rb);
    WriteSingleLevelPlotfile("pltfile",vfrac,name,geom_fine,0,0,"CartGrid-V2.0");
#endif

    // Level transfers
    // for (MFIter mfi(eblg_fine.getBoxArray(),eblg_fine.getDM()); mfi.isValid(); ++mfi)
    // {
    //     const Box& fbox = mfi.validbox();
    //     const Box cbox = amrex::coarsen(fbox,ratio);
    //     int cnum_fine = cnodes_fine[mfi.index()].size();
    //     average(fbox.loVect(), fbox.hiVect(),
    //             cnodes_fine[mfi.index()].data(), &cnum_fine,
    //             cnodes_fine[mfi.index()].data(), &cnum_fine,
    //             BL_TO_FORTRAN_N(ebmask_fine[mfi],0),
    //             BL_TO_FORTRAN_N(ebmask_fine[mfi],0));
    // }

    return 0;
}

int
main(int argc,char **argv)
{
    amrex::Initialize(argc,argv);
    {
        BL_PROFILE_VAR("main()", pmain);

        int retFlag = myTest();

        if (retFlag != 0)
        {
            amrex::Print() << "non zero return detected = " << retFlag << '\n';
            amrex::Print() << "sphere test failed" << '\n';
        }
        else
        {
            amrex::Print() << "test passed" << '\n';
        }

        BL_PROFILE_VAR_STOP(pmain);
    }
    amrex::Finalize();
    return 0;
}
