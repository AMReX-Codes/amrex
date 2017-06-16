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
#include <AMReX_BLFort.H>
#include "AMReX_EBCellFactory.H"
#include "AMReX_VisMF.H"
#include "AMReX_PlotFileUtil.H"
#include "AMReX_BaseUmap.H"

#include "AMReX_MeshRefine.H"
#include "AMReX_EBStruct.H"

using namespace amrex;

typedef CNode::CutCell CutCell;

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

    BoxArray ba(paramsCrse.coarsestDomain);
    int max_grid_size = 16;
    pp.query("max_grid_size",max_grid_size);
    ba.maxSize(max_grid_size);
    DistributionMapping dm(ba);
    std::vector<EBLevelGrid> eblg_tmp;
    getAllIrregEBLG(eblg_tmp, paramsFine);
    eblg_fine = eblg_tmp[0];
    getAllIrregEBLG(eblg_tmp, paramsCrse);
    eblg_crse = eblg_tmp[0];
}

static void
BuildFortranGraph(FabArray<BaseUmap<CutCell> >& cmap,
                  const EBLevelGrid& eblg)
{
    for (MFIter mfi(cmap); mfi.isValid(); ++mfi)
    {
        EBISBox ebis_box = eblg.getEBISL()[mfi];
        const Box gbox = cmap[mfi].box() & ebis_box.getDomain();

	int ebCellID = 0; // Reset to zero for each box
        if (!ebis_box.isAllRegular())
        {
            const EBGraph& ebgr = ebis_box.getEBGraph();
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
                            cc.ebCellID = ebCellID++;
                            const VolIndex& vof = gbox_vofs[icc];
                            for (int idir = 0; idir < SpaceDim; idir++)
                            {
                                for (SideIterator sit; sit.ok(); ++sit)
                                {
                                    int iside = sit() == Side::Lo ? 0 : 1;
                                    std::vector<FaceIndex> faces = ebgr.getFaces(vof,idir,sit());

				    IntVect iv_nbr = iv + sign(sit())*BASISV(idir);
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
                                            }
                                            nValid++;
					    BL_ASSERT(nValid<NCELLMAX);
                                        }
                                    }
                                }
                            }
                            cmap[mfi].setVal(cc,iv,0,icc);
                        }
                    }
                }
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
    std::cout << "CNode is POD: " << std::is_pod<CNode>::value << std::endl;
    std::cout << "FNode is POD: " << std::is_pod<FNode>::value << std::endl;
    std::cout << "Size of CNode: " << sizeof(CNode)/sizeof(int) << " ints" << std::endl;
    std::cout << "Size of FNode: " << sizeof(FNode)/sizeof(int) << " ints" << std::endl;

    EBLevelGrid eblg_fine, eblg_crse;
    int ratio;
    get_EBLGs(eblg_fine,eblg_crse,ratio);
    const BoxArray& ba_crse = eblg_crse.getBoxArray();
    const DistributionMapping& dm_crse = eblg_crse.getDM();

    int nComp =1;
    int nGrow =1;
    BaseUmapFactory<CutCell> umap_factory(5);
    FabArray<BaseUmap<CutCell> > cmap_crse(ba_crse,dm_crse,nComp, nGrow,MFInfo(), umap_factory);

    BuildFortranGraph(cmap_crse,eblg_crse);


    MultiFab vfrac(ba_crse, dm_crse,  nComp, 0);

    for(MFIter  mfi(vfrac); mfi.isValid(); ++mfi)
    {
        const EBISBox& ebis_box = eblg_crse.getEBISL()[mfi];
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
    Geometry geom_crse(eblg_crse.getDomain(),&rb);
    WriteSingleLevelPlotfile("pltfile",vfrac,name,geom_crse,0,0,"CartGrid-V2.0");

    // Level transfers
    // for (MFIter mfi(eblg_fine.getBoxArray(),eblg_fine.getDM()); mfi.isValid(); ++mfi)
    // {
    //     const Box& fbox = mfi.validbox();
    //     const Box cbox = amrex::coarsen(fbox,ratio);
    //     int cnum_crse = cnodes_crse[mfi.index()].size();
    //     average(fbox.loVect(), fbox.hiVect(),
    //             cnodes_crse[mfi.index()].data(), &cnum_crse,
    //             cnodes_fine[mfi.index()].data(), &cnum_fine,
    //             BL_TO_FORTRAN_N(ebmask_crse[mfi],0),
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
            amrex::Print() << "sphere test passed" << '\n';
        }

        BL_PROFILE_VAR_STOP(pmain);
    }
    amrex::Finalize();
    return 0;
}
