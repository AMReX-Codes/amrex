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
#include "umapTest.H"

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
    tface() : mL(-1), mR(-1), mParent(-1) {}
    tface(int L, int R, const SideIterator& s, int parent) : mL(L), mR(R), mParent(parent) {if (s()==Side::Lo) flip();}
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
    int mL, mR, mParent;
};

static void
BuildFortranGraph(FabArray<BaseUmap<CutCell> >& cmap,
                  std::array<FabArray<BaseUmap<CutFace> >, BL_SPACEDIM>& fmap,
                  FabArray<BaseUmap<VolIndex> >& cidx,
                  std::array<FabArray<BaseUmap<FaceIndex> >, BL_SPACEDIM>& fidx,
                  const EBLevelGrid& eblg_fine, EBLevelGrid& eblg_crsn, int ratio)
{
    coarsen(eblg_crsn, eblg_fine, ratio);
    eblg_crsn.setMaxRefinementRatio(ratio);
    int ng_crsn = 1;
    BL_ASSERT(eblg_crsn.getGhost() >= ng_crsn);

    for (MFIter mfi(cmap); mfi.isValid(); ++mfi)
    {
        BaseUmap<CutCell>& cmap_fab = cmap[mfi];
        const EBISBox& ebis_box_crsn = eblg_crsn.getEBISL()[mfi];
        const EBISBox& ebis_box_fine = eblg_fine.getEBISL()[mfi];
        const Box cbox = amrex::grow(eblg_crsn.getDBL()[mfi],ng_crsn) & eblg_crsn.getDomain();
        if (!ebis_box_crsn.isAllRegular())
        {
            const EBGraph& ebgr_crsn = ebis_box_crsn.getEBGraph();
            const EBGraph& ebgr_fine = ebis_box_fine.getEBGraph();
	    std::array<std::map<IntVect, std::set<tface> >, BL_SPACEDIM> tf;
            for (BoxIterator bit(cbox); bit.ok(); ++bit)
            {
                const IntVect& iv_crsn = bit();
                if (ebis_box_crsn.isIrregular(iv_crsn))
                {
		    // Set up cnode data structure
                    std::vector<VolIndex> vofs_crsn = ebis_box_crsn.getVoFs(iv_crsn);
                    int nCells_crsn = vofs_crsn.size();
                    for (int icc_crsn=0; icc_crsn<nCells_crsn; ++icc_crsn)
                    {
                        const VolIndex& vof_crsn = vofs_crsn[icc_crsn];
                        int idx_crsn = vof_crsn.cellIndex();

                        std::vector<VolIndex> vofs_fine = ebis_box_crsn.refine(vof_crsn);
                        int nCells_fine = vofs_fine.size();
                        
                        for (int icc_fine=0; icc_fine<nCells_fine; ++icc_fine)
                        {
                            CutCell cc;
                            cc.parent = idx_crsn;

                            const VolIndex& vof_fine = vofs_fine[icc_fine];
                            const IntVect& iv_fine = vof_fine.gridIndex();
                            int idx_fine = vof_fine.cellIndex();

                            cidx[mfi].setVal(vof_fine,iv_fine,0,idx_fine);
                            for (int idir = 0; idir < SpaceDim; idir++)
                            {
                                for (SideIterator sit; sit.ok(); ++sit)
                                {
                                    int iside = sit() == Side::Lo ? 0 : 1;
                                    
                                    std::vector<FaceIndex> faces_crsn = ebgr_crsn.getFaces(vof_crsn,idir,sit());
                                    std::vector<FaceIndex> faces_fine = ebgr_fine.getFaces(vof_fine,idir,sit());

				    IntVect iv_nbr_fine = iv_fine + sign(sit())*BASISV(idir);
				    IntVect iv_face_fine = iv_fine + iside*BASISV(idir);

                                    int nValid = 0;
                                    for (int iface_fine=0; iface_fine<faces_fine.size(); ++iface_fine)
                                    {
                                        const FaceIndex& face_fine = faces_fine[iface_fine];
                                        int crsn_parent = -1;
                                        if (ebgr_fine.coarsen(face_fine.getVoF(Side::Lo)).gridIndex()
                                            != ebgr_fine.coarsen(face_fine.getVoF(Side::Hi)).gridIndex()) 
                                        {
                                            FaceIndex face_crsnd_fine = ebis_box_fine.coarsen(face_fine);
                                            for (int iface_crsn=0; iface_crsn<faces_crsn.size() && crsn_parent==-1; ++iface_crsn)
                                            {
                                                if (faces_crsn[iface_crsn] == face_crsnd_fine)
                                                {
                                                    crsn_parent = iface_crsn;
                                                }
                                            }
                                            BL_ASSERT(crsn_parent>=0);
                                        }

                                        if (!ebis_box_fine.isCovered(iv_nbr_fine))
                                        {
                                            if (ebis_box_fine.isRegular(iv_nbr_fine))
                                            {
                                                cc.nbr[idir][iside][nValid] = REGULAR_CELL;
                                            }
                                            else 
                                            {
                                                cc.nbr[idir][iside][nValid] = face_fine.cellIndex(sit());
                                                tf[idir][iv_face_fine].insert(tface(idx_fine,cc.nbr[idir][iside][nValid],sit,crsn_parent));
                                                fidx[idir][mfi].setVal(face_fine,iv_face_fine,0,nValid);
                                            }
                                            nValid++;
					    BL_ASSERT(nValid<NCELLMAX);
                                        }
                                    }
                                    cc.Nnbr[idir][iside] = nValid;
                                }
                            }
                            cmap_fab.setVal(cc,iv_fine,0,idx_fine);
                        }
                    }
                }
            }

            const Box& vbox = mfi.validbox();
            for (int idir = 0; idir < SpaceDim; idir++)
            {
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
                            cf.parent = fit->mParent;
                            bool is_face_between_reg_and_cut_cells
                                = cf.cellLo==REGULAR_CELL
                                || cf.cellHi==REGULAR_CELL;
                                
                            if ( !is_face_between_reg_and_cut_cells )
                            {
                                bool found_lo=false;
                                IntVect iv_lo = iv_face - BASISV(idir);
                                int key_lo = cmap_fab.key(iv_lo,0,cf.cellLo);
                                BL_ASSERT(key_lo >= 0);
                                CutCell& cc_lo = cmap_fab.getVal(iv_lo,0,cf.cellLo);
                                for(int inbr=0; inbr<cc_lo.Nnbr[idir][1]; ++inbr) {
                                    if (cc_lo.nbr[idir][1][inbr] == cf.cellHi) {
                                        cc_lo.faceID[idir][1][inbr] = cnt;
                                        found_lo = true;
                                    }
                                }

                                bool found_hi = false;
                                IntVect iv_hi = iv_face;
                                int key_hi = cmap_fab.key(iv_hi,0,cf.cellHi);
                                BL_ASSERT(key_hi >= 0);
                                CutCell& cc_hi = cmap_fab.getVal(iv_hi,0,cf.cellHi);
                                for(int inbr=0; inbr<cc_hi.Nnbr[idir][0]; ++inbr) {
                                    if (cc_hi.nbr[idir][0][inbr] == cf.cellLo) {
                                        cc_hi.faceID[idir][0][inbr] = cnt;
                                        found_hi = true;
                                    }
                                }

                                if (found_lo || found_hi) 
                                {
                                    fmap[idir][mfi].setVal(cf,iv_face,0,cnt++);
                                }
                            }
                            else
                            {

                            }
                        }
                    }
                }
            }
        }
    }

#if 0
    // Dump graphs
    std::cout << "cmap: " << std::endl;
    for (MFIter mfi(cmap); mfi.isValid(); ++mfi)
    {
        const BaseUmap<CutCell>& cmap_fab = cmap[mfi];
        for (BaseUmap<CutCell>::const_iterator umi = cmap_fab.begin(); umi<cmap_fab.end(); ++umi )
        {
            const CutCell& cc = *umi;
            const BaseUmap<CutCell>::Tuple& tuple = umi.tuple();
            std::cout << " " << tuple.pos << " parent = " << cc.parent
                      <<  " L = " << tuple.l << " ncomp= " << tuple.ncomp;
            for (int idir = 0; idir < SpaceDim; idir++)
            {
                for (SideIterator sit; sit.ok(); ++sit)
                {
                    int iside = sit() == Side::Lo ? 0 : 1;
                    int num = cc.Nnbr[idir][iside];
                    std::cout << " (" << idir << ", " << iside << "): " << num;
                    if (num > 0)
                    {
                        std::cout << " [ ";
                        for (int L=0; L<num; ++L)
                        {
                            std::cout << cc.nbr[idir][iside][L] << " ";
                        }
                        std::cout << "]";
                    }
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
                const BaseUmap<CutFace>::Tuple& tuple = umi.tuple();
                std::cout << " " << tuple.pos << " parent = " << cf.parent
                          << " L = " << tuple.l << " ncomp= " << tuple.ncomp
                          << " (" << idir << "): " << cf.cellLo << ", " << cf.cellHi << std::endl;
            }
        }
    }
#endif
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

typedef std::array<Real,BL_SPACEDIM>  RealDIM;
static void
Copy(RealDIM& out, const RealVect& in)
{
    for (int d=0; d<BL_SPACEDIM; ++d) {
        out[d] = in[d];
    }    
}

#if 0
static void
Copy(RealVect& out, const RealDIM& in)
{
    for (int d=0; d<BL_SPACEDIM; ++d) {
        out[d] = in[d];
    }    
}
#endif

extern "C"
{
    amrex_real fort_umap_norm (const int* lo, const int* hi,
                               const amrex_real* src, const int* src_sz,
                               const amrex::key_table_type* kt, const int* ktlo, const int* kthi, 
                               const int* max_mv, const int* ncomp,
                               const int* p);

    void do_eb_work(const int* lo, const int* hi,
                    const EBBndryData* bd, const int* Nbd,
                    const CutCell* cc, const int* Ncc,
                    const amrex::key_table_type* ktc, const int* ktclo, const int* ktchi, 
                    const int* cmax_mv, const int* cncomp);
}

int myTest()
{
    std::cout << "CutCell is POD: " << std::is_pod<CutCell>::value << std::endl;
    std::cout << "CutFace is POD: " << std::is_pod<CutFace>::value << std::endl;
    std::cout << "FaceData is POD: " << std::is_pod<FaceData>::value << std::endl;
    std::cout << "EBBndryData is POD: " << std::is_pod<EBBndryData>::value << std::endl;
    std::cout << "Size of CutCell: " << sizeof(CutCell)/sizeof(int) << " ints" << std::endl;
    std::cout << "Size of CutFace: " << sizeof(CutFace)/sizeof(int) << " ints" << std::endl;
    std::cout << "Size of FaceData: " << sizeof(FaceData)/sizeof(Real) << " reals" << std::endl;
    std::cout << "Size of EBBndryData: " << sizeof(EBBndryData)/sizeof(Real) << " reals" << std::endl;

    EBLevelGrid eblg_fine, eblg_crse, eblg_crsn;
    int ratio;
    get_EBLGs(eblg_fine,eblg_crse,ratio);
    const BoxArray& ba_fine = eblg_fine.getBoxArray();
    const DistributionMapping& dm_fine = eblg_fine.getDM();

    int nComp =1;
    int nGrow =4;
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

    BaseUmapFactory<VolIndex> cidx_factory(NCELLMAX);
    FabArray<BaseUmap<VolIndex> > cidx_fine(ba_fine, dm_fine, nComp, nGrow, MFInfo(), cidx_factory);

    BaseUmapFactory<FaceIndex> fidx_factory(NFACEMAX);
    std::array<FabArray<BaseUmap<FaceIndex> >, BL_SPACEDIM> fidx_fine;
    for (int idir=0; idir<BL_SPACEDIM; ++idir)
    {
        fidx_fine[idir].define(fba_fine[idir], dm_fine, nComp, nGrow, MFInfo(), fidx_factory);
    }



    BuildFortranGraph(cmap_fine,fmap_fine,cidx_fine,fidx_fine,eblg_fine,eblg_crsn,ratio);

    Real boundary_value = 10;
    Real internal_value = 5;

    // Build structure of EBBndryData to hold value, volume, and area, centroid and normal of boundary of cut cells
    BaseUmapFactory<EBBndryData> EBBD_factory(NCELLMAX);
    FabArray<BaseUmap<EBBndryData> > ebbd_fine(ba_fine, dm_fine, nComp, nGrow, MFInfo(), EBBD_factory);

    for (MFIter mfi(cidx_fine); mfi.isValid(); ++mfi)
    {
        const BaseUmap<VolIndex>& cidx_fab = cidx_fine[mfi];
        const EBISBox& ebis_box = eblg_fine.getEBISL()[mfi];
        BaseUmap<EBBndryData>& ebbd = ebbd_fine[mfi];
        for (BaseUmap<VolIndex>::const_iterator cit = cidx_fab.begin(); cit<cidx_fab.end(); ++cit )
        {
            const VolIndex& vof = *cit;
            const BaseUmap<VolIndex>::Tuple& tuple = cit.tuple();

            EBBndryData bd;
            Copy(bd.m_bndry_centroid,ebis_box.bndryCentroid(vof));
            Copy(bd.m_normal,ebis_box.normal(vof));
            bd.m_value = boundary_value;
            bd.m_bndry_area = ebis_box.bndryArea(vof);
            bd.m_vol_frac = ebis_box.volFrac(vof);
            ebbd.setVal(bd,tuple.pos,0,tuple.l);
        }
    }

    // Build structure of FaceData to hold aperature and centroid of faces
    BaseUmapFactory<FaceData> fd_factory(NFACEMAX);
    std::array<FabArray<BaseUmap<FaceData> >, BL_SPACEDIM> fd_fine;

    for (int idir=0; idir<BL_SPACEDIM; ++idir)
    {
        fd_fine[idir].define(fba_fine[idir], dm_fine, nComp, nGrow, MFInfo(), fd_factory);
    }

    for (int idir = 0; idir < SpaceDim; idir++)
    {
        for (MFIter mfi(fidx_fine[idir]); mfi.isValid(); ++mfi)
        {
            const BaseUmap<FaceIndex>& fidx_fab = fidx_fine[idir][mfi];
            const EBISBox& ebis_box = eblg_fine.getEBISL()[mfi];
            BaseUmap<FaceData>& fd = fd_fine[idir][mfi];
            for (BaseUmap<FaceIndex>::const_iterator fdi = fidx_fab.begin(); fdi<fidx_fab.end(); ++fdi )
            {
                const FaceIndex& fi = *fdi;
                const BaseUmap<FaceIndex>::Tuple& tuple = fdi.tuple();

                FaceData f;
                f.m_aperature = ebis_box.areaFrac(fi);
                Copy(f.m_centroid,ebis_box.centroid(fi));
                fd.setVal(f,tuple.pos,0,tuple.l);
            }
        }
    }

    // Build a structure to hold a Real value at each cell.  Set value to a constant.
    BaseUmapFactory<Real> CellData_factory(NCELLMAX);
    FabArray<BaseUmap<Real> > cd_fine(ba_fine, dm_fine, nComp, nGrow, MFInfo(), CellData_factory);

    for (MFIter mfi(cidx_fine); mfi.isValid(); ++mfi)
    {
        const BaseUmap<CutCell>& cmap_fab = cmap_fine[mfi];
        BaseUmap<Real>& c = cd_fine[mfi];
        for (BaseUmap<CutCell>::const_iterator cc = cmap_fab.begin(); cc<cmap_fab.end(); ++cc )
        {
            const BaseUmap<CutCell>::Tuple& tuple = cc.tuple();
            c.setVal(internal_value,tuple.pos,0,tuple.l);
        }
    }

    // Call a fortran routine to compute max norm of cell data over a BoxArray
    Real norm = 0;
    for (MFIter mfi(cd_fine); mfi.isValid(); ++mfi)
    {
        const BaseUmap<Real>& cd_fab = cd_fine[mfi];
        int npts = cd_fab.nPts();
        int max_mv = cd_fab.MaxMV();
        int ncomp = cd_fab.nComp();
        Real this_norm =  fort_umap_norm(ARLIM_3D(cd_fab.box().loVect()), ARLIM_3D(cd_fab.box().hiVect()),
                                         cd_fab.dataPtr(),&npts, 
                                         cd_fab.keyTablePtr(), ARLIM_3D(cd_fab.box().loVect()), ARLIM_3D(cd_fab.box().hiVect()),
                                         &max_mv, &ncomp, 0);
        norm = std::max(norm, this_norm);
    }
    ParallelDescriptor::ReduceRealMax(norm);
    std::cout << "Norm is " << norm << std::endl;


    for (MFIter mfi(ebbd_fine); mfi.isValid(); ++mfi)
    {
        const BaseUmap<EBBndryData>& ebbd_fab = ebbd_fine[mfi];
        const BaseUmap<CutCell>& cc_fab = cmap_fine[mfi];
        const Box& box = mfi.validbox();
        int Nbd = ebbd_fab.nPts();
        int Ncc = cc_fab.nPts();
        int max_mv = ebbd_fab.MaxMV();
        int ncomp = ebbd_fab.nComp();

        do_eb_work(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
                   ebbd_fab.dataPtr(), &Nbd,
                   cc_fab.dataPtr(),   &Ncc,
                   ebbd_fab.keyTablePtr(), ARLIM_3D(ebbd_fab.box().loVect()), ARLIM_3D(ebbd_fab.box().hiVect()), &max_mv, &ncomp);
    }




#if 1
    MultiFab vfrac(ba_fine, dm_fine,  nComp, 0);

    for(MFIter  mfi(vfrac); mfi.isValid(); ++mfi)
    {
        const EBISBox& ebis_box = eblg_fine.getEBISL()[mfi];
        const Box& vbox = mfi.validbox();
        FArrayBox& vfrac_fab = vfrac[mfi];
        vfrac[mfi].setVal(-1);
        if ( !(ebis_box.isAllRegular()) )
        {
            if (ebis_box.isCovered(vbox))
            {
                vfrac_fab.setVal(0);
            }
            else if (ebis_box.isRegular(vbox))
            {
                vfrac_fab.setVal(1);
            }
            else
            {
                // Only cut cells
                IntVectSet isIrreg = ebis_box.getIrregIVS(vbox);
                for (IVSIterator it(isIrreg); it.ok(); ++it)
                {
                    const IntVect& iv=it();
                    std::vector<VolIndex> vofs = ebis_box.getVoFs(iv);
                    Real vfrac_tot = 0;
                    for (int ivof=0; ivof<vofs.size(); ++ivof)
                    {
                        vfrac_tot += ebis_box.volFrac(vofs[ivof]);
                    }
                    vfrac_fab(iv,0) = vfrac_tot;
                }

                // Only regular cells
                for (IntVect iv(vbox.smallEnd()); iv<=vbox.bigEnd(); vbox.next(iv))
                {
                    if (ebis_box.isRegular(iv))
                    {
                        vfrac_fab(iv,0) = 1;                        
                    }
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
