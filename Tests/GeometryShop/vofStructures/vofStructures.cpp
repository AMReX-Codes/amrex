#include <cmath>
#include <type_traits>

#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_GeometryShop.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_ParmParse.H"
#include "AMReX_RealVect.H"
#include "AMReX_SphereIF.H"
#include "AMReX_EBGraph.H"
#include "AMReX_EBISBox.H"
#include "AMReX_MultiFab.H"
#include "AMReX_VisMF.H"
#include "AMReX_BLProfiler.H"
using namespace amrex;

static const int NmvMAX = 4;

//
// Structure to contain all necessary graph and geometry data for cut cells at an IntVect
//

typedef std::array<Real,BL_SPACEDIM> RealDIM;
typedef std::array<int,BL_SPACEDIM>  intDIM;

struct CutCell
{
    struct CutCellInfo
    {
        int      Nnbr[6];                    // [i] = number of neighbors in ith orientation
        int      nbr[6][NmvMAX];             // [i][j]  = index of jth neighbor in ith orientation
        Real     vol_frac;                   // volume fraction of this cut cull
        RealDIM  centroid;                   // coordinates of cell centroid
        RealDIM  face_centroid[6][NmvMAX];   // [i][j] coordinates of face_centroid of jth neighbor in ith orientation 
        Real     area_frac[6][NmvMAX];       // [i][j] area fraction of jth face in ith orientation
        RealDIM  bndry_centroid;             // coordinates of centroid of EB
        RealDIM  bndry_normal;               // normal vector of EB
        Real     bndry_area;                 // area of EB
        int      offset;                     // Offset into eb structure corresponding to this cut cell
    };
    int         Nvols;                       // number of cut cells at this IntVect
    CutCellInfo cells[NmvMAX];               // array of cut cells at this IntVect
    intDIM      iv;                          // IntVect where these cut cells live
};

static void
Copy(RealDIM& out, const RealVect& in)
{
    for (int d=0; d<BL_SPACEDIM; ++d) {
        out[d] = in[d];
    }    
}

static void
Copy(intDIM& out, const IntVect& in)
{
    for (int d=0; d<BL_SPACEDIM; ++d) {
        out[d] = in[d];
    }    
}

static void
dumpReals(const RealDIM& rv)
{
    for (int d=0; d<BL_SPACEDIM; ++d) {
        std::cout << rv[d] << " ";
    }
    std::cout << '\n';
}

static void
dumpInts(const intDIM& iv)
{
    for (int d=0; d<BL_SPACEDIM; ++d) {
        std::cout << iv[d] << " ";
    }
    std::cout << '\n';
}

void FillEBMask(const EBLevelGrid& eblg,
                MultiFab&          mf)
{
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        EBISBox ebis_box = eblg.getEBISL()[mfi];
        const EBGraph& ebg = ebis_box.getEBGraph();
        FArrayBox& fab = mf[mfi];
        const Box& gbox = fab.box();
        BoxIterator bit(gbox);
        for (bit.reset(); bit.ok(); ++bit)
        {
            Real this_val = -1; // In fluid, by default
            const IntVect& iv = bit();
            if (ebg.isIrregular(iv)) {
                this_val = 0;
            }
            else if (ebg.isCovered(iv)) {
                this_val = 1;
            }
            fab(iv,0) = this_val;
        }
    }
}

void
get_EGLG(EBLevelGrid& eblg)
{
    Real radius = 0.5;
    Real domlen = 1;
    std::vector<Real> centervec(SpaceDim);
    std::vector<int>  ncellsvec(SpaceDim);

    ParmParse pp;
    pp.getarr(  "n_cell"       , ncellsvec, 0, SpaceDim);
    pp.get(   "sphere_radius", radius);
    pp.getarr("sphere_center", centervec, 0, SpaceDim);
    pp.get("domain_length", domlen);                     
    RealVect center;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
        center[idir] = centervec[idir];
    }
    bool insideRegular = false;
    SphereIF sphere(radius, center, insideRegular);
    int verbosity = 0;

    pp.get("verbosity", verbosity);
    GeometryShop gshop(sphere, verbosity);
    BaseFab<int> regIrregCovered;
    std::vector<IrregNode> nodes;

    IntVect ivlo = IntVect::TheZeroVector();
    IntVect ivhi;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
        ivhi[idir] = ncellsvec[idir] - 1;
    }


    Box domain(ivlo, ivhi);
    RealVect origin = RealVect::Zero;
    Real dx = domlen/ncellsvec[0];

    EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
    ebisPtr->define(domain, origin, dx, gshop);

    BoxArray ba(domain);
    int maxboxsize = 16;
    pp.query("maxboxsize",maxboxsize);
    ba.maxSize(maxboxsize);
    DistributionMapping dm(ba);
    int nGrowEBLG = 2;
    eblg.define(ba, dm, domain, nGrowEBLG);
}

struct TC
{
    struct TB
    {
        int c;
    };
    std::array<int,3> b;
    Real a;
};

int myTest()
{

    std::cout << "Is POD: " << std::is_pod<CutCell>::value << std::endl;

#if 1
    EBLevelGrid eblg;
    get_EGLG(eblg);
    const BoxArray& ba = eblg.getBoxArray();
    const DistributionMapping& dm = eblg.getDM();

    MultiFab ebmask(ba,dm,1,0);
    FillEBMask(eblg,ebmask);
    VisMF::Write(ebmask,"ebmask");

    IntVect tilesize(AMREX_D_DECL(10240,8,32));
    MultiFab mf(ba,dm,1,0);
    std::map<int,std::vector<CutCell> > cutCells;

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        const Box& gbox = mf[mfi].box();
        EBISBox ebis_box = eblg.getEBISL()[mfi];

        int offset = 0;

        if (!ebis_box.isAllRegular())
        {
            const EBGraph& ebgr = ebis_box.getEBGraph();
            for (BoxIterator bit(gbox); bit.ok(); ++bit)
            {
                const IntVect iv = bit();
                if (ebis_box.isIrregular(iv)) {
                    std::vector<VolIndex> gbox_vofs = ebis_box.getVoFs(iv);
                    if (gbox_vofs.size() > 0)
                    {
                        cutCells[gid].push_back(CutCell());
                        CutCell& cc = cutCells[gid].back();
                        cc.Nvols = gbox_vofs.size();
                        Copy(cc.iv,iv);

                        for (int icc=0; icc<cc.Nvols; ++icc)
                        {
                            const VolIndex& vof = gbox_vofs[icc];
                            cc.cells[icc].offset = offset + icc;
                            cc.cells[icc].vol_frac = ebis_box.volFrac(vof);
                            Copy(cc.cells[icc].centroid, ebis_box.centroid(vof));
                            Copy(cc.cells[icc].bndry_centroid, ebis_box.bndryCentroid(vof));
                            Copy(cc.cells[icc].bndry_normal, ebis_box.normal(vof));
                            cc.cells[icc].bndry_area = ebis_box.bndryArea(vof);

                            for (int idir = 0; idir < SpaceDim; idir++)
                            {
                                for (SideIterator sit; sit.ok(); ++sit)
                                {
                                    int idx1D = IrregNode::index(idir, sit());
                                    std::vector<FaceIndex> faces = ebgr.getFaces(vof,idir,sit());
                                    cc.cells[icc].Nnbr[idx1D] = faces.size();
                                    for (int iface=0; iface<faces.size(); ++iface)
                                    {
                                        cc.cells[icc].nbr[idx1D][iface] = faces[iface].cellIndex(sit());
                                        Copy(cc.cells[icc].face_centroid[idx1D][iface], ebis_box.centroid(faces[iface]));
                                        cc.cells[icc].area_frac[idx1D][iface] = ebis_box.areaFrac(faces[iface]);
                                    }
                                }
                            }

                        }
                        offset += cc.Nvols;
#if 1
                        std::cout << "MyVOF: \n";
                        std::cout << "  Nvols: " << cc.Nvols << '\n';
                        std::cout << "  iv: "; dumpInts(cc.iv); 
                        std::cout << "  offset:  ";
                        for (int icc=0; icc<cc.Nvols; ++icc) {
                            std::cout << cc.cells[icc].offset << " ";
                        }
                        std::cout << '\n';
                        for (int icc=0; icc<cc.Nvols; ++icc) {
                            std::cout << "  i,vol: " << icc << "," << cc.cells[icc].vol_frac << '\n';
                            std::cout << "    centroid: "; dumpReals(cc.cells[icc].centroid);
                            std::cout << "    bndry_centroid: "; dumpReals(cc.cells[icc].bndry_centroid);
                            std::cout << "    bndry_normal: "; dumpReals(cc.cells[icc].bndry_normal);
                            std::cout << "    bndry_area: " << cc.cells[icc].bndry_area << '\n';
                            std::cout << "    FACES:\n";
                            for (SideIterator sit; sit.ok(); ++sit) {
                                for (int idir = 0; idir < SpaceDim; idir++) {
                                    std::string side = sit()==Side::Lo ? "LO" : "HI";
                                    int d = IrregNode::index(idir, sit());
                                    std::cout << "       d: " << d << " (dir,face = " << idir << ", " << side << ")\n";
                                    if (cc.cells[icc].Nnbr[d] == 0) {
                                        std::cout << "              NO FACES \n\n";
                                    } else {
                                        for (int j=0; j<cc.cells[icc].Nnbr[d]; ++j) {
                                            std::cout << "          nbr = " << cc.cells[icc].nbr[d][j];
                                            std::cout << " face area: " << cc.cells[icc].area_frac[d][j];
                                            std::cout << " face cent: "; dumpReals(cc.cells[icc].face_centroid[d][j]);
                                            std::cout << std::endl;
                                        }
                                    }
                                }
                            }
                        }
#endif

                    }
                }
            }
        }
    }
#endif
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
