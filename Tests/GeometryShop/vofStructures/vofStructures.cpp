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
#include "AMReX_iMultiFab.H"
#include "AMReX_BLProfiler.H"
#include <AMReX_BLFort.H>
using namespace amrex;

static const int NmvMAX = 4;
static const int COVERED_CELL = -1;
static const int REGULAR_CELL = -2;

//
// Structure to contain all necessary graph and geometry data for cut cells at an IntVect
//

typedef std::array<int,BL_SPACEDIM>  intDIM;

struct Node
{
    struct CutCellGraphNode
    {
        int Nnbr[BL_SPACEDIM][2];
        int nbr[BL_SPACEDIM][2][NmvMAX];
        int ebID;
    };
    int nCells;
    intDIM iv;
    CutCellGraphNode cells[NmvMAX];
};

extern "C"
{
    void do_eb_work(const int* lo, const int* hi,
                    Node* nodes, const int* num,
                    const BL_FORT_IFAB_ARG(mask));
}

static void
Copy(intDIM& out, const IntVect& in)
{
    for (int d=0; d<BL_SPACEDIM; ++d) {
        out[d] = in[d];
    }    
}

static void
Copy(IntVect& out, const intDIM& in)
{
    for (int d=0; d<BL_SPACEDIM; ++d) {
        out[d] = in[d];
    }    
}

#if 0
static void
dumpInts(const intDIM& iv)
{
    for (int d=0; d<BL_SPACEDIM; ++d) {
        std::cout << iv[d] << " ";
    }
    std::cout << '\n';
}
#endif

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
    //bool insideRegular = false;
    bool insideRegular = true;
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

int myTest()
{

    std::cout << "Is POD: " << std::is_pod<Node>::value << std::endl;
    std::cout << "Size: " << sizeof(Node)/sizeof(int) << " ints" << std::endl;

    EBLevelGrid eblg;
    get_EGLG(eblg);
    const BoxArray& ba = eblg.getBoxArray();
    const DistributionMapping& dm = eblg.getDM();

    IntVect tilesize(AMREX_D_DECL(10240,8,32));
    std::map<int,std::vector<Node> > graphNodes;

    int nGrow = 1;
    iMultiFab ebmask(ba,dm,1,nGrow); // Will contain location of Node in graphNodes vector

    for (MFIter mfi(ebmask); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        EBISBox ebis_box = eblg.getEBISL()[mfi];
        const Box gbox = ebmask[mfi].box() & ebis_box.getDomain();

        if (!ebis_box.isAllRegular())
        {
            const EBGraph& ebgr = ebis_box.getEBGraph();
            for (BoxIterator bit(gbox); bit.ok(); ++bit)
            {
                const IntVect iv = bit();
                if (ebis_box.isIrregular(iv))
                {
                    std::vector<VolIndex> gbox_vofs = ebis_box.getVoFs(iv);
                    if (gbox_vofs.size() > 0)
                    {
                        graphNodes[gid].push_back(Node());
                        Node& gn = graphNodes[gid].back();
                        Copy(gn.iv,iv);
                        gn.nCells = gbox_vofs.size();

                        for (int icc=0; icc<gn.nCells; ++icc)
                        {
                            const VolIndex& vof = gbox_vofs[icc];
                            for (int idir = 0; idir < SpaceDim; idir++)
                            {
                                for (SideIterator sit; sit.ok(); ++sit)
                                {
                                    int iside = sit() == Side::Lo ? 0 : 1;
                                    std::vector<FaceIndex> faces = ebgr.getFaces(vof,idir,sit());

                                    int nValid = 0;
                                    for (int iface=0; iface<faces.size(); ++iface)
                                    {
                                        IntVect iv_nbr = iv + sign(sit())*BASISV(idir);
                                        if (ebis_box.isRegular(iv_nbr))
                                        {
                                            gn.cells[icc].nbr[idir][iside][nValid++] = REGULAR_CELL;
                                        }
                                        else if (!ebis_box.isCovered(iv_nbr))
                                        {
                                            gn.cells[icc].nbr[idir][iside][nValid++] = faces[iface].cellIndex(sit());
                                        }
                                    }
                                    gn.cells[icc].Nnbr[idir][iside] = nValid;
                                }
                            }
                        }
                    }
                }
            }
        }

        // TODO: Fix up nbr ptrs to reference vector<Node> ordering

        BaseFab<int>& mask_fab = ebmask[mfi];
        mask_fab.setVal(REGULAR_CELL);
        for (int inode=0; inode<graphNodes[gid].size(); ++inode)
        {
            IntVect iv; 
            Copy(iv,graphNodes[gid][inode].iv);
            mask_fab(iv,0) = graphNodes[gid][inode].nCells;
        }
        if (ebis_box.isAllCovered()) {
            mask_fab.setVal(COVERED_CELL);
        }
        else if (!ebis_box.isAllRegular()) {
            // Not all covered, or all regular, there might be covered cells...
            for (BoxIterator bit(gbox); bit.ok(); ++bit)
            {
                if (ebis_box.isCovered(bit())) {
                    mask_fab(bit(),0) = COVERED_CELL;
                }
            }
        }

        const Box& vbox = mfi.validbox();
        int s = graphNodes[mfi.index()].size();
        do_eb_work(vbox.loVect(), vbox.hiVect(),
                   graphNodes[mfi.index()].data(), &s,
                   BL_TO_FORTRAN_N(mask_fab,0));


    }

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
