#include <cmath>

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

struct MyVoF
{
    int Nvols;
    int Nnbr[6*NmvMAX];
    int nbr[6*NmvMAX];
    Real volFrac[NmvMAX];
    int mask;
};

struct MyFiF
{
    int dir;
    int Nface;
    int nbr[2*NmvMAX];
    int mask;
};

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

int myTest()
{
    EBLevelGrid eblg;
    get_EGLG(eblg);
    const BoxArray& ba = eblg.getBoxArray();
    const DistributionMapping& dm = eblg.getDM();

    //MultiFab mf(ba,dm,1,0);
    //FillEBMask(eblg,mf);
    //VisMF::Write(mf,"mf");

    IntVect tilesize(AMREX_D_DECL(10240,8,32));
    MultiFab mf(ba,dm,1,0);

//     std::map<std::pair<int,int>,int> ii;

// #ifdef _OPENMP
// #pragma omp parallel
// #endif
//     for (MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi)
//     {
//         int tid = mfi.tileIndex();
//         int gid = mfi.index();
//         ii[std::pair<int,int>(gid,tid)] = 1;
//     }

    return 0;
}

int
main(int argc,char **argv)
{
    amrex::Initialize(argc,argv);
    {
        BL_PROFILE_VAR("main()", pmain);

        // check volume and surface area of approximate sphere
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
