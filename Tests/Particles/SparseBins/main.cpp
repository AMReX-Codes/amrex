#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SparseBins.H>

using namespace amrex;

struct TestParams
{
    IntVect size;
    int max_grid_size;
    int is_periodic;
    int nlevs;
};

void get_test_params(TestParams& params, const std::string& prefix)
{
    ParmParse pp(prefix);
    pp.get("size", params.size);
    pp.get("max_grid_size", params.max_grid_size);
    pp.get("is_periodic", params.is_periodic);
    pp.get("nlevs", params.nlevs);
}

void testIntersection()
{
    TestParams params;
    get_test_params(params, "intersect");
    
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = params.is_periodic;
    
    Vector<IntVect> rr(params.nlevs-1);
    for (int lev = 1; lev < params.nlevs; lev++)
        rr[lev-1] = IntVect(D_DECL(2,2,2));
    
    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++)
    {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    IntVect domain_lo(D_DECL(0 , 0, 0));
    IntVect domain_hi(D_DECL(params.size[0]-1, params.size[1]-1, params.size[2]-1));
    const Box base_domain(domain_lo, domain_hi);
    
    Vector<Geometry> geom(params.nlevs);
    geom[0].define(base_domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < params.nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, CoordSys::cartesian, is_per);
    }
    
    Vector<BoxArray> ba(params.nlevs);
    IntVect lo = IntVect(D_DECL(0, 0, 0));
    IntVect size = params.size;
    for (int lev = 0; lev < params.nlevs; ++lev)
    {
        ba[lev].define(Box(lo, lo+params.size-1));
        ba[lev].maxSize(params.max_grid_size);
        lo += size/2;
        size *= 2;
    }

    {
        auto num_boxes = ba[0].size();
        
    }
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    testIntersection();
    
    amrex::Finalize();
}
