#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

#include "CheckPair.H"

#include "MDParticleContainer.H"

#include <string>

namespace amrex
{
    template <typename T, typename S> 
    std::ostream& operator<<(std::ostream& os, const std::pair<T, S>& v) 
    { 
        os << "("; 
        os << v.first << ", " 
           << v.second << ")"; 
        return os; 
    } 
}

using namespace amrex;

struct TestParams
{
    IntVect size;
    int max_grid_size;
    int num_ppc;
    int is_periodic;
};

void testNeighborParticles();

void testNeighborList();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    amrex::PrintToFile("neighbor_test") << "Running neighbor particles test \n";
    testNeighborParticles();

    amrex::PrintToFile("neighbor_test") << "Running neighbor list test \n";
    testNeighborList();

    amrex::Finalize();
}

void get_test_params(TestParams& params, const std::string& prefix)
{
    ParmParse pp(prefix);
    pp.get("size", params.size);
    pp.get("max_grid_size", params.max_grid_size);
    pp.get("num_ppc", params.num_ppc);
    pp.get("is_periodic", params.is_periodic);
}

void testNeighborParticles ()
{
    BL_PROFILE("testNeighborParticles");
    TestParams params;
    get_test_params(params, "nbor_parts");

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, params.size[n]);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(params.size[0]-1,params.size[1]-1,params.size[2]-1));
    const Box domain(domain_lo, domain_hi);

    int coord = 0;
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = params.is_periodic;
    Geometry geom(domain, &real_box, coord, is_per);
    
    BoxArray ba(domain);
    ba.maxSize(params.max_grid_size);
    DistributionMapping dm(ba);

    const int ncells = 1;
    MDParticleContainer pc(geom, dm, ba, ncells);

    int npc = params.num_ppc;
    IntVect nppc = IntVect(AMREX_D_DECL(npc, npc, npc));

    if (ParallelDescriptor::MyProc() == dm[0])
        amrex::PrintToFile("neighbor_test") << "About to initialize particles \n";

    pc.InitParticles(nppc, 1.0, 0.0);

    if (ParallelDescriptor::MyProc() == dm[0])
        amrex::PrintToFile("neighbor_test") << "Check neighbors after init ... \n";
    pc.checkNeighborParticles();

    pc.fillNeighbors();

    if (ParallelDescriptor::MyProc() == dm[0])
        amrex::PrintToFile("neighbor_test") << "Check neighbors after fill ... \n";
    pc.checkNeighborParticles();

    pc.updateNeighbors();

    if (ParallelDescriptor::MyProc() == dm[0])
        amrex::PrintToFile("neighbor_test") << "Check neighbors after update ... \n";
    pc.checkNeighborParticles();

    if (ParallelDescriptor::MyProc() == dm[0])
        amrex::PrintToFile("neighbor_test") << "Now resetting the particle test_id values  \n";
    pc.reset_test_id();

    if (ParallelDescriptor::MyProc() == dm[0])
        amrex::PrintToFile("neighbor_test") << "Check neighbors after reset ... \n";
    pc.checkNeighborParticles();

    if (ParallelDescriptor::MyProc() == dm[0])
        amrex::PrintToFile("neighbor_test") << "Now updateNeighbors again ...  \n";
    pc.updateNeighbors();

    if (ParallelDescriptor::MyProc() == dm[0])
        amrex::PrintToFile("neighbor_test") << "Check neighbors after update ... \n";
    pc.checkNeighborParticles();

    ParallelDescriptor::Barrier();

    amrex::PrintToFile("neighbor_test") << "Testing neighbor particles after move \n";

    // so we can call minAndMaxDistance
    pc.buildNeighborList(CheckPair());

    amrex::PrintToFile("neighbor_test") << "Min distance is " << pc.minAndMaxDistance() << ", should be (1, 1) \n"; 

    amrex::PrintToFile("neighbor_test") << "Moving particles and updating neighbors \n";
    pc.moveParticles(0.1);
    pc.updateNeighbors();

    amrex::PrintToFile("neighbor_test") << "Min distance is " << pc.minAndMaxDistance() << ", should be (1, 1) \n"; 

    amrex::PrintToFile("neighbor_test") << "Moving particles and updating neighbors again \n";
    pc.moveParticles(0.1);
    pc.updateNeighbors();

    amrex::PrintToFile("neighbor_test") << "Min distance is " << pc.minAndMaxDistance() << ", should be (1, 1) \n"; 

    amrex::PrintToFile("neighbor_test") << "Moving particles and updating neighbors yet again \n";
    pc.moveParticles(0.1);
    pc.updateNeighbors();

    amrex::PrintToFile("neighbor_test") << "Min distance is " << pc.minAndMaxDistance() << ", should be (1, 1) \n";     
}

void testNeighborList ()
{
    BL_PROFILE("main::main()");
    TestParams params;
    get_test_params(params, "nbor_list");

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, params.size[n]);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(params.size[0]-1,params.size[1]-1,params.size[2]-1));
    const Box domain(domain_lo, domain_hi);

    int coord = 0;
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = params.is_periodic;
    Geometry geom(domain, &real_box, coord, is_per);
    
    BoxArray ba(domain);
    ba.maxSize(params.max_grid_size);
    DistributionMapping dm(ba);

    const int ncells = 1;
    MDParticleContainer pc(geom, dm, ba, ncells);

    int npc = params.num_ppc;
    IntVect nppc = IntVect(AMREX_D_DECL(npc, npc, npc));

    amrex::PrintToFile("neighbor_test") << "About to initialize particles" << std::endl;

    pc.InitParticles(nppc, 1.0, 0.0);
    pc.fillNeighbors();

    pc.buildNeighborList(CheckPair());

    pc.checkNeighborList();
}
