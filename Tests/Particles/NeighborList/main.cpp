#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <string>

using namespace amrex;

struct TestParams
{
    IntVect size;
    int max_grid_size;
    int is_periodic;
};

void testNeighborList();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    testNeighborList();

    amrex::Finalize();
}

void get_test_params(TestParams& params, const std::string& prefix)
{
    ParmParse pp(prefix);
    pp.get("size", params.size);
    pp.get("max_grid_size", params.max_grid_size);
    pp.get("is_periodic", params.is_periodic);
}

namespace Params
{
    static constexpr amrex::Real cutoff = 0.2;
}

struct CheckPair
{
    template <class P1, class P2>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    bool operator() (const P1& p1, const P2& p2) const
    {
        AMREX_D_TERM(amrex::Real d0 = (p1.pos(0) - p2.pos(0));,
                     amrex::Real d1 = (p1.pos(1) - p2.pos(1));,
                     amrex::Real d2 = (p1.pos(2) - p2.pos(2));)
        amrex::Real dsquared = AMREX_D_TERM(d0*d0, + d1*d1, + d2*d2);
        return (dsquared <= 25.0*Params::cutoff*Params::cutoff);
    }
};

using PCType1 = amrex::NeighborParticleContainer<0, 0>;
using PCType2 = amrex::NeighborParticleContainer<0, 1>;
using PType1 = PCType1::ParticleType;
using PType2 = PCType2::ParticleType;
using MyParIter = PCType1::ParIterType;
using NeighborListContainer1 = Vector<std::map<std::pair<int, int>, amrex::NeighborList<PType1> > >;
using NeighborListContainer2 = Vector<std::map<std::pair<int, int>, amrex::NeighborList<PType2> > >;

void addParticles (PCType1& pc1, PCType2& pc2)
{
    auto& ptile1 = pc1.DefineAndReturnParticleTile(0, 0, 0);
    {
        PType1 p;
        AMREX_D_TERM(p.pos(0) = 12.0;,
                     p.pos(1) = 12.0;,
                     p.pos(2) = 12.0;)
        p.id() = PType1::NextID();
        p.cpu() = ParallelDescriptor::MyProc();
        ptile1.push_back(p);
    }

    auto& ptile2 = pc2.DefineAndReturnParticleTile(0, 0, 0);
    {
        PType2 p;
        AMREX_D_TERM(p.pos(0) = 12.0;,
                     p.pos(1) = 13.0;,
                     p.pos(2) = 12.0;)
        p.id() = PType2::NextID();
        p.cpu() = ParallelDescriptor::MyProc();
        p.idata(0) = 1;
        ptile2.push_back(p);
    }

    pc1.Redistribute();
    pc2.Redistribute();
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
    int is_per[] = {AMREX_D_DECL(params.is_periodic,
                                 params.is_periodic,
                                 params.is_periodic)};
    Geometry geom(domain, &real_box, coord, is_per);

    BoxArray ba(domain);
    ba.maxSize(params.max_grid_size);
    DistributionMapping dm(ba);

    const int ncells = 1;
    PCType1 pc1(geom, dm, ba, ncells);
    PCType2 pc2(geom, dm, ba, ncells);

    addParticles(pc1, pc2);

    pc1.fillNeighbors();
    pc2.fillNeighbors();

    NeighborListContainer1 neighbors_type_1;
    NeighborListContainer2 neighbors_type_2;
    pc1.buildNeighborList(CheckPair(), pc1, neighbors_type_1);
    pc1.buildNeighborList(CheckPair(), pc2, neighbors_type_2);

    for (MyParIter pti(pc1, 0); pti.isValid(); ++pti) {
        int gid = pti.index();
        int tid = pti.LocalTileIndex();
        auto& nlist1 = neighbors_type_1[0][std::make_pair(gid, tid)];
        nlist1.print();
        auto& nlist2 = neighbors_type_2[0][std::make_pair(gid, tid)];
        nlist2.print();
    }
}
