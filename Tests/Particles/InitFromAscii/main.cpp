#include <cstdlib>
#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;

struct TestParams {
  int nx;
  int ny;
  int nz;
  int max_grid_size;
  bool verbose;
};

void test_init_ascii (TestParams& parms)
{
    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++)
    {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1));
    const Box domain(domain_lo, domain_hi);

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[] = {AMREX_D_DECL(1,1,1)};
    Geometry geom(domain, &real_box, CoordSys::cartesian, is_per);

    BoxArray ba(domain);
    ba.maxSize(parms.max_grid_size);

    DistributionMapping dmap(ba);

    using MyParticleContainer = ParticleContainer<4, 0>;
    MyParticleContainer myPC(geom, dmap, ba);

    myPC.InitFromAsciiFile("particles.txt", 1 + AMREX_SPACEDIM);

    // should be 8
    AMREX_ALWAYS_ASSERT(myPC.TotalNumberOfParticles() == 8);

    // should be 8010.0
    using PType = MyParticleContainer::SuperParticleType;
    auto tot = amrex::ReduceSum(myPC,
                                [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                                {
                                    Real total = 0.0;
                                    for (int i = 0; i < 1 + AMREX_SPACEDIM; ++i)
                                    {
                                        total += p.rdata(i);
                                    }
                                    return total;
                                });

    AMREX_ALWAYS_ASSERT(tot == 8010);

    // The same file, with the extra values split between the struct, a compile-time
    // SoA component and a runtime one. The int components are not in the file and
    // should be zero.
    using SoAParticleContainer = ParticleContainer<2, 0, 1, 1>;
    SoAParticleContainer soaPC(geom, dmap, ba);
    soaPC.AddRealComp("runtime_real");
    soaPC.AddIntComp("runtime_int");

    soaPC.InitFromAsciiFile("particles.txt", 1 + AMREX_SPACEDIM);

    AMREX_ALWAYS_ASSERT(soaPC.TotalNumberOfParticles() == 8);

    // sums of the four extra columns in particles.txt
    Vector<Real> const expected{8000, -1, 4, 7};
    Vector<Real> sums(4, 0);
    Vector<Long> int_sums(2, 0);
    for (int lev = 0; lev <= soaPC.finestLevel(); ++lev) {
        for (auto& kv : soaPC.GetParticles(lev)) {
            auto const& ptile = kv.second;
            auto const np = ptile.numParticles();
            auto const& soa = ptile.GetStructOfArrays();
            AMREX_ALWAYS_ASSERT(soa.NumRealComps() == 2 && soa.NumIntComps() == 2);

            Vector<SoAParticleContainer::ParticleType> aos(np);
            Gpu::copy(Gpu::deviceToHost, ptile.GetArrayOfStructs().begin(),
                      ptile.GetArrayOfStructs().end(), aos.begin());
            for (auto const& p : aos) {
                sums[0] += p.rdata(0);
                sums[1] += p.rdata(1);
            }
            for (int comp = 0; comp < 2; ++comp) {
                AMREX_ALWAYS_ASSERT(Long(soa.GetRealData(comp).size()) == np);
                AMREX_ALWAYS_ASSERT(Long(soa.GetIntData(comp).size()) == np);
                Vector<ParticleReal> r(np);
                Gpu::copy(Gpu::deviceToHost, soa.GetRealData(comp).begin(),
                          soa.GetRealData(comp).end(), r.begin());
                for (auto v : r) { sums[2+comp] += v; }
                Vector<int> iv(np);
                Gpu::copy(Gpu::deviceToHost, soa.GetIntData(comp).begin(),
                          soa.GetIntData(comp).end(), iv.begin());
                for (auto v : iv) { int_sums[comp] += std::abs(v); }
            }
        }
    }
    ParallelDescriptor::ReduceRealSum(sums.data(), int(sums.size()));
    ParallelDescriptor::ReduceLongSum(int_sums.data(), int(int_sums.size()));
    for (int i = 0; i < 4; ++i) {
        AMREX_ALWAYS_ASSERT(sums[i] == expected[i]);
    }
    AMREX_ALWAYS_ASSERT(int_sums[0] == 0 && int_sums[1] == 0);
}

int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);

  ParmParse pp;

  TestParams parms;
  pp.get("nx", parms.nx);
  pp.get("ny", parms.ny);
  pp.get("nz", parms.nz);
  pp.get("max_grid_size", parms.max_grid_size);

  test_init_ascii(parms);

  amrex::Finalize();
}
