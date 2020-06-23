#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_PlotFileUtil.H>

#include <thread>

using namespace amrex;

static constexpr int NSR = 1 + AMREX_SPACEDIM;
static constexpr int NSI = 0;
static constexpr int NAR = 0;
static constexpr int NAI = 0;

struct TestParams {
  int nx;
  int ny;
  int nz;
  int max_grid_size;
  int nppc;
  int nlevs;
  bool verbose;
};

void get_position_unit_cell(Real* r, const IntVect& nppc, int i_part)
{
    int nx = nppc[0];
#if AMREX_SPACEDIM > 1
    int ny = nppc[1];
#else
    int ny = 1;
#endif
#if AMREX_SPACEDIM > 2
    int nz = nppc[2];
#else
    int nz = 1;
#endif

    int ix_part = i_part/(ny * nz);
    int iy_part = (i_part % (ny * nz)) % ny;
    int iz_part = (i_part % (ny * nz)) / ny;

    r[0] = (0.5+ix_part)/nx;
    r[1] = (0.5+iy_part)/ny;
    r[2] = (0.5+iz_part)/nz;
}

class MyParticleContainer
    : public amrex::ParticleContainer<NSR, NSI, NAR, NAI>
{

public:

    MyParticleContainer (const Vector<amrex::Geometry>            & a_geom,
                         const Vector<amrex::DistributionMapping> & a_dmap,
                         const Vector<amrex::BoxArray>            & a_ba,
                         const Vector<amrex::IntVect>             & a_rr)
        : amrex::ParticleContainer<NSR, NSI, NAR, NAI>(a_geom, a_dmap, a_ba, a_rr)
    {}

    void InitParticles (const amrex::IntVect& a_num_particles_per_cell)
    {
        BL_PROFILE("InitParticles");

        const int lev = 0;  // only add particles on level 0
        const Real* dx = Geom(lev).CellSize();
        const Real* plo = Geom(lev).ProbLo();

        const int num_ppc = AMREX_D_TERM( a_num_particles_per_cell[0],
                                         *a_num_particles_per_cell[1],
                                         *a_num_particles_per_cell[2]);

        for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            const Box& tile_box  = mfi.tilebox();

            Gpu::HostVector<ParticleType> host_particles;
            std::array<Gpu::HostVector<Real>, NAR> host_real;
            std::array<Gpu::HostVector<int>, NAI> host_int;

            std::vector<Gpu::HostVector<Real> > host_runtime_real(NumRuntimeRealComps());
            std::vector<Gpu::HostVector<int> > host_runtime_int(NumRuntimeIntComps());

            for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
            {
                for (int i_part=0; i_part<num_ppc;i_part++) {
                    Real r[3];
                    get_position_unit_cell(r, a_num_particles_per_cell, i_part);

                    ParticleType p;
                    p.id()  = ParticleType::NextID();
                    p.cpu() = ParallelDescriptor::MyProc();
                    p.pos(0) = plo[0] + (iv[0] + r[0])*dx[0];
#if AMREX_SPACEDIM > 1
                    p.pos(1) = plo[1] + (iv[1] + r[1])*dx[1];
#endif
#if AMREX_SPACEDIM > 2
                    p.pos(2) = plo[2] + (iv[2] + r[2])*dx[2];
#endif

                    for (int i = 0; i < NSR; ++i) p.rdata(i) = p.id();
                    for (int i = 0; i < NSI; ++i) p.idata(i) = p.id();

                    host_particles.push_back(p);
                    for (int i = 0; i < NAR; ++i)
                        host_real[i].push_back(p.id());
                    for (int i = 0; i < NAI; ++i)
                        host_int[i].push_back(p.id());
                    for (int i = 0; i < NumRuntimeRealComps(); ++i)
                        host_runtime_real[i].push_back(p.id());
                    for (int i = 0; i < NumRuntimeIntComps(); ++i)
                        host_runtime_int[i].push_back(p.id());
                }
            }

            auto& particle_tile = DefineAndReturnParticleTile(lev, mfi.index(), mfi.LocalTileIndex());
            auto old_size = particle_tile.GetArrayOfStructs().size();
            auto new_size = old_size + host_particles.size();
            particle_tile.resize(new_size);

            Gpu::copy(Gpu::hostToDevice,
                      host_particles.begin(),
                      host_particles.end(),
                      particle_tile.GetArrayOfStructs().begin() + old_size);

            auto& soa = particle_tile.GetStructOfArrays();
            for (int i = 0; i < NAR; ++i)
            {
                Gpu::copy(Gpu::hostToDevice,
                          host_real[i].begin(),
                          host_real[i].end(),
                          soa.GetRealData(i).begin() + old_size);
            }

            for (int i = 0; i < NAI; ++i)
            {
                Gpu::copy(Gpu::hostToDevice,
                          host_int[i].begin(),
                          host_int[i].end(),
                          soa.GetIntData(i).begin() + old_size);
            }
            for (int i = 0; i < NumRuntimeRealComps(); ++i)
            {
                Gpu::copy(Gpu::hostToDevice,
                          host_runtime_real[i].begin(),
                          host_runtime_real[i].end(),
                          soa.GetRealData(NAR+i).begin() + old_size);
            }

            for (int i = 0; i < NumRuntimeIntComps(); ++i)
            {
                Gpu::copy(Gpu::hostToDevice,
                          host_runtime_int[i].begin(),
                          host_runtime_int[i].end(),
                          soa.GetIntData(NAI+i).begin() + old_size);
            }
        }

        Redistribute();
    }
};

void test_async_io(TestParams& parms)
{
    int nlevs = parms.nlevs;

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    RealBox fine_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
       fine_box.setLo(n,0.25);
       fine_box.setHi(n,0.75);
    }

    IntVect domain_lo(D_DECL(0 , 0, 0));
    IntVect domain_hi(D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1));
    const Box domain(domain_lo, domain_hi);

    // Define the refinement ratio
    Vector<IntVect> rr(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++)
        rr[lev-1] = IntVect(AMREX_D_DECL(2, 2, 2));

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = 1;

    // This defines a Geometry object which is useful for writing the plotfiles
    Vector<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
	geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
			 &real_box, CoordSys::cartesian, is_per);
    }

    Vector<BoxArray> ba(nlevs);
    ba[0].define(domain);

    if (nlevs > 1) {
        int n_fine = parms.nx*rr[0][0];
        IntVect refined_lo(D_DECL(n_fine/4,n_fine/4,n_fine/4));
        IntVect refined_hi(D_DECL(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba[1].define(refined_patch);
    }

    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(parms.max_grid_size);
    }

    Vector<DistributionMapping> dmap(nlevs);

    Vector<std::unique_ptr<MultiFab> > partMF(nlevs);
    Vector<std::unique_ptr<MultiFab> > density(nlevs);
    Vector<std::unique_ptr<MultiFab> > acceleration(nlevs);
    for (int lev = 0; lev < nlevs; lev++) {
        dmap[lev] = DistributionMapping{ba[lev]};
        density[lev].reset(new MultiFab(ba[lev], dmap[lev], 1, 0));
        density[lev]->setVal(0.0);
        acceleration[lev].reset(new MultiFab(ba[lev], dmap[lev], 3, 1));
        acceleration[lev]->setVal(5.0, 1);
    }

    MyParticleContainer myPC(geom, dmap, ba, rr);
    myPC.SetVerbose(false);

    myPC.InitParticles(IntVect(2, 2, 2));

    for (int step = 0; step < 4000; ++step)
    {
        myPC.AssignDensity(0, partMF, 0, 1, nlevs-1);

        for (int lev = 0; lev < nlevs; ++lev) {
            MultiFab::Copy(*density[lev], *partMF[lev], 0, 0, 1, 0);
        }

        if (step % 1000 == 0) {
            Vector<std::string> varnames;
            varnames.push_back("density");

            Vector<std::string> particle_varnames;
            particle_varnames.push_back("mass");

            Vector<int> level_steps;
            level_steps.push_back(0);
            level_steps.push_back(0);

            int output_levs = nlevs;

            Vector<const MultiFab*> outputMF(output_levs);
            Vector<IntVect> outputRR(output_levs);
            for (int lev = 0; lev < output_levs; ++lev) {
                outputMF[lev] = density[lev].get();
                outputRR[lev] = IntVect(D_DECL(2, 2, 2));
            }

            std::string fn = amrex::Concatenate("plt", step, 5);

            WriteMultiLevelPlotfile(fn, output_levs, outputMF,
                                    varnames, geom, 0.0, level_steps, outputRR);

            myPC.WritePlotFile(fn, "particle0");
        }
    }
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
  pp.get("nlevs", parms.nlevs);
  pp.get("nppc", parms.nppc);
  if (parms.nppc < 1 && ParallelDescriptor::IOProcessor())
    amrex::Abort("Must specify at least one particle per cell");

  parms.verbose = false;
  pp.query("verbose", parms.verbose);

  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << std::endl;
    std::cout << "Number of particles per cell : ";
    std::cout << parms.nppc  << std::endl;
    std::cout << "Size of domain               : ";
    std::cout << "Num levels: ";
    std::cout << parms.nlevs << std::endl;
    std::cout << parms.nx << " " << parms.ny << " " << parms.nz << std::endl;
  }

  test_async_io(parms);

  amrex::Finalize();
}
