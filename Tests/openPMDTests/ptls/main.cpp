#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"
#include "AMReX_PlotFileUtil.H"
#include <AMReX_AmrParticles.H>

#include "mypc.H"
#include "trilinear_deposition_K.H"

#include "warpxBTD.H"
#include "warpxWriter.H"

using namespace amrex;
#define TEST_BTD 1

struct TestParams {
    int nx;
    int ny;
    int nz;
    int nlevs;
    int max_grid_size;
    int nppc;
    bool verbose;

    int nplotfile=1;
};

void checkMFBox(const TestParams& parms,
                Vector<const MultiFab*> outputMF)
{
      for (int lev=0; lev < parms.nlevs; lev++)
        {
          auto const* curr_mf = outputMF[lev];
          int const ncomp = curr_mf->nComp();
          amrex::Print()<<" checking boxes, lev="<<lev<<" ncomp="<<ncomp<<std::endl;

          for ( int icomp=0; icomp<ncomp; icomp++ )
            {
              for( amrex::MFIter mfi(*curr_mf); mfi.isValid(); ++mfi )
                {
                  amrex::FArrayBox const& fab = (*curr_mf)[mfi];
                  amrex::Box const& local_box = fab.box();

                  amrex::Print()<<"  .. checking:   icomp="<<icomp<< " local box="<<local_box;
                  amrex::Print()<<"   "<<*(fab.dataPtr())<<std::endl;
                }
            }
        }
}



void testParticleMesh (TestParams& parms, int nghost)
{
    Vector<IntVect> rr(parms.nlevs-1);
    for (int lev = 1; lev < parms.nlevs; lev++) {
        rr[lev-1] = IntVect(AMREX_D_DECL(2,2,2));
    }

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1));
    const Box base_domain(domain_lo, domain_hi);

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM] = {AMREX_D_DECL(1,1,1)};

    Vector<Geometry> geom(parms.nlevs);
    geom[0].define(base_domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < parms.nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, CoordSys::cartesian, is_per);
    }

    Vector<BoxArray> ba(parms.nlevs);
    Vector<DistributionMapping> dm(parms.nlevs);

    Box domain = base_domain;
    IntVect size = IntVect(AMREX_D_DECL(parms.nx, parms.ny, parms.nz));
    for (int lev = 0; lev < parms.nlevs; ++lev)
    {
        ba[lev].define(domain);
        ba[lev].maxSize(parms.max_grid_size);
        dm[lev].define(ba[lev]);
        domain.grow(-size/4);   // fine level cover the middle of the coarse domain
        domain.refine(2);
    }

    Vector<MultiFab> density1(parms.nlevs);
    Vector<MultiFab> density2(parms.nlevs);

    for (int lev = 0; lev < parms.nlevs; lev++) {
        // one field comp for density1
        density1[lev].define(ba[lev], dm[lev], 1, nghost);
        density1[lev].setVal(0.0);
        // and two field comp for density2
        density2[lev].define(ba[lev], dm[lev], 2, nghost);
        density2[lev].setVal(2.0);
    }

    MyParticleContainer myPC(geom, dm, ba, rr);
    myPC.SetVerbose(false);

    bool serialize = true;
    if (ParallelDescriptor::NProcs() > 1) {
        serialize = false;
    }

    amrex::Long num_particles = (amrex::Long)parms.nppc * parms.nx * parms.ny * parms.nz;
    amrex::Print() << serialize<< " Total number of particles   :" << num_particles << '\n';

    int iseed = 451;
    double mass = 10.0;

    MyParticleContainer::ParticleInitData pdata = {{mass}, {}, {}, {}};
    myPC.InitRandom(num_particles, iseed, pdata, serialize);

    //
    // Here we provide an example of one way to call ParticleToMesh
    //
    amrex::ParticleToMesh(myPC, GetVecOfPtrs(density1), 0, parms.nlevs-1,
        [=] AMREX_GPU_DEVICE (const MyParticleContainer::ParticleType& p,
                              amrex::Array4<amrex::Real> const& rho,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi) noexcept
        {
            ParticleInterpolator::Linear interp(p, plo, dxi);

            interp.ParticleToMesh(p, rho, 0, 0, 1,
                [=] AMREX_GPU_DEVICE (const MyParticleContainer::ParticleType& part, int comp)
                {
                    return part.rdata(comp);  // no weighting
                });
        });

    //
    // Here we provide an example of another way to call ParticleToMesh
    //
    int start_part_comp = 0;
    int start_mesh_comp = 0;
    int        num_comp = 1;

    amrex::ParticleToMesh(myPC,GetVecOfPtrs(density2),0,parms.nlevs-1,
                          TrilinearDeposition{start_part_comp,start_mesh_comp,num_comp});

    //
    // Now write the output from each into separate plotfiles for comparison
    //

    Vector<std::string> varnames1, varnames2;
    // varname1 is for density1
    varnames1.push_back("density");  // has openPMD component scalar

    // varname2 is for density2, has two components, so assigned two names
    varnames2.push_back("velosity"); // has openPMD component scalar
    varnames2.push_back("Ex");       // has openPMD component E/x

    Vector<std::string> particle_varnames;
    particle_varnames.push_back("mass");

    int output_levs = parms.nlevs;

    Vector<const MultiFab*> outputMF1(output_levs);
    Vector<const MultiFab*> outputMF2(output_levs);
    Vector<IntVect> outputRR(output_levs);
    for (int lev = 0; lev < output_levs; ++lev) {
        outputMF1[lev] = &density1[lev];
        outputMF2[lev] = &density2[lev];
        outputRR[lev] = IntVect(AMREX_D_DECL(2, 2, 2));
    }

    //checkMFBox(parms, outputMF1);
    //checkMFBox(parms, outputMF2);

    // call count ptls to prepare ahead of time
    myPC.CountParticles();

    std::string fname;
    openpmd_api::InitHandler(fname);

    // one specie per particle container
    std::string specieName="ptlSpecie";

    int nsteps = 3;
    for (int ts = 0; ts < nsteps; ts++)
      {
        //Vector<int> level_steps;
        //level_steps.push_back(ts);
        //level_steps.push_back(ts);

        openpmd_api::SetStep(ts);
        openpmd_api::WriteParticles(myPC, specieName /*ts*/);  //  with default component names

        char name[512];
        std::snprintf(name, sizeof name, "plotfile_%05d",  ts);

        if ( 1 == parms.nlevs )
          {
          // example to store coarse level
          openpmd_api::WriteSingleLevel(*(outputMF1[0]), varnames1, geom[0], 0.0);
          openpmd_api::WriteSingleLevel(*(outputMF2[0]), varnames2, geom[0], 0.0);
          }
        else
          {
            // store multi mesh levels
            openpmd_api::WriteMultiLevel(//parms.nlevs,
                                         outputMF1,
                                         varnames1,
                                         geom,
                                         0.0,
                                         //level_steps,
                                         outputRR
                                         );
            // store multi mesh levels
            openpmd_api::WriteMultiLevel(//parms.nlevs,
                                         outputMF2,
                                         varnames2,
                                         geom,
                                         0.0,
                                         //level_steps,
                                         outputRR
                                         );
          }
        openpmd_api::CloseStep(ts);

        amrex::Print()<<"Timestep: "<<ts<<" done \n";
      }

    openpmd_api::CloseHandler();
}

void testBTD (TestParams& parms, int nghost)
{
    Vector<IntVect> rr(parms.nlevs-1);
    for (int lev = 1; lev < parms.nlevs; lev++) {
        rr[lev-1] = IntVect(AMREX_D_DECL(2,2,2));
    }

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1));
    const Box base_domain(domain_lo, domain_hi);

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM] = {AMREX_D_DECL(1,1,1)};

    Vector<Geometry> geom(parms.nlevs);
    geom[0].define(base_domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < parms.nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, CoordSys::cartesian, is_per);
    }

    Vector<BoxArray> ba(parms.nlevs);
    Vector<DistributionMapping> dm(parms.nlevs);

    Box domain = base_domain;
    IntVect size = IntVect(AMREX_D_DECL(parms.nx, parms.ny, parms.nz));
    for (int lev = 0; lev < parms.nlevs; ++lev)
    {
        ba[lev].define(domain);
        ba[lev].maxSize(parms.max_grid_size);
        dm[lev].define(ba[lev]);
        domain.grow(-size/4);   // fine level cover the middle of the coarse domain
        domain.refine(2);
    }

    Vector<MultiFab> density1(parms.nlevs);

    for (int lev = 0; lev < parms.nlevs; lev++) {
        // one field comp for density1
        density1[lev].define(ba[lev], dm[lev], 1, nghost);
        density1[lev].setVal(0.0);
    }

    MyParticleContainer myPC(geom, dm, ba, rr);
    myPC.SetVerbose(false);

    bool serialize = true;
    if (ParallelDescriptor::NProcs() > 1) {
        serialize = false;
    }

    amrex::Long num_particles = (amrex::Long)parms.nppc * parms.nx * parms.ny * parms.nz;
    amrex::Print() << serialize<< " Total number of particles   :" << num_particles << '\n';

    int iseed = 451;
    double mass = 10.0;

    MyParticleContainer::ParticleInitData pdata = {{mass}, {}, {}, {}};
    myPC.InitRandom(num_particles, iseed, pdata, serialize);

    //
    // Here we provide an example of one way to call ParticleToMesh
    //
    amrex::ParticleToMesh(myPC, GetVecOfPtrs(density1), 0, parms.nlevs-1,
        [=] AMREX_GPU_DEVICE (const MyParticleContainer::ParticleType& p,
                              amrex::Array4<amrex::Real> const& rho,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi) noexcept
        {
            ParticleInterpolator::Linear interp(p, plo, dxi);

            interp.ParticleToMesh(p, rho, 0, 0, 1,
                [=] AMREX_GPU_DEVICE (const MyParticleContainer::ParticleType& part, int comp)
                {
                    return part.rdata(comp);  // no weighting
                });
        });

    //
    // Now write the output from each into separate plotfiles for comparison
    //

    Vector<std::string> varnames1;
    // varname1 is for density1
    varnames1.push_back("density");  // has openPMD component scalar

    Vector<std::string> particle_varnames;
    particle_varnames.push_back("mass");

    int output_levs = parms.nlevs;

    Vector<const MultiFab*> outputMF1(output_levs);
    Vector<IntVect> outputRR(output_levs);
    for (int lev = 0; lev < output_levs; ++lev) {
        outputMF1[lev] = &density1[lev];
        outputRR[lev] = IntVect(AMREX_D_DECL(2, 2, 2));
    }

    //checkMFBox(parms, outputMF1);
    //checkMFBox(parms, outputMF2);

    // call count ptls to prepare ahead of time
    myPC.CountParticles();

    std::string fname;
    openpmd_api::InitHandler(fname);

    std::vector<bool> warpxPMLs(10);
    auto* testWriter = new AMReX_warpxBTDWriter(warpxPMLs); // xxxxx is there a memory leak?
    amrex::openpmd_api::UseCustomWriter(testWriter);

    // one specie per particle container
    std::string specieName="ptlSpecie";

    int nsteps = 4;
    //
    // mimicking  BTD behavior. Total is 2 actual steps, each steps is written twice
    // step writing order is 0 1 0 1, the last two writes are the final flushes
    // To make things simple, at the end of the second write, we should see double the ptls.
    // AsssignPtlOffsets(num_ptls) makes sure the second write starts off correctly
    //
    for (int its = 0; its < nsteps; its++)
      {
        int ts = its % 2;
        openpmd_api::SetStep(ts);

        if ( (its - 2) == ts )
          {
            //call AssignPtloffset() to assign the right starting offset of ptl batch
            testWriter->AssignPtlOffset(num_particles);
            testWriter->SetLastFlush();
          }
#if 0
        if (0)
          { //  test with default component names
            openpmd_api::WriteParticles(myPC, specieName);
          }
        else
#endif
          { // test with RZ style pos id
            openpmd_api::WriteParticles(myPC,
                                        specieName,
                                        [=] (auto& /*pc*/, openPMD::ParticleSpecies& currSpecies, unsigned long long localTotal)
                                        {
                                          amrex::ParticleReal lcharge = 0.01; // warpx: pc->getCharge()
                                          amrex::ParticleReal lmass = 0.5; // warpx: pc->getMass();

                                          testWriter->SetConstantMassCharge(currSpecies, localTotal, lcharge,  lmass);
                                        },
                                        [=] (auto& pti, openPMD::ParticleSpecies& currSpecies, unsigned long long offset)
                                        {
                                          testWriter->SavePosId_RZ(pti, currSpecies, offset); // also supports RZ
                                        });
          }

        {
          // store multi mesh levels
          openpmd_api::WriteMultiLevel(//parms.nlevs,
                                       outputMF1,
                                       varnames1,
                                       geom,
                                       0.0,
                                       outputRR
                                       );
        }

        openpmd_api::CloseStep(ts);

        amrex::Print()<<"Timestep: "<<ts<<" done \n";
      }

    openpmd_api::CloseHandler();
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
  pp.get("nppc", parms.nppc);
  pp.get("nlevs", parms.nlevs);
  //pp.get("nplotfile", parms.nplotfile);

  if (parms.nppc < 1 && ParallelDescriptor::IOProcessor()) {
    amrex::Abort("Must specify at least one particle per cell");
  }

  parms.verbose = false;
  pp.query("verbose", parms.verbose);

  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << std::endl;
    std::cout << "Number of particles per cell : ";
    std::cout << parms.nppc  << std::endl;
    std::cout << "Size of domain               : ";
    std::cout << parms.nx << " " << parms.ny << " " << parms.nz << std::endl;
  }


  int nghost = 1 ;
  std::cout<<"  TODO: RESOLVE!!!  if nghost=1  there is tile offset be  at -1  "<<std::endl;

#ifdef TEST_BTD
  testBTD(parms, nghost);
#else
  testParticleMesh(parms, nghost);
#endif

  amrex::Finalize();
}
