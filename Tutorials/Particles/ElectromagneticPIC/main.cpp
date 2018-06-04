#include <iostream>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

#include "ElectromagneticParticleContainer.H"
#include "Evolve.H"
#include "NodalFlags.H"
#include "Constants.H"
#include "IO.H"

#include "em_pic_F.H"

using namespace amrex;

struct TestParams
{
    IntVect ncell;      // num cells in domain
    IntVect nppc;       // number of particles per cell in each dim
    int max_grid_size;
    int nsteps;
    bool verbose;
};

void check_solution(const MultiFab& jx, const Geometry& geom, const Real time);

void test_em_pic(const TestParams& parms)
{
    BL_PROFILE("test_em_pic");
    BL_PROFILE_VAR_NS("evolve_time", blp_evolve);
    BL_PROFILE_VAR_NS("init_time", blp_init);
    
    BL_PROFILE_VAR_START(blp_init);
    
    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        real_box.setLo(n, -20e-6);
        real_box.setHi(n,  20e-6);
    }
    
    IntVect domain_lo(AMREX_D_DECL(0, 0, 0)); 
    IntVect domain_hi(AMREX_D_DECL(parms.ncell[0]-1,parms.ncell[1]-1,parms.ncell[2]-1));
    const Box domain(domain_lo, domain_hi);
    
    int coord = 0;
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) 
        is_per[i] = 1; 
    Geometry geom(domain, &real_box, coord, is_per);
    
    BoxArray ba(domain);
    ba.maxSize(parms.max_grid_size);
    DistributionMapping dm(ba);
    
    const int ng = 1;
    
    MultiFab Bx(amrex::convert(ba, YeeGrid::Bx_nodal_flag), dm, 1, ng);
    MultiFab By(amrex::convert(ba, YeeGrid::By_nodal_flag), dm, 1, ng);
    MultiFab Bz(amrex::convert(ba, YeeGrid::Bz_nodal_flag), dm, 1, ng);
    
    MultiFab Ex(amrex::convert(ba, YeeGrid::Ex_nodal_flag), dm, 1, ng);
    MultiFab Ey(amrex::convert(ba, YeeGrid::Ey_nodal_flag), dm, 1, ng);
    MultiFab Ez(amrex::convert(ba, YeeGrid::Ez_nodal_flag), dm, 1, ng);
    
    MultiFab jx(amrex::convert(ba, YeeGrid::jx_nodal_flag), dm, 1, ng);
    MultiFab jy(amrex::convert(ba, YeeGrid::jy_nodal_flag), dm, 1, ng);
    MultiFab jz(amrex::convert(ba, YeeGrid::jz_nodal_flag), dm, 1, ng);
    
    Bx.setVal(0.0); By.setVal(0.0); Bz.setVal(0.0);
    Ex.setVal(0.0); Ey.setVal(0.0); Ez.setVal(0.0);
    jx.setVal(0.0); jy.setVal(0.0); jz.setVal(0.0);
    
    std::cout << "Initializing particles... ";

    int num_species;
    Vector<ElectromagneticParticleContainer*> particles(2);
    int problem = 2;

    if (problem == 1) {
        num_species = 2;
    
        RealBox electron_bounds = RealBox(AMREX_D_DECL(-20e-6, -20e-6, -20e-6),
                                          AMREX_D_DECL( 20e-6, 20e-6, 20e-6));
        ElectromagneticParticleContainer* electrons;
        electrons = new ElectromagneticParticleContainer(geom, dm, ba, 
                                                   0, -PhysConst::q_e, PhysConst::m_e);
        electrons->InitParticles(parms.nppc, 0.01, 10.0, 1e25, real_box, problem);

        RealBox H_ions_bounds = RealBox(AMREX_D_DECL(-20e-6, -20e-6, -20e-6),
                                        AMREX_D_DECL( 20e-6,  20e-6,  20e-6));        
        ElectromagneticParticleContainer* H_ions;
        H_ions = new ElectromagneticParticleContainer(geom, dm, ba, 
                                                      1, PhysConst::q_e, PhysConst::m_p);
        H_ions->InitParticles(parms.nppc, 0.01, 10.0, 1e25, H_ions_bounds, problem);

        particles[0] = electrons;
        particles[1] = H_ions;
    }

    else if (problem == 2) {
        num_species = 1;

        RealBox electron_bounds = RealBox(AMREX_D_DECL(-20e-6, -20e-6, -20e-6),
                                          AMREX_D_DECL( 0.0,    20e-6,  20e-6));
        ElectromagneticParticleContainer* electrons;
        electrons = new ElectromagneticParticleContainer(geom, dm, ba,
                                                         0, -PhysConst::q_e, PhysConst::m_e);
        electrons->InitParticles(parms.nppc, 0.01, 10.0, 1e25, electron_bounds, problem);
        
        particles[0] = electrons;
    }

    std::cout << "Done. " << std::endl;
    
    for (int i = 0; i < num_species; ++i) particles[i]->OK();
    
    std::cout << "Starting main PIC loop... " << std::endl;
    
    int nsteps = parms.nsteps;
    const Real dt = compute_dt(geom);
    bool synchronized = true;
    
    BL_PROFILE_VAR_STOP(blp_init);
    
    BL_PROFILE_VAR_START(blp_evolve);
    
    Real time = 0.0;
    for (int step = 0; step < nsteps; ++step)
    {         
        std::cout << "    Time step: " <<  step << std::endl;
        
        if (synchronized)
        {
            for (int i = 0; i < num_species; ++i)
            {
                particles[i]->PushParticlesOnly(-0.5*dt);
            }
            synchronized = false;
        }
        else
        {
            fill_boundary_electric_field(Ex, Ey, Ez, geom);
            evolve_magnetic_field(Ex, Ey, Ez, Bx, By, Bz, geom, 0.5*dt);
            fill_boundary_magnetic_field(Bx, By, Bz, geom);
        }
        
        jx.setVal(0.0); jy.setVal(0.0), jz.setVal(0.0);

        amrex::Print() << jx.max(0) << std::endl;

        for (int i = 0; i < num_species; ++i) {
            particles[i]->PushAndDeposeParticles(Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, dt);
        }
        jx.SumBoundary(geom.periodicity());
        jy.SumBoundary(geom.periodicity());
        jz.SumBoundary(geom.periodicity());
        
        evolve_magnetic_field(Ex, Ey, Ez, Bx, By, Bz, geom, 0.5*dt);
        fill_boundary_magnetic_field(Bx, By, Bz, geom);

        evolve_electric_field(Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, geom, dt);
        
        if (step == nsteps - 1) {
            evolve_magnetic_field(Ex, Ey, Ez, Bx, By, Bz, geom, 0.5*dt);
            for (int i = 0; i < num_species; ++i) {
                particles[i]->PushParticlesOnly(0.5*dt);
            }
            synchronized = true;
        }
        
        for (int i = 0; i < num_species; ++i) {
            particles[i]->Redistribute();
            particles[i]->EnforcePeriodicBCs();
            particles[i]->OK();            
        }

        time += dt;
    }
    
    std::cout << "Done. " << std::endl;

    BL_PROFILE_VAR_STOP(blp_evolve);

    if (problem == 2) check_solution(jx, geom, time);

    WritePlotFile(Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, geom, time, nsteps);
}

void check_solution(const MultiFab& jx, const Geometry& geom, Real time)
{
    BL_PROFILE("ElectromagneticParticleContainer::check_solution");
    
    const Real* dx = geom.CellSize();

    Box test_box = geom.Domain();
    test_box.setSmall(IntVect(AMREX_D_DECL(2, 2, 2)));
    test_box.setBig(IntVect(AMREX_D_DECL(30, 30, 30)));

    Real max_error;
    for(MFIter mfi(jx); mfi.isValid(); ++mfi)
    {
        const Box& tbx  = mfi.tilebox();
        check_langmuir_solution(BL_TO_FORTRAN_BOX(tbx),
                                BL_TO_FORTRAN_BOX(test_box),
                                BL_TO_FORTRAN_3D(jx[mfi]), time, &max_error);
    }

    amrex::Print() << "Max error is: " << max_error << std::endl;
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    amrex::InitRandom(451);

    ParmParse pp;
    
    TestParams parms;
    
    pp.get("ncell", parms.ncell);
    pp.get("nppc",  parms.nppc);
    pp.get("max_grid_size", parms.max_grid_size);
    pp.get("nsteps", parms.nsteps);

    parms.verbose = false;
    pp.query("verbose", parms.verbose);
    
    test_em_pic(parms);
    
    amrex::Finalize();
}
