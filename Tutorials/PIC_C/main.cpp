#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <BLFort.H>
#include <MacBndry.H>
#include <MGT_Solver.H>
#include <mg_cpp_f.h>
#include <stencil_types.H>

#include "Particles.H"

// declare routines below
void solve_for_phi(MultiFab& rhs, MultiFab& grad_phi, const Geometry& geom);

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    // define the lower and upper corner of a 3D domain
    IntVect domain_lo(0 , 0, 0); 
    IntVect domain_hi(63,63,63); 
 
    // build a box for the domain
    Box domain(domain_lo, domain_hi);

    // build a box array from the 64^3 domain box
    BoxArray ba(domain);
    // break the box array into 32^3 boxes
    ba.maxSize(32);

    // build a multifab for the rhs on the box array with 1 component, 0 ghost cells
    MultiFab rhs(ba, 1, 0);  

    // build a multifab for grad_phi on the box array with BL_SPACEDIM components, 1 ghost cell
    MultiFab grad_phi(ba, BL_SPACEDIM, 1);  

    // Initialize phi to 0 everywhere.
    grad_phi.setVal(0.0);

    // This defines the physical size of the box.  Right now the box is [-1,1] in each direction.
    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
       real_box.setLo(n,0.0);
       real_box.setHi(n,1.0);
    }

    // This says we are using Cartesian coordinates
    int coord = 0;

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) is_per[i] = 1; 

    // This defines a Geometry object which is useful for writing the plotfiles  
    Geometry geom(domain,&real_box,coord,is_per);

    DistributionMapping dmap = rhs.DistributionMap();

    // Define a new particle container to hold my particles.
    // This holds a charge as well as three velocity components, three acceleration components  and three position components.
    typedef ParticleContainer<1+2*BL_SPACEDIM> MyParticleContainer;
    
    // Build a new particle container to hold my particles.
    MyParticleContainer* MyPC = new MyParticleContainer(geom,dmap,ba);

    int count = 1;
    int iseed = 10;
    Real mass  = 10.0;

    // Initialize "count" number of particles, each with mass/charge "mass"
    MyPC->InitRandom(count,iseed,mass);

    // Use the PIC approach to deposit the "mass" onto the grid
    MyPC->AssignDensitySingleLevel(rhs,0);

    // Write out the positions, masses and accelerations of each particle.
    MyPC->WriteAsciiFile("Particles_before");

    // Use multigrid to solve Lap(phi) = rhs with periodic boundary conditions (set above)
    solve_for_phi(rhs,grad_phi,geom);

    // Fill the particle data with the acceleration at the particle location
    int start_comp = BL_SPACEDIM+1;
    Real dummy_dt = 0.0;
    MyPC->moveKick(grad_phi,0,dummy_dt,1.0,1.0,start_comp);

    // Write out the positions, masses and accelerations of each particle.
    MyPC->WriteAsciiFile("Particles_after");

    BoxLib::Finalize();
}

void 
solve_for_phi(MultiFab& rhs, MultiFab& grad_phi, const Geometry& geom)
{
    Real     strt    = ParallelDescriptor::second();

    BCRec* phys_bc;
    int mg_bc[2*BL_SPACEDIM];

    int MAX_LEV = 10;
    Array< PArray<MultiFab> > grad_phi_edge;
    grad_phi_edge.resize(MAX_LEV);
    grad_phi_edge[0].resize(BL_SPACEDIM, PArrayManage);
    for (int n = 0; n < BL_SPACEDIM; ++n)
        grad_phi_edge[0].set(n, new MultiFab(BoxArray(rhs.boxArray()).surroundingNodes(n), 1, 1));

    // build (and initialize) a multifab for the solution on the box array with 1 component, 1 ghost cell
    MultiFab phi(BoxArray(rhs.boxArray()), 1, 1);
    phi.setVal(0.);

    std::cout << "Max of rhs in solve_for_phi before correction " << rhs.norm0() << std::endl;

    // This is a correction for fully periodic domains only
    if (Geometry::isAllPeriodic())
    {
        Real numpts = rhs.boxArray().d_numPts();
        Real sum = rhs.sum(0) / numpts;
        rhs.plus(-sum, 0, 1, 0);
    }

    std::cout << "Max of rhs in solve_for_phi  after correction " << rhs.norm0() << std::endl;

    // This tells the solver that we are using periodic bc's
    if (Geometry::isAllPeriodic())
    {
        for (int dir = 0; dir < BL_SPACEDIM; ++dir)
        {
            mg_bc[2*dir + 0] = 0;
            mg_bc[2*dir + 1] = 0;
        }
    }

    MacBndry bndry(rhs.boxArray(), 1, geom);
    //
    // Set Dirichlet boundary condition for phi in phi grow cells, use to
    // initialize bndry.
    //
    const int src_comp  = 0;
    const int dest_comp = 0;
    const int num_comp  = 1;

    // Need to set the boundary values here so they can get copied into "bndry"  
    {
        bndry.setBndryValues(phi, src_comp, dest_comp, num_comp, *phys_bc);
    }

    std::vector<BoxArray> bav(1);
    bav[0] = phi.boxArray();
    std::vector<DistributionMapping> dmv(1);
    dmv[0] = rhs.DistributionMap();
    std::vector<Geometry> fgeom(1);
    fgeom[0] = geom;

    Array< Array<Real> > xa(1);
    Array< Array<Real> > xb(1);

    xa[0].resize(BL_SPACEDIM);
    xb[0].resize(BL_SPACEDIM);

//  if (level == 0)
    {
        for (int i = 0; i < BL_SPACEDIM; ++i)
        {
            xa[0][i] = 0;
            xb[0][i] = 0;
        }
    }
#if 0
    else
    {
        const Real* dx_crse = parent->Geom(level-1).CellSize();
        for (int i = 0; i < BL_SPACEDIM; ++i)
        {
            xa[0][i] = 0.5 * dx_crse[i];
            xb[0][i] = 0.5 * dx_crse[i];
        }
    }
#endif

    MultiFab* phi_p[1] = {&phi};
    MultiFab* rhs_p[1] = {&rhs};

    const Real  tol     = 1.e-10;
    const Real  abs_tol = 0.;

#ifdef USEHPGMG
    if (solve_with_hpgmg)
    {
        solve_with_HPGMG(level, phi, grad_phi_edge, rhs, tol, abs_tol);
    }
    else
#endif
    {
        int stencil_type = CC_CROSS_STENCIL;
        int always_use_bnorm = 1;
        int need_grad_phi = 1;
        Real final_resnorm;
        MGT_Solver mgt_solver(fgeom, mg_bc, bav, dmv, false, stencil_type, false, 0, 1, 2);
        mgt_solver.set_const_gravity_coeffs(xa, xb);
        mgt_solver.solve(phi_p, rhs_p, bndry, tol, abs_tol, always_use_bnorm, final_resnorm,need_grad_phi); 

        const Real* dx = geom.CellSize();
        mgt_solver.get_fluxes(0,grad_phi_edge[0],dx);
    }

    // Average edge-centered gradients to cell centers.
    BoxLib::average_face_to_cellcenter(grad_phi, grad_phi_edge[0], geom);
    geom.FillPeriodicBoundary(grad_phi,true);

    VisMF::Write(grad_phi,"GradPhi");

    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt; 

#ifdef BL_LAZY 
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "solve_for_phi() time = " << end << std::endl;
#ifdef BL_LAZY 
        });
#endif
    }
}

