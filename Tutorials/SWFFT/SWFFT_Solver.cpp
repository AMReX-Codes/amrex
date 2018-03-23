
#include <SWFFT_Solver.H>
#include <SWFFT_Solver_F.H>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

// These are for SWFFT
#include <Distribution.H>
#include <AlignedAllocator.h>
#include <Dfft.H>

#include <string>

#define ALIGN 16

using namespace amrex;

SWFFT_Solver::SWFFT_Solver ()
{
    static_assert(AMREX_SPACEDIM == 3, "3D only");

    // runtime parameters
    {
        ParmParse pp;

        pp.query("n_cell", n_cell);
        pp.query("max_grid_size", max_grid_size);

        pp.query("verbose", verbose);
    }
    
    BoxArray ba;
    {
        IntVect dom_lo(0,0,0);
        IntVect dom_hi(n_cell-1,n_cell-1,n_cell-1);
        Box domain(dom_lo,dom_hi);
        ba.define(domain);
        ba.maxSize(max_grid_size);

        // We assume a unit box (just for convenience)
        RealBox real_box({0.0,0.0,0.0}, {1.0,1.0,1.0});

        // The FFT assumes fully periodic boundaries
        std::array<int,3> is_periodic {1,1,1};

        geom.define(domain, &real_box, CoordSys::cartesian, is_periodic.data());
    }

    DistributionMapping dmap{ba};

    // Note that we are defining rhs with NO ghost cells
    rhs.define(ba, dmap, 1, 0);
    init_rhs();

    // Note that we are defining soln with NO ghost cells
    soln.define(ba, dmap, 1, 0);
    the_soln.define(ba, dmap, 1, 0);

    comp_the_solution();
}

void
SWFFT_Solver::init_rhs ()
{
    const Real* dx = geom.CellSize();

    for (MFIter mfi(rhs,true); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();
        fort_init_rhs(BL_TO_FORTRAN_BOX(tbx),
                      BL_TO_FORTRAN_ANYD(rhs[mfi]),
                      dx, &a, &b, &sigma, &w);
    }
}

void
SWFFT_Solver::comp_the_solution ()
{
    const Real* dx = geom.CellSize();

    for (MFIter mfi(the_soln); mfi.isValid(); ++mfi)
    {
        fort_comp_asol(BL_TO_FORTRAN_ANYD(the_soln[mfi]),dx);
    }
}

void
SWFFT_Solver::solve ()
{
    const BoxArray& ba = soln.boxArray();
    const DistributionMapping& dm = soln.DistributionMap();

    // Define (two pi) here
    Real  pi = 4 * std::atan(1.0);
    Real tpi = 2 * pi;

    // Assume one box for now
    int nx = ba[0].size()[0];
    int ny = ba[0].size()[1];
    int nz = ba[0].size()[2];

    // Assume for now that nx = ny = nz
    hacc::Distribution d(MPI_COMM_WORLD,nx);
    hacc::Dfft dfft(d);

    size_t local_size  = dfft.local_size();

    std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > a;
    std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > b;

    a.resize(nx*ny*nz);
    b.resize(nx*ny*nz);

    dfft.makePlans(&a[0],&b[0],&a[0],&b[0]);

    // Copy real data from Rhs into real part of a -- no ghost cells and
    // put into C++ ordering (not Fortran)
    complex_t zero(0.0, 0.0);
    complex_t one(1.0, 0.0);
    size_t local_indx = 0;
    for(int i=0; i< nx; i++) {
     for(int j=0; j< ny; j++) {
      for(int k=0; k< nz; k++) {

        complex_t temp(rhs[0].dataPtr()[local_indx],0.);
        a[local_indx] = temp;
	local_indx++;

      }
    }
   }

    VisMF::Write(rhs,"RHS_BEFORE");

//  *******************************************
//  Compute the forward transform
//  *******************************************
    dfft.forward(&a[0]);

//  *******************************************
//  Now divide the coefficients of the transform
//  *******************************************
    local_indx = 0;
    for(size_t i=0; i<(size_t)nx; i++) {
     for(size_t j=0; j<(size_t)ny; j++) {
      for(size_t k=0; k<(size_t)nz; k++) {

        if (i == 0 && j == 0 & k == 0) {
           a[local_indx] = 0;
        } else {
           double ksq = 2. * (
                        (cos(tpi*double(i)/double(nx)) - 1.) + 
                        (cos(tpi*double(j)/double(ny)) - 1.) + 
                        (cos(tpi*double(k)/double(nz)) - 1.) );

           a[local_indx] = a[local_indx] / ksq;
        }
	local_indx++;

      }
     }
    }

//  *******************************************
//  Compute the backward transform
//  *******************************************
    dfft.backward(&a[0]);

    double divisor = local_size * local_size / pi;

    local_indx = 0;
    for(size_t i=0; i<(size_t)nx; i++) {
     for(size_t j=0; j<(size_t)ny; j++) {
      for(size_t k=0; k<(size_t)nz; k++) {

        // Divide by 2 pi N
        soln[0].dataPtr()[local_indx] = -std::real(a[local_indx]) / divisor;
	local_indx++;

      }
     }
    }

    VisMF::Write(soln,"SOLN");
    VisMF::Write(the_soln,"THE_SOLN");

    std::cout << "MAX / MIN VALUE OF SOLN " <<     soln.max(0) << " " <<     soln.min(0) << std::endl;
    std::cout << "MAX / MIN VALUE OF THES " << the_soln.max(0) << " " << the_soln.min(0) << std::endl;
    std::cout << "RATIO OF MAX            " << soln.max(0) / the_soln.max(0) << std::endl;
    std::cout << "RATIO OF MIN            " << soln.min(0) / the_soln.min(0) << std::endl;

    {
        MultiFab diff(ba, dm, 1, 0);
        MultiFab::Copy(diff, soln, 0, 0, 1, 0);
        MultiFab::Subtract(diff, the_soln, 0, 0, 1, 0);
        amrex::Print() << "\nMax-norm of the error is " << diff.norm0()
                       << "\nL2-norm  of the error is " << diff.norm2() << "\n";
        amrex::Print() << "\nMaximum absolute value of the rhs is " << rhs.norm0()
                       << "\nMaximum absolute value of the solution is " << the_soln.norm0()
                       << "\n";
    }
}
