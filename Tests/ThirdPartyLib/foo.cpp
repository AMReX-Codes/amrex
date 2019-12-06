#include <mpi.h>
#include <AMReX.H>
#include <AMReX_Print.H>

void csolve ();

extern "C"
{
    void foo(MPI_Comm comm)
    {
        amrex::Initialize(comm);
        
        amrex::Print() << " foo: AMReX C++ has been initialized.\n";

        csolve();
        
        amrex::Finalize();

        // After amrex::Finalize(), amrex can no longer be used.
        int rank;
        MPI_Comm_rank(comm, &rank);
        if (rank == 0) {
            std::cout << " foo: AMReX C++ has been finalized.\n";
        }
    }
}

#include <AMReX_MLPoisson.H>
#include <AMReX_MLMG.H>

using namespace amrex;

void csolve ()
{
    RealBox rb({0.,0.,0.}, {1.,1.,1.}); // physical domain size
    std::array<int,3> is_periodic{1,1,1}; // periodic bc
    Geometry::Setup(&rb, 0, is_periodic.data());

    Box domain({0,0,0}, {63,63,63});  // # of cells

    Geometry geom(domain);

    BoxArray grids(domain);
    grids.maxSize(32);

    DistributionMapping dm(grids);

    MultiFab rhs(grids, dm, 1, 0);
    MultiFab phi(grids, dm, 1, 1);

    // set right hand side to some random numbers
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
    {
        auto const& fab = rhs.array(mfi);
        Box const& bx = mfi.fabbox();
        amrex::For(bx, [=] (int i, int j, int k) noexcept
        {
            fab(i,j,k) = Random();
        });
    }

    // set initial guess of potential to zero
    phi.setVal(0.0);

    MLPoisson mlpoisson({geom}, {grids}, {dm});

    mlpoisson.setDomainBC({LinOpBCType::Periodic,LinOpBCType::Periodic,LinOpBCType::Periodic},
                          {LinOpBCType::Periodic,LinOpBCType::Periodic,LinOpBCType::Periodic});
    mlpoisson.setLevelBC(0,nullptr);
   
    MLMG mlmg(mlpoisson);
    mlmg.setVerbose(1);
    mlmg.solve({&phi}, {&rhs}, 1.e-10, 0.0);
}
