
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_KernelTimer.H>

using namespace amrex;

void main_main();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    main_main();
    amrex::Finalize();
}


void main_main ()
{
    BoxArray ba;
    {
        // This parameter can control GPU compute work
        int n_cell = 256;
        int max_grid_size = 64;
        ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.query("max_grid_size", max_grid_size);
        Box domain_box(IntVect(0), IntVect(n_cell-1));
        ba.define(domain_box);
        ba.maxSize(max_grid_size);
    }

    MultiFab mf(ba,DistributionMapping{ba},1,0);
    mf.setVal(0.0);

    // Sample KernelTimer instrumentation of amrex::ParallelFor function
    {
	// We pass this to KernelTimer to store accumulated thread cycles,
	// as a proxy for GPU compute work
	amrex::Real* cost = new amrex::Real(0.);
	
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
	        // GPU timer instrumentation; constructor takes two
		// arguments, first controls whether timer is active
		// and second is a pointer to the Real that stores 
		// accumulated GPU thread cycles
		amrex::KernelTimer KnlTimer(true, cost);
                fab(i,j,k) += 1.;
            });
        }
	// Now cost is filled with thread-wise summed cycles
	printf("I measured %.0f cycles over all threads.\n", *cost);
	delete cost;
    }

}
