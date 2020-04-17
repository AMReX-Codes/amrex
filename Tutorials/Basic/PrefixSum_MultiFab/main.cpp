#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

using namespace amrex;

void main_main();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    main_main();

    amrex::Finalize();
    return 0;
}

void main_main ()
{
    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = amrex::second();

    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size;
    Vector<int> is_periodic(AMREX_SPACEDIM, 1);
    int reverse = 0;

    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the number of cells on each side of
        //   a square (or cubic) domain.
        pp.get("n_cell",n_cell);

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size",max_grid_size);

        // Should we reverse the prefix sum from hi to lo instead of lo to hi?
        // Default is to sum from lo to hi.
        pp.query("reverse",reverse);
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        ba.maxSize(max_grid_size);

       // This defines the physical box, [0,1] in each direction.
        RealBox real_box({AMREX_D_DECL(0.0,0.0,0.0)},
                         {AMREX_D_DECL(1.0,1.0,1.0)});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    const Real* dx = geom.CellSize();

    // Nghost = number of ghost cells for each array
    int Nghost = 0;

    // Ncomp = number of components for each array
    int Ncomp  = 1;

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // we allocate two multifabs; one will store the field phi, one will store the prefix sum of phi.
    MultiFab phi(ba, dm, Ncomp, Nghost);
    MultiFab phi_psum(ba, dm, Ncomp, Nghost);

    // init
    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);
        auto phi_fab = phi.array(mfi);

        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            phi_fab(i, j, k) = i;
        }}}
    }

    phi_psum = 0.0;
    Real time = 0.0;

    // Parallel reduce sum the partial sums on each grid

    // Write a plotfile of the prefix sum initial data
    {
        int n = 0;
        const std::string& pltfile = amrex::Concatenate("plt",n,5);
        WriteSingleLevelPlotfile(pltfile, phi_psum, {"phi_psum"}, geom, time, 0);
    }

    // prefix sum along the dimension 0 axis
    const int axis = 0;
    Long num_grids = ba.size();

    struct GidLo {
        Long gid;
        int lo;
        int reverse;
        Real exclusive_sum;

        GidLo(Long id, int smallend, int reverse_order) : gid(id), lo(smallend), reverse(reverse_order), exclusive_sum(0.0) {};

        // sort grids either increasing (forward) or decreasing (reverse) by their lo index
        bool operator<(const GidLo& other) const
        {
            if (reverse == 1)
                return lo > other.lo;
            else
                return lo < other.lo;
        }
    };

    // make a vector of GidLo structs with the grid ID, lo index, and exclusive sum for each grid in the domain
    Vector<GidLo> grids;

    for (Long i = 0; i < num_grids; ++i) {
        grids.emplace_back(i, ba[i].smallEnd(axis), reverse);
    }

    // sort the GidLo structs by their lo indexes
    std::sort(grids.begin(), grids.end());

    // make inverse map of grid IDs to indexes in the sorted grids
    Vector<Long> grid_inv(num_grids);
    for (Long i = 0; i < num_grids; ++i) grid_inv[grids[i].gid] = i;

    // array of prefix sums for each grid
    Vector<Real> prefix_sums(num_grids, 0.0);

    // compute the partial prefix sums for the grids on this MPI rank
    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        Long gidi = mfi.index();
        const Box& box = mfi.validbox();
        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);
        auto phi_fab = phi.array(mfi);
        auto phi_psum_fab = phi_psum.array(mfi);

        Real fab_prefix_sum = 0.0;
        if (reverse == 1) {
            for (int k = hi.z; k >= lo.z; --k) {
            for (int j = hi.y; j >= lo.y; --j) {
            for (int i = hi.x; i >= lo.x; --i) {
                fab_prefix_sum += phi_fab(i, j, k);
                phi_psum_fab(i, j, k) = fab_prefix_sum;
            }}}
        } else {
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                fab_prefix_sum += phi_fab(i, j, k);
                phi_psum_fab(i, j, k) = fab_prefix_sum;
            }}}
        }

        prefix_sums[gidi] = fab_prefix_sum;
    }

    // Parallel reduce sum the partial sums on each grid
    //   to communicate the per-grid partial sums across the domain.
    ParallelDescriptor::ReduceRealSum(prefix_sums.dataPtr(), num_grids);

    // The partial prefix sums are not yet accumulated between grids,
    //   so now we cumulative sum the partial sums in GidLo::exclusive_sum.
    //   This gets us the quantity we need to add to the partial sums of each grid.
    //   It is exclusive because it does not include the partial sum of its corresponding grid.
    for (Long i = 1; i < num_grids; ++i) {
        grids[i].exclusive_sum = grids[i-1].exclusive_sum + prefix_sums[grids[i-1].gid];
    }

    // Add the exclusive sums to the data on our corresponding local grids
    for (MFIter mfi(phi_psum); mfi.isValid(); ++mfi)
    {
        Long gidi = mfi.index();
        auto& phi_psum_fab = phi_psum[mfi];
        phi_psum_fab.plus<RunOn::Host>(grids[grid_inv[gidi]].exclusive_sum, 0, 1);
    }

    // Write a plotfile of the prefix sum data
    {
        int n = 1;
        const std::string& pltfile = amrex::Concatenate("plt",n,5);
        WriteSingleLevelPlotfile(pltfile, phi_psum, {"phi_psum"}, geom, time, 0);
    }

    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    Real stop_time = amrex::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;
}
