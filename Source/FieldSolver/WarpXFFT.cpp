
#include <WarpX.H>
#include <WarpX_f.H>
#include <AMReX_iMultiFab.H>

using namespace amrex;

constexpr int WarpX::FFTData::N;

namespace {
static std::unique_ptr<WarpX::FFTData> nullfftdata; // This for process with nz_fft=0

/** \brief Returns an "owner mask" which 1 for all cells, except
 *  for the duplicated (physical) cells of a nodal grid.
 *
 *  More precisely, for these cells (which are represented on several grids)
 *  the owner mask is 1 only if these cells are at the lower left end of
 *  the local grid - or if these cells are at the end of the physical domain
 *  Therefore, there for these cells, there will be only one grid for
 *  which the owner mask is non-zero.
 */
static iMultiFab
BuildFFTOwnerMask (const MultiFab& mf, const Geometry& geom)
{
    const BoxArray& ba = mf.boxArray();
    const DistributionMapping& dm = mf.DistributionMap();
    iMultiFab mask(ba, dm, 1, 0);
    const int owner = 1;
    const int nonowner = 0;
    mask.setVal(owner);

    const Box& domain_box = amrex::convert(geom.Domain(), ba.ixType());

    AMREX_ASSERT(ba.complementIn(domain_box).isEmpty());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mask); mfi.isValid(); ++mfi)
    {
        IArrayBox& fab = mask[mfi];
        const Box& bx = fab.box();
        Box bx2 = bx;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            // Detect nodal dimensions
            if (bx2.type(idim) == IndexType::NODE) {
                // Make sure that this grid does not touch the end of
                // the physical domain.
                if (bx2.bigEnd(idim) < domain_box.bigEnd(idim)) {
                    bx2.growHi(idim, -1);
                }
            }
        }
        const BoxList& bl = amrex::boxDiff(bx, bx2);
        // Set owner mask in these cells
        for (const auto& b : bl) {
            fab.setVal(nonowner, b, 0, 1);
        }

    }

    return mask;
}

/** \brief Copy the data from the FFT grid to the regular grid
 *
 * Because, for nodal grid, some cells are duplicated on several boxes,
 * special care has to be taken in order to have consistent values on
 * each boxes when copying this data. Here this is done by setting a
 * mask, where, for these duplicated cells, the mask is non-zero on only
 * one box.
 */
static void
CopyDataFromFFTToValid (MultiFab& mf, const MultiFab& mf_fft, const BoxArray& ba_valid_fft, const Geometry& geom)
{
    auto idx_type = mf_fft.ixType();
    MultiFab mftmp(amrex::convert(ba_valid_fft,idx_type), mf_fft.DistributionMap(), 1, 0);

    const iMultiFab& mask = BuildFFTOwnerMask(mftmp, geom);

    // Local copy: whenever an MPI rank owns both the data from the FFT
    // grid and from the regular grid, for overlapping region, copy it locally
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mftmp,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        FArrayBox& dstfab = mftmp[mfi];

        const FArrayBox& srcfab = mf_fft[mfi];
        const Box& srcbox = srcfab.box();

        if (srcbox.contains(bx))
        {
            // Copy the interior region (without guard cells)
            dstfab.copy(srcfab, bx, 0, bx, 0, 1);
            // Set the value to 0 whenever the mask is 0
            // (i.e. for nodal duplicated cells, there is a single box
            // for which the mask is different than 0)
            // if mask == 0, set value to zero
            dstfab.setValIfNot(0.0, bx, mask[mfi], 0, 1);
        }
    }

    // Global copy: Get the remaining the data from other procs
    // Use ParallelAdd instead of ParallelCopy, so that the value from
    // the cell that has non-zero mask is the one which is retained.
    mf.setVal(0.0, 0);
    mf.ParallelAdd(mftmp);


}

}

void
WarpX::AllocLevelDataFFT (int lev)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev == 0, "PSATD doesn't work with mesh refinement yet");

    static_assert(std::is_standard_layout<FFTData>::value, "FFTData must have standard layout");
    static_assert(sizeof(FFTData) == sizeof(void*)*FFTData::N, "sizeof FFTData is wrong");

    InitFFTComm(lev);

    BoxArray ba_fp_fft;
    DistributionMapping dm_fp_fft;
    FFTDomainDecomposition(lev, ba_fp_fft, dm_fp_fft, ba_valid_fp_fft[lev], domain_fp_fft[lev],
                           geom[lev].Domain());

    // rho2 has one extra ghost cell, so that it's safe to deposit charge density after
    // pushing particle.

    Efield_fp_fft[lev][0].reset(new MultiFab(amrex::convert(ba_fp_fft,Ex_nodal_flag),
                                             dm_fp_fft, 1, 0));
    Efield_fp_fft[lev][1].reset(new MultiFab(amrex::convert(ba_fp_fft,Ey_nodal_flag),
                                             dm_fp_fft, 1, 0));
    Efield_fp_fft[lev][2].reset(new MultiFab(amrex::convert(ba_fp_fft,Ez_nodal_flag),
                                             dm_fp_fft, 1, 0));
    Bfield_fp_fft[lev][0].reset(new MultiFab(amrex::convert(ba_fp_fft,Bx_nodal_flag),
                                             dm_fp_fft, 1, 0));
    Bfield_fp_fft[lev][1].reset(new MultiFab(amrex::convert(ba_fp_fft,By_nodal_flag),
                                             dm_fp_fft, 1, 0));
    Bfield_fp_fft[lev][2].reset(new MultiFab(amrex::convert(ba_fp_fft,Bz_nodal_flag),
                                             dm_fp_fft, 1, 0));
    current_fp_fft[lev][0].reset(new MultiFab(amrex::convert(ba_fp_fft,jx_nodal_flag),
                                              dm_fp_fft, 1, 0));
    current_fp_fft[lev][1].reset(new MultiFab(amrex::convert(ba_fp_fft,jy_nodal_flag),
                                              dm_fp_fft, 1, 0));
    current_fp_fft[lev][2].reset(new MultiFab(amrex::convert(ba_fp_fft,jz_nodal_flag),
                                              dm_fp_fft, 1, 0));
    rho_fp_fft[lev].reset(new MultiFab(amrex::convert(ba_fp_fft,IntVect::TheNodeVector()),
                                       dm_fp_fft, 2, 0));

    dataptr_fp_fft[lev].reset(new LayoutData<FFTData>(ba_fp_fft, dm_fp_fft));

    if (lev > 0)
    {
        BoxArray ba_cp_fft;
        DistributionMapping dm_cp_fft;
        FFTDomainDecomposition(lev, ba_cp_fft, dm_cp_fft, ba_valid_cp_fft[lev], domain_cp_fft[lev],
                               amrex::coarsen(geom[lev].Domain(),2));

        Efield_cp_fft[lev][0].reset(new MultiFab(amrex::convert(ba_cp_fft,Ex_nodal_flag),
                                                 dm_cp_fft, 1, 0));
        Efield_cp_fft[lev][1].reset(new MultiFab(amrex::convert(ba_cp_fft,Ey_nodal_flag),
                                                 dm_cp_fft, 1, 0));
        Efield_cp_fft[lev][2].reset(new MultiFab(amrex::convert(ba_cp_fft,Ez_nodal_flag),
                                                 dm_cp_fft, 1, 0));
        Bfield_cp_fft[lev][0].reset(new MultiFab(amrex::convert(ba_cp_fft,Bx_nodal_flag),
                                                 dm_cp_fft, 1, 0));
        Bfield_cp_fft[lev][1].reset(new MultiFab(amrex::convert(ba_cp_fft,By_nodal_flag),
                                                 dm_cp_fft, 1, 0));
        Bfield_cp_fft[lev][2].reset(new MultiFab(amrex::convert(ba_cp_fft,Bz_nodal_flag),
                                                 dm_cp_fft, 1, 0));
        current_cp_fft[lev][0].reset(new MultiFab(amrex::convert(ba_cp_fft,jx_nodal_flag),
                                                  dm_cp_fft, 1, 0));
        current_cp_fft[lev][1].reset(new MultiFab(amrex::convert(ba_cp_fft,jy_nodal_flag),
                                                  dm_cp_fft, 1, 0));
        current_cp_fft[lev][2].reset(new MultiFab(amrex::convert(ba_cp_fft,jz_nodal_flag),
                                                  dm_cp_fft, 1, 0));
        rho_cp_fft[lev].reset(new MultiFab(amrex::convert(ba_cp_fft,IntVect::TheNodeVector()),
                                           dm_cp_fft, 2, 0));

        dataptr_cp_fft[lev].reset(new LayoutData<FFTData>(ba_cp_fft, dm_cp_fft));
    }

    InitFFTDataPlan(lev);
}

/** \brief Create MPI sub-communicators for each FFT group,
 *         and put them in PICSAR module
 *
 * These communicators are passed to the parallel FFTW library, in order
 * to perform a global FFT within each FFT group.
 */
void
WarpX::InitFFTComm (int lev)
{
    int nprocs = ParallelDescriptor::NProcs();
    ngroups_fft = std::min(ngroups_fft, nprocs);

    // # of processes in the subcommunicator
    int np_fft = nprocs / ngroups_fft;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(np_fft*ngroups_fft == nprocs,
        "Number of processes must be divisible by number of FFT groups");

    int myproc = ParallelDescriptor::MyProc();
    // my color in ngroups_fft subcommunicators.  0 <= color_fft < ngroups_fft
    color_fft[lev] = myproc / np_fft;
    MPI_Comm_split(ParallelDescriptor::Communicator(), color_fft[lev], myproc, &comm_fft[lev]);

    int fcomm = MPI_Comm_c2f(comm_fft[lev]);
    // Set the communicator of the PICSAR module to the one we just created
    warpx_fft_mpi_init(fcomm);
}

/** \brief Perform domain decomposition for the FFTW
 *
 *  Attribute one (unique) box to each proc, in such a way that:
 *    - The global domain is divided among FFT groups,
 *      with additional guard cells around each FFT group
 *    - The domain associated to an FFT group (with its guard cells)
 *      is further divided in sub-subdomains along z, so as to distribute
 *      it among the procs within an FFT group
 *
 *  The attribution is done by setting (within this function):
 *  - ba_fft: the BoxArray representing the final set of sub-domains for the FFT
 *            (includes/covers the guard cells of the FFT groups)
 *  - dm_fft: the mapping between these sub-domains and the corresponding proc
 *              (imposes one unique box for each proc)
 *  - ba_valid: the BoxArray that contains valid part of the sub-domains of ba_fft
 *            (i.e. does not include/cover the guard cells of the FFT groups)
 *  - domain_fft: a Box that represent the domain of the FFT group for the current proc
 */
void
WarpX::FFTDomainDecomposition (int lev, BoxArray& ba_fft, DistributionMapping& dm_fft,
                               BoxArray& ba_valid, Box& domain_fft, const Box& domain)
{

    IntVect nguards_fft(AMREX_D_DECL(nox_fft/2,noy_fft/2,noz_fft/2));

    int nprocs = ParallelDescriptor::NProcs();

    BoxList bl(domain, ngroups_fft);  // This does a multi-D domain decomposition for groups
    AMREX_ALWAYS_ASSERT(bl.size() == ngroups_fft);
    const Vector<Box>& bldata = bl.data();

    // This is the domain for the FFT sub-group (including guard cells)
    domain_fft = amrex::grow(bldata[color_fft[lev]], nguards_fft);
    // Ask FFTW to chop the current FFT sub-group domain in the z-direction
    // and give a chunk to each MPI rank in the current sub-group.
    int nz_fft, z0_fft;

    warpx_fft_domain_decomp(&nz_fft, &z0_fft, WARPX_TO_FORTRAN_BOX(domain_fft));
    // Each MPI rank adds a box with its chunk of the FFT grid
    // (given by the above decomposition) to the list `bx_fft`,
    // then list is shared among all MPI ranks via AllGather
    Vector<Box> bx_fft;
    if (nz_fft > 0) {
        Box b = domain_fft;
        b.setRange(AMREX_SPACEDIM-1, z0_fft+domain_fft.smallEnd(AMREX_SPACEDIM-1), nz_fft);
        bx_fft.push_back(b);
    } else {
        // Add empty box for the AllGather call
        bx_fft.push_back(Box());
    }
    amrex::AllGatherBoxes(bx_fft);
    AMREX_ASSERT(bx_fft.size() == ParallelDescriptor::NProcs());
    // Build pmap and bx_fft without the empty boxes
    Vector<int> pmap;
    for (int i = 0; i < bx_fft.size(); ++i) {
        if (bx_fft[i].ok()) {
            pmap.push_back(i);
        }
    }
    bx_fft.erase(std::remove_if(bx_fft.begin(),bx_fft.end(),
                                [](Box const& b) { return b.isEmpty(); }),
                 bx_fft.end());
    AMREX_ASSERT(bx_fft.size() == pmap.size());

    // Define the AMReX objects for the FFT grid: BoxArray and DistributionMapping
    ba_fft.define(BoxList(std::move(bx_fft)));
    dm_fft.define(std::move(pmap));

    // For communication between WarpX normal domain and FFT domain, we need to create a
    // special BoxArray ba_valid
    const Box foobox(-nguards_fft-2, -nguards_fft-2);

    BoxList bl_valid; // List of boxes: will be filled by the valid part of the subdomains of ba_fft
    bl_valid.reserve(ba_fft.size());
    int np_fft = nprocs / ngroups_fft;
    for (int i = 0; i < ba_fft.size(); ++i)
    {
        int igroup = dm_fft[i] / np_fft; // This should be consistent with InitFFTComm
        const Box& bx = ba_fft[i] & bldata[igroup]; // Intersection with the domain of
                                                    // the FFT group *without* guard cells
        if (bx.ok())
        {
            bl_valid.push_back(bx);
        }
        else
        {
            bl_valid.push_back(foobox);
        }
    }

    ba_valid.define(std::move(bl_valid));
}

/** /brief Set all the flags and metadata of the PICSAR FFT module.
 *         Allocate the auxiliary arrays of `fft_data`
 *
 * Note: dataptr_data is a stuct containing 22 pointers to arrays
 * 1-11: padded arrays in real space ; 12-22 arrays for the fields in Fourier space
 */
void
WarpX::InitFFTDataPlan (int lev)
{
    auto dx_fp = CellSize(lev);

    if (Efield_fp_fft[lev][0]->local_size() == 1)
      //Only one FFT patch on this MPI
    {
        for (MFIter mfi(*Efield_fp_fft[lev][0]); mfi.isValid(); ++mfi)
        {
            warpx_fft_dataplan_init(&nox_fft, &noy_fft, &noz_fft,
                                    (*dataptr_fp_fft[lev])[mfi].data, &FFTData::N,
                                    dx_fp.data(), &dt[lev], &fftw_plan_measure, &WarpX::do_nodal );
        }
    }
    else if (Efield_fp_fft[lev][0]->local_size() == 0)
      // No FFT patch on this MPI rank (may happen with FFTW)
      // Still need to call the MPI-FFT initialization routines
    {
	nullfftdata.reset(new FFTData());
        warpx_fft_dataplan_init(&nox_fft, &noy_fft, &noz_fft,
                                nullfftdata->data, &FFTData::N,
                                dx_fp.data(), &dt[lev], &fftw_plan_measure,
                                &WarpX::do_nodal );
    }
    else
    {
      // Multiple FFT patches on this MPI rank
        amrex::Abort("WarpX::InitFFTDataPlan: TODO");
    }

    if (lev > 0)
    {
        amrex::Abort("WarpX::InitFFTDataPlan: TODO");
    }
}

void
WarpX::FreeFFT (int lev)
{
    nullfftdata.reset();

    warpx_fft_nullify();

    if (comm_fft[lev] != MPI_COMM_NULL) {
        MPI_Comm_free(&comm_fft[lev]);
    }
    comm_fft[lev] = MPI_COMM_NULL;
}

void
WarpX::PushPSATD (amrex::Real a_dt)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt[lev] == a_dt, "dt must be consistent");
        if (fft_hybrid_mpi_decomposition){
            PushPSATD_hybridFFT(lev, a_dt);
        } else {
            PushPSATD_localFFT(lev, a_dt);
        }
    }
}

void WarpX::PushPSATD_localFFT (int lev, amrex::Real /* dt */)
{
    auto& solver = *spectral_solver_fp[lev];

    // Perform forward Fourier transform
    solver.ForwardTransform(*Efield_fp[lev][0], SpectralFieldIndex::Ex);
    solver.ForwardTransform(*Efield_fp[lev][1], SpectralFieldIndex::Ey);
    solver.ForwardTransform(*Efield_fp[lev][2], SpectralFieldIndex::Ez);
    solver.ForwardTransform(*Bfield_fp[lev][0], SpectralFieldIndex::Bx);
    solver.ForwardTransform(*Bfield_fp[lev][1], SpectralFieldIndex::By);
    solver.ForwardTransform(*Bfield_fp[lev][2], SpectralFieldIndex::Bz);
    solver.ForwardTransform(*current_fp[lev][0], SpectralFieldIndex::Jx);
    solver.ForwardTransform(*current_fp[lev][1], SpectralFieldIndex::Jy);
    solver.ForwardTransform(*current_fp[lev][2], SpectralFieldIndex::Jz);
    solver.ForwardTransform(*rho_fp[lev], SpectralFieldIndex::rho_old, 0);
    solver.ForwardTransform(*rho_fp[lev], SpectralFieldIndex::rho_new, 1);

    // Advance fields in spectral space
    solver.pushSpectralFields();

    // Perform backward Fourier Transform
    solver.BackwardTransform(*Efield_fp[lev][0], SpectralFieldIndex::Ex);
    solver.BackwardTransform(*Efield_fp[lev][1], SpectralFieldIndex::Ey);
    solver.BackwardTransform(*Efield_fp[lev][2], SpectralFieldIndex::Ez);
    solver.BackwardTransform(*Bfield_fp[lev][0], SpectralFieldIndex::Bx);
    solver.BackwardTransform(*Bfield_fp[lev][1], SpectralFieldIndex::By);
    solver.BackwardTransform(*Bfield_fp[lev][2], SpectralFieldIndex::Bz);
}

void
WarpX::PushPSATD_hybridFFT (int lev, amrex::Real /* dt */)
{
#ifndef AMREX_USE_CUDA // Running on CPU ; use PICSAR code for the hybrid FFT

    BL_PROFILE_VAR_NS("WarpXFFT::CopyDualGrid", blp_copy);
    BL_PROFILE_VAR_NS("PICSAR::FftPushEB", blp_push_eb);

    auto period_fp = geom[lev].periodicity();

    BL_PROFILE_VAR_START(blp_copy);
    Efield_fp_fft[lev][0]->ParallelCopy(*Efield_fp[lev][0], 0, 0, 1, Efield_fp[lev][0]->nGrow(), 0, period_fp);
    Efield_fp_fft[lev][1]->ParallelCopy(*Efield_fp[lev][1], 0, 0, 1, Efield_fp[lev][1]->nGrow(), 0, period_fp);
    Efield_fp_fft[lev][2]->ParallelCopy(*Efield_fp[lev][2], 0, 0, 1, Efield_fp[lev][2]->nGrow(), 0, period_fp);
    Bfield_fp_fft[lev][0]->ParallelCopy(*Bfield_fp[lev][0], 0, 0, 1, Bfield_fp[lev][0]->nGrow(), 0, period_fp);
    Bfield_fp_fft[lev][1]->ParallelCopy(*Bfield_fp[lev][1], 0, 0, 1, Bfield_fp[lev][1]->nGrow(), 0, period_fp);
    Bfield_fp_fft[lev][2]->ParallelCopy(*Bfield_fp[lev][2], 0, 0, 1, Bfield_fp[lev][2]->nGrow(), 0, period_fp);
    current_fp_fft[lev][0]->ParallelCopy(*current_fp[lev][0], 0, 0, 1, current_fp[lev][0]->nGrow(), 0, period_fp);
    current_fp_fft[lev][1]->ParallelCopy(*current_fp[lev][1], 0, 0, 1, current_fp[lev][1]->nGrow(), 0, period_fp);
    current_fp_fft[lev][2]->ParallelCopy(*current_fp[lev][2], 0, 0, 1, current_fp[lev][2]->nGrow(), 0, period_fp);
    rho_fp_fft[lev]->ParallelCopy(*rho_fp[lev], 0, 0, 2, rho_fp[lev]->nGrow(), 0, period_fp);
    BL_PROFILE_VAR_STOP(blp_copy);

    BL_PROFILE_VAR_START(blp_push_eb);
    if (Efield_fp_fft[lev][0]->local_size() == 1)
       //Only one FFT patch on this MPI
    {
    	for (MFIter mfi(*Efield_fp_fft[lev][0]); mfi.isValid(); ++mfi)
    	{
                warpx_fft_push_eb(WARPX_TO_FORTRAN_ANYD((*Efield_fp_fft[lev][0])[mfi]),
                                  WARPX_TO_FORTRAN_ANYD((*Efield_fp_fft[lev][1])[mfi]),
                                  WARPX_TO_FORTRAN_ANYD((*Efield_fp_fft[lev][2])[mfi]),
                                  WARPX_TO_FORTRAN_ANYD((*Bfield_fp_fft[lev][0])[mfi]),
                                  WARPX_TO_FORTRAN_ANYD((*Bfield_fp_fft[lev][1])[mfi]),
                                  WARPX_TO_FORTRAN_ANYD((*Bfield_fp_fft[lev][2])[mfi]),
                                  WARPX_TO_FORTRAN_ANYD((*current_fp_fft[lev][0])[mfi]),
                                  WARPX_TO_FORTRAN_ANYD((*current_fp_fft[lev][1])[mfi]),
                                  WARPX_TO_FORTRAN_ANYD((*current_fp_fft[lev][2])[mfi]),
                                  WARPX_TO_FORTRAN_N_ANYD((*rho_fp_fft[lev])[mfi],0),
                                  WARPX_TO_FORTRAN_N_ANYD((*rho_fp_fft[lev])[mfi],1));
    	}
    }
    else if (Efield_fp_fft[lev][0]->local_size() == 0)
      // No FFT patch on this MPI rank
      // Still need to call the MPI-FFT routine.
    {
    	FArrayBox fab(Box(IntVect::TheZeroVector(), IntVect::TheUnitVector()));
    	warpx_fft_push_eb(WARPX_TO_FORTRAN_ANYD(fab),
    			  WARPX_TO_FORTRAN_ANYD(fab),
    			  WARPX_TO_FORTRAN_ANYD(fab),
    			  WARPX_TO_FORTRAN_ANYD(fab),
    			  WARPX_TO_FORTRAN_ANYD(fab),
    			  WARPX_TO_FORTRAN_ANYD(fab),
    			  WARPX_TO_FORTRAN_ANYD(fab),
    			  WARPX_TO_FORTRAN_ANYD(fab),
    			  WARPX_TO_FORTRAN_ANYD(fab),
    			  WARPX_TO_FORTRAN_ANYD(fab),
    			  WARPX_TO_FORTRAN_ANYD(fab));
    }
    else
      // Multiple FFT patches on this MPI rank
    {
	amrex::Abort("WarpX::PushPSATD: TODO");
    }
    BL_PROFILE_VAR_STOP(blp_push_eb);

    BL_PROFILE_VAR_START(blp_copy);
    CopyDataFromFFTToValid(*Efield_fp[lev][0], *Efield_fp_fft[lev][0], ba_valid_fp_fft[lev], geom[lev]);
    CopyDataFromFFTToValid(*Efield_fp[lev][1], *Efield_fp_fft[lev][1], ba_valid_fp_fft[lev], geom[lev]);
    CopyDataFromFFTToValid(*Efield_fp[lev][2], *Efield_fp_fft[lev][2], ba_valid_fp_fft[lev], geom[lev]);
    CopyDataFromFFTToValid(*Bfield_fp[lev][0], *Bfield_fp_fft[lev][0], ba_valid_fp_fft[lev], geom[lev]);
    CopyDataFromFFTToValid(*Bfield_fp[lev][1], *Bfield_fp_fft[lev][1], ba_valid_fp_fft[lev], geom[lev]);
    CopyDataFromFFTToValid(*Bfield_fp[lev][2], *Bfield_fp_fft[lev][2], ba_valid_fp_fft[lev], geom[lev]);
    BL_PROFILE_VAR_STOP(blp_copy);

    if (lev > 0)
    {
        amrex::Abort("WarpX::PushPSATD: TODO");
    }
#else // AMREX_USE_CUDA is defined ; running on GPU
        amrex::Abort("The option `psatd.fft_hybrid_mpi_decomposition` does not work on GPU.");
#endif

}
