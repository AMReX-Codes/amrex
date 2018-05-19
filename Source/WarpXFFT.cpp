
#include <WarpX.H>
#include <WarpX_f.H>
#include <AMReX_BaseFab_f.H>
#include <AMReX_iMultiFab.H>

using namespace amrex;

constexpr int WarpX::FFTData::N;

namespace {

static iMultiFab
BuildFFTOwnerMask (const MultiFab& mf)
{
    const BoxArray& ba = mf.boxArray();
    const DistributionMapping& dm = mf.DistributionMap();
    iMultiFab mask(ba, dm, 1, 0);
    const int owner = 1;
    const int nonowner = 0;
    mask.setVal(owner);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector< std::pair<int,Box> > isects;
        for (MFIter mfi(mask); mfi.isValid(); ++mfi)
        {
            IArrayBox& fab = mask[mfi];
            const Box& bx = fab.box();
            const int igrid = mfi.index();
            ba.intersections(bx, isects);
            for (const auto& is : isects)
            {
                if (is.first != igrid)
                {
                    const Box& ibx = is.second;
                    bool on_high_end = false;
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                        Box bhi = bx;
                        bhi.setSmall(idim, bx.bigEnd(idim));
                        if (bhi.contains(ibx)) {
                            on_high_end = true;
                            break;
                        }
                    }
                    if (on_high_end) {
                        fab.setVal(nonowner, ibx, 0, 1);
                    }
                }
            }
        }
    }

    return mask;
}

static void
CopyDataFromFFTToValid (MultiFab& mf, const MultiFab& mf_fft, const BoxArray& ba_valid_fft)
{
    auto idx_type = mf_fft.ixType();
    MultiFab mftmp(amrex::convert(ba_valid_fft,idx_type), mf_fft.DistributionMap(), 1, 0);

    const iMultiFab& mask = BuildFFTOwnerMask(mftmp);

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
            dstfab.copy(srcfab, bx, 0, bx, 0, 1);
            amrex_fab_setval_ifnot (BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_FAB(dstfab),
                                    BL_TO_FORTRAN_ANYD(mask[mfi]),
                                    0.0); // if mask == 0, set value to zero
        }
    }

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
    FFTDomainDecompsition(lev, ba_fp_fft, dm_fp_fft, ba_valid_fp_fft[lev], domain_fp_fft[lev],
                          geom[lev].Domain());

    int ngRho = Efield_fp[lev][0]->nGrow();
    if (rho_fp[lev] == nullptr)
    {
        const BoxArray& ba = Efield_fp[lev][0]->boxArray();
        const DistributionMapping& dm = Efield_fp[lev][0]->DistributionMap();
        rho_fp[lev].reset(new MultiFab(amrex::convert(ba,IntVect::TheUnitVector()),dm,1,ngRho));
    }

    rho2_fp[lev].reset(new MultiFab(rho_fp[lev]->boxArray(),
                                    rho_fp[lev]->DistributionMap(),
                                    1, ngRho+1));
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
    rho_prev_fp_fft[lev].reset(new MultiFab(amrex::convert(ba_fp_fft,IntVect::TheNodeVector()),
                                            dm_fp_fft, 1, 0));
    rho_next_fp_fft[lev].reset(new MultiFab(amrex::convert(ba_fp_fft,IntVect::TheNodeVector()),
                                            dm_fp_fft, 1, 0));

    dataptr_fp_fft[lev].reset(new LayoutData<FFTData>(ba_fp_fft, dm_fp_fft));

    if (lev > 0)
    {
        BoxArray ba_cp_fft;
        DistributionMapping dm_cp_fft;
        FFTDomainDecompsition(lev, ba_cp_fft, dm_cp_fft, ba_valid_cp_fft[lev], domain_cp_fft[lev],
                              amrex::coarsen(geom[lev].Domain(),2));

        int ngRho = Efield_cp[lev][0]->nGrow();
        if (rho_cp[lev] == nullptr)
        {
            const BoxArray& ba = Efield_cp[lev][0]->boxArray();
            const DistributionMapping& dm = Efield_cp[lev][0]->DistributionMap();
            rho_cp[lev].reset(new MultiFab(amrex::convert(ba,IntVect::TheUnitVector()),dm,1,ngRho));
        }

        rho2_cp[lev].reset(new MultiFab(rho_cp[lev]->boxArray(),
                                        rho_cp[lev]->DistributionMap(),
                                        1, ngRho));
        // rho2 has one extra ghost cell, so that it's safe to deposit charge density after
        // pushing particle.

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
        rho_prev_cp_fft[lev].reset(new MultiFab(amrex::convert(ba_cp_fft,IntVect::TheNodeVector()),
                                                dm_cp_fft, 1, 0));
        rho_next_cp_fft[lev].reset(new MultiFab(amrex::convert(ba_cp_fft,IntVect::TheNodeVector()),
                                                dm_cp_fft, 1, 0));

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
                                     "Number of processes must be divisilbe by number of FFT groups");

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
WarpX::FFTDomainDecompsition (int lev, BoxArray& ba_fft, DistributionMapping& dm_fft,
                              BoxArray& ba_valid, Box& domain_fft, const Box& domain)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMREX_SPACEDIM == 3, "PSATD only works in 3D");

    IntVect nguards_fft(AMREX_D_DECL(nox_fft/2,noy_fft/2,noz_fft/2));

    int nprocs = ParallelDescriptor::NProcs();

    BoxList bl(domain, ngroups_fft);  // This does a multi-D domain decomposition for groups
    AMREX_ALWAYS_ASSERT(bl.size() == ngroups_fft);
    const Vector<Box>& bldata = bl.data();

    // This is the domain for the group.
    domain_fft = amrex::grow(bldata[color_fft[lev]], nguards_fft);
    // Ask FFTW to chop in z-direction into pieces
    int nz_fft, z0_fft;
    warpx_fft_domain_decomp(&nz_fft, &z0_fft, BL_TO_FORTRAN_BOX(domain_fft));

    Vector<Box> bx_fft;
    if (nz_fft > 0) {
        Box b = domain_fft;
        b.setRange(2, z0_fft+domain_fft.smallEnd(2), nz_fft);
        bx_fft.push_back(b);
    }
    amrex::AllGatherBoxes(bx_fft);

    ba_fft.define(BoxList(std::move(bx_fft)));
    AMREX_ALWAYS_ASSERT(ba_fft.size() == ParallelDescriptor::NProcs());

    Vector<int> pmap(ba_fft.size());
    std::iota(pmap.begin(), pmap.end(), 0);
    dm_fft.define(std::move(pmap));

    //
    // For communication between WarpX normal domain and FFT domain, we need to create a
    // special BoxArray ba_valid
    //

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
    AMREX_ALWAYS_ASSERT(Efield_fp_fft[lev][0]->local_size() == 1);

    auto dx_fp = CellSize(lev);

    for (MFIter mfi(*Efield_fp_fft[lev][0]); mfi.isValid(); ++mfi)
    {
        const Box& local_domain = amrex::enclosedCells(mfi.fabbox());
        warpx_fft_dataplan_init(&nox_fft, &noy_fft, &noz_fft,
                                (*dataptr_fp_fft[lev])[mfi].data, &FFTData::N,
                                dx_fp.data(), &dt[lev], &fftw_plan_measure );
    }

    if (lev > 0)
    {
        amrex::Abort("WarpX::InitFFTDataPlan: TODO");
    }
}

void
WarpX::FreeFFT (int lev)
{
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
        PushPSATD(lev, a_dt);
    }
}

void
WarpX::PushPSATD (int lev, amrex::Real /* dt */)
{
    BL_PROFILE_VAR_NS("WarpXFFT::CopyDualGrid", blp_copy);
    BL_PROFILE_VAR_NS("PICSAR::FftPushEB", blp_push_eb);

    auto period_fp = geom[lev].periodicity();

    BL_PROFILE_VAR_START(blp_copy);
    Efield_fp_fft[lev][0]->ParallelCopy(*Efield_fp[lev][0], 0, 0, 1, 0, 0, period_fp);
    Efield_fp_fft[lev][1]->ParallelCopy(*Efield_fp[lev][1], 0, 0, 1, 0, 0, period_fp);
    Efield_fp_fft[lev][2]->ParallelCopy(*Efield_fp[lev][2], 0, 0, 1, 0, 0, period_fp);
    Bfield_fp_fft[lev][0]->ParallelCopy(*Bfield_fp[lev][0], 0, 0, 1, 0, 0, period_fp);
    Bfield_fp_fft[lev][1]->ParallelCopy(*Bfield_fp[lev][1], 0, 0, 1, 0, 0, period_fp);
    Bfield_fp_fft[lev][2]->ParallelCopy(*Bfield_fp[lev][2], 0, 0, 1, 0, 0, period_fp);
    current_fp_fft[lev][0]->ParallelCopy(*current_fp[lev][0], 0, 0, 1, 0, 0, period_fp);
    current_fp_fft[lev][1]->ParallelCopy(*current_fp[lev][1], 0, 0, 1, 0, 0, period_fp);
    current_fp_fft[lev][2]->ParallelCopy(*current_fp[lev][2], 0, 0, 1, 0, 0, period_fp);
    rho_prev_fp_fft[lev]->ParallelCopy(*rho_fp[lev], 0, 0, 1, 0, 0, period_fp);
    rho_next_fp_fft[lev]->ParallelCopy(*rho2_fp[lev], 0, 0, 1, 0, 0, period_fp);
    BL_PROFILE_VAR_STOP(blp_copy);

    BL_PROFILE_VAR_START(blp_push_eb);
    for (MFIter mfi(*Efield_fp_fft[lev][0]); mfi.isValid(); ++mfi)
    {
        warpx_fft_push_eb(BL_TO_FORTRAN_ANYD((*Efield_fp_fft[lev][0])[mfi]),
                          BL_TO_FORTRAN_ANYD((*Efield_fp_fft[lev][1])[mfi]),
                          BL_TO_FORTRAN_ANYD((*Efield_fp_fft[lev][2])[mfi]),
                          BL_TO_FORTRAN_ANYD((*Bfield_fp_fft[lev][0])[mfi]),
                          BL_TO_FORTRAN_ANYD((*Bfield_fp_fft[lev][1])[mfi]),
                          BL_TO_FORTRAN_ANYD((*Bfield_fp_fft[lev][2])[mfi]),
                          BL_TO_FORTRAN_ANYD((*current_fp_fft[lev][0])[mfi]),
                          BL_TO_FORTRAN_ANYD((*current_fp_fft[lev][1])[mfi]),
                          BL_TO_FORTRAN_ANYD((*current_fp_fft[lev][2])[mfi]),
                          BL_TO_FORTRAN_ANYD((*rho_prev_fp_fft[lev])[mfi]),
                          BL_TO_FORTRAN_ANYD((*rho_next_fp_fft[lev])[mfi]));
    }
    BL_PROFILE_VAR_STOP(blp_push_eb);

    BL_PROFILE_VAR_START(blp_copy);
    CopyDataFromFFTToValid(*Efield_fp[lev][0], *Efield_fp_fft[lev][0], ba_valid_fp_fft[lev]);
    CopyDataFromFFTToValid(*Efield_fp[lev][1], *Efield_fp_fft[lev][1], ba_valid_fp_fft[lev]);
    CopyDataFromFFTToValid(*Efield_fp[lev][2], *Efield_fp_fft[lev][2], ba_valid_fp_fft[lev]);
    CopyDataFromFFTToValid(*Bfield_fp[lev][0], *Bfield_fp_fft[lev][0], ba_valid_fp_fft[lev]);
    CopyDataFromFFTToValid(*Bfield_fp[lev][1], *Bfield_fp_fft[lev][1], ba_valid_fp_fft[lev]);
    CopyDataFromFFTToValid(*Bfield_fp[lev][2], *Bfield_fp_fft[lev][2], ba_valid_fp_fft[lev]);
    BL_PROFILE_VAR_STOP(blp_copy);

    if (lev > 0)
    {
        amrex::Abort("WarpX::PushPSATD: TODO");
    }
}
