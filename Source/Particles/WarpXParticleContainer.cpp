#include <limits>

#include <MultiParticleContainer.H>
#include <WarpXParticleContainer.H>
#include <AMReX_AmrParGDB.H>
#include <WarpXComm.H>
#include <WarpX_f.H>
#include <WarpX.H>
#include <WarpXAlgorithmSelection.H>
#include <WarpXComm.H>
// Import low-level single-particle kernels
#include <GetAndSetPosition.H>
#include <UpdatePosition.H>
#include <CurrentDeposition.H>
#include <ChargeDeposition.H>

using namespace amrex;

int WarpXParticleContainer::do_not_push = 0;
int WarpXParticleContainer::do_not_deposit = 0;

WarpXParIter::WarpXParIter (ContainerType& pc, int level)
    : ParIter(pc, level, MFItInfo().SetDynamic(WarpX::do_dynamic_scheduling))
{
}

#if (AMREX_SPACEDIM == 2)
void
WarpXParIter::GetPosition (Gpu::ManagedDeviceVector<ParticleReal>& xp,
                           Gpu::ManagedDeviceVector<ParticleReal>& yp,
                           Gpu::ManagedDeviceVector<ParticleReal>& zp) const
{
    amrex::ParIter<0,0,PIdx::nattribs>::GetPosition(xp, zp);
#ifdef WARPX_DIM_RZ
    const auto& attribs = GetAttribs();
    const auto& thetap = attribs[PIdx::theta];
    yp.resize(xp.size());
    ParticleReal* const AMREX_RESTRICT x = xp.dataPtr();
    ParticleReal* const AMREX_RESTRICT y = yp.dataPtr();
    const ParticleReal* const AMREX_RESTRICT theta = thetap.dataPtr();
    amrex::ParallelFor( xp.size(),
        [=] AMREX_GPU_DEVICE (long i) {
        // The x stored in the particles is actually the radius
        y[i] = x[i]*std::sin(theta[i]);
        x[i] = x[i]*std::cos(theta[i]);
    });
#else
    yp.resize(xp.size(), std::numeric_limits<ParticleReal>::quiet_NaN());
#endif
}

void
WarpXParIter::SetPosition (const Gpu::ManagedDeviceVector<ParticleReal>& xp,
                           const Gpu::ManagedDeviceVector<ParticleReal>& yp,
                           const Gpu::ManagedDeviceVector<ParticleReal>& zp)
{
#ifdef WARPX_DIM_RZ
    auto& attribs = GetAttribs();
    auto& thetap = attribs[PIdx::theta];
    Gpu::ManagedDeviceVector<ParticleReal> rp(xp.size());
    const ParticleReal* const AMREX_RESTRICT x = xp.dataPtr();
    const ParticleReal* const AMREX_RESTRICT y = yp.dataPtr();
    ParticleReal* const AMREX_RESTRICT r = rp.dataPtr();
    ParticleReal* const AMREX_RESTRICT theta = thetap.dataPtr();
    amrex::ParallelFor( xp.size(),
        [=] AMREX_GPU_DEVICE (long i) {
        theta[i] = std::atan2(y[i], x[i]);
        r[i] = std::sqrt(x[i]*x[i] + y[i]*y[i]);
    });
    amrex::ParIter<0,0,PIdx::nattribs>::SetPosition(rp, zp);
#else
    amrex::ParIter<0,0,PIdx::nattribs>::SetPosition(xp, zp);
#endif
}
#endif

WarpXParticleContainer::WarpXParticleContainer (AmrCore* amr_core, int ispecies)
    : ParticleContainer<0,0,PIdx::nattribs>(amr_core->GetParGDB())
    , species_id(ispecies)
{
    for (unsigned int i = PIdx::Ex; i <= PIdx::Bz; ++i) {
        communicate_real_comp[i] = false; // Don't need to communicate E and B.
    }
    SetParticleSize();
    ReadParameters();

    // build up the map of string names to particle component numbers
    particle_comps["w"]  = PIdx::w;
    particle_comps["ux"] = PIdx::ux;
    particle_comps["uy"] = PIdx::uy;
    particle_comps["uz"] = PIdx::uz;
    particle_comps["Ex"] = PIdx::Ex;
    particle_comps["Ey"] = PIdx::Ey;
    particle_comps["Ez"] = PIdx::Ez;
    particle_comps["Bx"] = PIdx::Bx;
    particle_comps["By"] = PIdx::By;
    particle_comps["Bz"] = PIdx::Bz;
#ifdef WARPX_DIM_RZ
    particle_comps["theta"] = PIdx::theta;
#endif

    // Initialize temporary local arrays for charge/current deposition
    int num_threads = 1;
    #ifdef _OPENMP
    #pragma omp parallel
    #pragma omp single
    num_threads = omp_get_num_threads();
    #endif
    local_rho.resize(num_threads);
    local_jx.resize(num_threads);
    local_jy.resize(num_threads);
    local_jz.resize(num_threads);
    m_xp.resize(num_threads);
    m_yp.resize(num_threads);
    m_zp.resize(num_threads);
}

void
WarpXParticleContainer::ReadParameters ()
{
    static bool initialized = false;
    if (!initialized)
    {
        ParmParse pp("particles");

#ifdef AMREX_USE_GPU
        do_tiling = false; // By default, tiling is off on GPU
#else
        do_tiling = true;
#endif
        pp.query("do_tiling",  do_tiling);
        pp.query("do_not_push", do_not_push);

        initialized = true;
    }
}

void
WarpXParticleContainer::AllocData ()
{
    // have to resize here, not in the constructor because grids have not
    // been built when constructor was called.
    reserveData();
    resizeData();
}

void
WarpXParticleContainer::AddNParticles (int lev,
                                       int n, const ParticleReal* x, const ParticleReal* y, const ParticleReal* z,
                                       const ParticleReal* vx, const ParticleReal* vy, const ParticleReal* vz,
                                       int nattr, const ParticleReal* attr, int uniqueparticles, int id)
{
    BL_ASSERT(nattr == 1);
    const ParticleReal* weight = attr;

    int ibegin, iend;
    if (uniqueparticles) {
        ibegin = 0;
        iend = n;
    } else {
        int myproc = ParallelDescriptor::MyProc();
        int nprocs = ParallelDescriptor::NProcs();
        int navg = n/nprocs;
        int nleft = n - navg * nprocs;
        if (myproc < nleft) {
            ibegin = myproc*(navg+1);
            iend = ibegin + navg+1;
        } else {
            ibegin = myproc*navg + nleft;
            iend = ibegin + navg;
        }
    }

    //  Add to grid 0 and tile 0
    // Redistribute() will move them to proper places.
    std::pair<int,int> key {0,0};
    auto& particle_tile = DefineAndReturnParticleTile(0, 0, 0);

    std::size_t np = iend-ibegin;

#ifdef WARPX_DIM_RZ
    Vector<ParticleReal> theta(np);
#endif

    for (int i = ibegin; i < iend; ++i)
    {
        ParticleType p;
        if (id==-1)
        {
            p.id() = ParticleType::NextID();
        } else {
            p.id() = id;
        }
        p.cpu() = ParallelDescriptor::MyProc();
#if (AMREX_SPACEDIM == 3)
        p.pos(0) = x[i];
        p.pos(1) = y[i];
        p.pos(2) = z[i];
#elif (AMREX_SPACEDIM == 2)
#ifdef WARPX_DIM_RZ
        theta[i-ibegin] = std::atan2(y[i], x[i]);
        p.pos(0) = std::sqrt(x[i]*x[i] + y[i]*y[i]);
#else
        p.pos(0) = x[i];
#endif
        p.pos(1) = z[i];
#endif

        if ( (NumRuntimeRealComps()>0) || (NumRuntimeIntComps()>0) ){
            auto& ptile = DefineAndReturnParticleTile(0, 0, 0);
        }

        particle_tile.push_back(p);
    }

    if (np > 0)
    {
        particle_tile.push_back_real(PIdx::w , weight + ibegin, weight + iend);
        particle_tile.push_back_real(PIdx::ux,     vx + ibegin,     vx + iend);
        particle_tile.push_back_real(PIdx::uy,     vy + ibegin,     vy + iend);
        particle_tile.push_back_real(PIdx::uz,     vz + ibegin,     vz + iend);

        if ( (NumRuntimeRealComps()>0) || (NumRuntimeIntComps()>0) ){
            auto& ptile = DefineAndReturnParticleTile(0, 0, 0);
        }

        for (int comp = PIdx::uz+1; comp < PIdx::nattribs; ++comp)
        {
#ifdef WARPX_DIM_RZ
            if (comp == PIdx::theta) {
                particle_tile.push_back_real(comp, theta.data(), theta.data() + np);
            }
            else {
                particle_tile.push_back_real(comp, np, 0.0);
            }
#else
            particle_tile.push_back_real(comp, np, 0.0);
#endif
        }

        for (int i = PIdx::nattribs; i < NumRealComps(); ++i)
        {
            particle_tile.push_back_real(i, 0.0);
        }
    }

    Redistribute();
}

/* \brief Current Deposition for thread thread_num
 * \param pti         : Particle iterator
 * \param wp          : Array of particle weights
 * \param uxp uyp uzp : Array of particle
 * \param ion_lev      : Pointer to array of particle ionization level. This is
                         required to have the charge of each macroparticle
                         since q is a scalar. For non-ionizable species,
                         ion_lev is a null pointer.
 * \param jx jy jz    : Full array of current density
 * \param offset      : Index of first particle for which current is deposited
 * \param np_to_depose: Number of particles for which current is deposited.
                        Particles [offset,offset+np_tp_depose] deposit current
 * \param thread_num  : Thread number (if tiling)
 * \param lev         : Level of box that contains particles
 * \param depos_lev   : Level on which particles deposit (if buffers are used)
 * \param dt          : Time step for particle level
 */
void
WarpXParticleContainer::DepositCurrent(WarpXParIter& pti,
                                       RealVector& wp, RealVector& uxp,
                                       RealVector& uyp, RealVector& uzp,
                                       const int * const ion_lev,
                                       MultiFab* jx, MultiFab* jy, MultiFab* jz,
                                       const long offset, const long np_to_depose,
                                       int thread_num, int lev, int depos_lev,
                                       Real dt)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE((depos_lev==(lev-1)) ||
                                     (depos_lev==(lev  )),
                                     "Deposition buffers only work for lev-1");
    // If no particles, do not do anything
    if (np_to_depose == 0) return;

    // If user decides not to deposit
    if (do_not_deposit) return;

    const long ngJ = jx->nGrow();
    const std::array<Real,3>& dx = WarpX::CellSize(std::max(depos_lev,0));
    Real q = this->charge;

    BL_PROFILE_VAR_NS("PPC::Evolve::Accumulate", blp_accumulate);
    BL_PROFILE_VAR_NS("PPC::CurrentDeposition", blp_deposit);


    // Get tile box where current is deposited.
    // The tile box is different when depositing in the buffers (depos_lev<lev)
    // or when depositing inside the level (depos_lev=lev)
    Box tilebox;
    if (lev == depos_lev) {
        tilebox = pti.tilebox();
    } else {
        const IntVect& ref_ratio = WarpX::RefRatio(depos_lev);
        tilebox = amrex::coarsen(pti.tilebox(),ref_ratio);
    }

    // Staggered tile boxes (different in each direction)
    Box tbx = convert(tilebox, WarpX::jx_nodal_flag);
    Box tby = convert(tilebox, WarpX::jy_nodal_flag);
    Box tbz = convert(tilebox, WarpX::jz_nodal_flag);
    tilebox.grow(ngJ);

#ifdef AMREX_USE_GPU
    // No tiling on GPU: jx_ptr points to the full
    // jx array (same for jy_ptr and jz_ptr).
    auto & jx_fab = jx->get(pti);
    auto & jy_fab = jy->get(pti);
    auto & jz_fab = jz->get(pti);
    Array4<Real> const& jx_arr = jx->array(pti);
    Array4<Real> const& jy_arr = jy->array(pti);
    Array4<Real> const& jz_arr = jz->array(pti);
#else
    // Tiling is on: jx_ptr points to local_jx[thread_num]
    // (same for jy_ptr and jz_ptr)
    tbx.grow(ngJ);
    tby.grow(ngJ);
    tbz.grow(ngJ);

    local_jx[thread_num].resize(tbx, jx->nComp());
    local_jy[thread_num].resize(tby, jy->nComp());
    local_jz[thread_num].resize(tbz, jz->nComp());

    // local_jx[thread_num] is set to zero
    local_jx[thread_num].setVal(0.0);
    local_jy[thread_num].setVal(0.0);
    local_jz[thread_num].setVal(0.0);

    auto & jx_fab = local_jx[thread_num];
    auto & jy_fab = local_jy[thread_num];
    auto & jz_fab = local_jz[thread_num];
    Array4<Real> const& jx_arr = local_jx[thread_num].array();
    Array4<Real> const& jy_arr = local_jy[thread_num].array();
    Array4<Real> const& jz_arr = local_jz[thread_num].array();
#endif
    // GPU, no tiling: deposit directly in jx
    // CPU, tiling: deposit into local_jx
    // (same for jx and jz)

    ParticleReal* AMREX_RESTRICT xp = m_xp[thread_num].dataPtr() + offset;
    ParticleReal* AMREX_RESTRICT zp = m_zp[thread_num].dataPtr() + offset;
    ParticleReal* AMREX_RESTRICT yp = m_yp[thread_num].dataPtr() + offset;

    // Lower corner of tile box physical domain
    // Note that this includes guard cells since it is after tilebox.ngrow
    const Dim3 lo = lbound(tilebox);
    const std::array<Real, 3>& xyzmin = WarpX::LowerCorner(tilebox, depos_lev);

    BL_PROFILE_VAR_START(blp_deposit);
    if (WarpX::current_deposition_algo == CurrentDepositionAlgo::Esirkepov) {
        if        (WarpX::nox == 1){
            doEsirkepovDepositionShapeN<1>(
                xp, yp, zp, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                jx_arr, jy_arr, jz_arr, np_to_depose, dt, dx, xyzmin, lo, q,
                WarpX::n_rz_azimuthal_modes);
        } else if (WarpX::nox == 2){
            doEsirkepovDepositionShapeN<2>(
                xp, yp, zp, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                jx_arr, jy_arr, jz_arr, np_to_depose, dt, dx, xyzmin, lo, q,
                WarpX::n_rz_azimuthal_modes);
        } else if (WarpX::nox == 3){
            doEsirkepovDepositionShapeN<3>(
                xp, yp, zp, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                jx_arr, jy_arr, jz_arr, np_to_depose, dt, dx, xyzmin, lo, q,
                WarpX::n_rz_azimuthal_modes);
        }
    } else {
        if        (WarpX::nox == 1){
            doDepositionShapeN<1>(
                xp, yp, zp, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                jx_fab, jy_fab, jz_fab, np_to_depose, dt, dx,
                xyzmin, lo, q);
        } else if (WarpX::nox == 2){
            doDepositionShapeN<2>(
                xp, yp, zp, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                jx_fab, jy_fab, jz_fab, np_to_depose, dt, dx,
                xyzmin, lo, q);
        } else if (WarpX::nox == 3){
            doDepositionShapeN<3>(
                xp, yp, zp, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                jx_fab, jy_fab, jz_fab, np_to_depose, dt, dx,
                xyzmin, lo, q);
        }
    }
    BL_PROFILE_VAR_STOP(blp_deposit);

#ifndef AMREX_USE_GPU
    BL_PROFILE_VAR_START(blp_accumulate);
    // CPU, tiling: atomicAdd local_jx into jx
    // (same for jx and jz)
    (*jx)[pti].atomicAdd(local_jx[thread_num], tbx, tbx, 0, 0, jx->nComp());
    (*jy)[pti].atomicAdd(local_jy[thread_num], tby, tby, 0, 0, jy->nComp());
    (*jz)[pti].atomicAdd(local_jz[thread_num], tbz, tbz, 0, 0, jz->nComp());
    BL_PROFILE_VAR_STOP(blp_accumulate);
#endif
}

/* \brief Charge Deposition for thread thread_num
 * \param pti         : Particle iterator
 * \param wp          : Array of particle weights
 * \param ion_lev     : Pointer to array of particle ionization level. This is
                         required to have the charge of each macroparticle
                         since q is a scalar. For non-ionizable species,
                         ion_lev is a null pointer.
 * \param rho         : Full array of charge density
 * \param icomp       : Component of rho into which charge is deposited.
                        0: old value (before particle push).
                        1: new value (after particle push).
 * \param offset      : Index of first particle for which charge is deposited
 * \param np_to_depose: Number of particles for which charge is deposited.
                        Particles [offset,offset+np_tp_depose] deposit charge
 * \param thread_num  : Thread number (if tiling)
 * \param lev         : Level of box that contains particles
 * \param depos_lev   : Level on which particles deposit (if buffers are used)
 */
void
WarpXParticleContainer::DepositCharge (WarpXParIter& pti, RealVector& wp,
                                       const int * const ion_lev,
                                       amrex::MultiFab* rho, int icomp,
                                       const long offset, const long np_to_depose,
                                       int thread_num, int lev, int depos_lev)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE((depos_lev==(lev-1)) ||
                                     (depos_lev==(lev  )),
                                     "Deposition buffers only work for lev-1");

    // If no particles, do not do anything
    if (np_to_depose == 0) return;

    // If user decides not to deposit
    if (do_not_deposit) return;

    const long ngRho = rho->nGrow();
    const std::array<Real,3>& dx = WarpX::CellSize(std::max(depos_lev,0));
    const Real q = this->charge;

    BL_PROFILE_VAR_NS("PPC::ChargeDeposition", blp_ppc_chd);
    BL_PROFILE_VAR_NS("PPC::Evolve::Accumulate", blp_accumulate);

    // Get tile box where charge is deposited.
    // The tile box is different when depositing in the buffers (depos_lev<lev)
    // or when depositing inside the level (depos_lev=lev)
    Box tilebox;
    if (lev == depos_lev) {
        tilebox = pti.tilebox();
    } else {
        const IntVect& ref_ratio = WarpX::RefRatio(depos_lev);
        tilebox = amrex::coarsen(pti.tilebox(),ref_ratio);
    }

    tilebox.grow(ngRho);

    const int nc = (rho->nComp() == 1 ? 1 : rho->nComp()/2);

#ifdef AMREX_USE_GPU
    // No tiling on GPU: rho_arr points to the full rho array.
    MultiFab rhoi(*rho, amrex::make_alias, icomp*nc, nc);
    Array4<Real> const& rho_arr = rhoi.array(pti);
#else
    // Tiling is on: rho_arr points to local_rho[thread_num]
    const Box tb = amrex::convert(tilebox, IntVect::TheUnitVector());

    local_rho[thread_num].resize(tb, nc);

    // local_rho[thread_num] is set to zero
    local_rho[thread_num].setVal(0.0);

    Array4<Real> const& rho_arr = local_rho[thread_num].array();
#endif
    // GPU, no tiling: deposit directly in rho
    // CPU, tiling: deposit into local_rho

    ParticleReal* AMREX_RESTRICT xp = m_xp[thread_num].dataPtr() + offset;
    ParticleReal* AMREX_RESTRICT zp = m_zp[thread_num].dataPtr() + offset;
    ParticleReal* AMREX_RESTRICT yp = m_yp[thread_num].dataPtr() + offset;

    // Lower corner of tile box physical domain
    // Note that this includes guard cells since it is after tilebox.ngrow
    const std::array<Real, 3>& xyzmin = WarpX::LowerCorner(tilebox, depos_lev);
    // Indices of the lower bound
    const Dim3 lo = lbound(tilebox);

    BL_PROFILE_VAR_START(blp_ppc_chd);
    if        (WarpX::nox == 1){
        doChargeDepositionShapeN<1>(xp, yp, zp, wp.dataPtr()+offset, ion_lev,
                                    rho_arr, np_to_depose, dx, xyzmin, lo, q);
    } else if (WarpX::nox == 2){
        doChargeDepositionShapeN<2>(xp, yp, zp, wp.dataPtr()+offset, ion_lev,
                                    rho_arr, np_to_depose, dx, xyzmin, lo, q);
    } else if (WarpX::nox == 3){
        doChargeDepositionShapeN<3>(xp, yp, zp, wp.dataPtr()+offset, ion_lev,
                                    rho_arr, np_to_depose, dx, xyzmin, lo, q);
    }
    BL_PROFILE_VAR_STOP(blp_ppc_chd);

#ifndef AMREX_USE_GPU
    BL_PROFILE_VAR_START(blp_accumulate);

    (*rho)[pti].atomicAdd(local_rho[thread_num], tb, tb, 0, icomp*nc, nc);

    BL_PROFILE_VAR_STOP(blp_accumulate);
#endif
}

void
WarpXParticleContainer::DepositCharge (amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                                        bool local, bool reset,
                                        bool do_rz_volume_scaling)
{
    // Loop over the refinement levels
    int const finest_level = rho.size() - 1;
    for (int lev = 0; lev <= finest_level; ++lev) {

        // Reset the `rho` array if `reset` is True
        if (reset) rho[lev]->setVal(0.0, rho[lev]->nGrow());

        // Loop over particle tiles and deposit charge on each level
#ifdef _OPENMP
        #pragma omp parallel
        {
        int thread_num = omp_get_thread_num();
#else
        int thread_num = 0;
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            const long np = pti.numParticles();
            auto& wp = pti.GetAttribs(PIdx::w);

            pti.GetPosition(m_xp[thread_num], m_yp[thread_num], m_zp[thread_num]);

            int* AMREX_RESTRICT ion_lev;
            if (do_field_ionization){
                ion_lev = pti.GetiAttribs(particle_icomps["ionization_level"]).dataPtr();
            } else {
                ion_lev = nullptr;
            }

            DepositCharge(pti, wp, ion_lev, rho[lev].get(), 0, 0, np, thread_num, lev, lev);
        }
#ifdef _OPENMP
        }
#endif

#ifdef WARPX_DIM_RZ
        if (do_rz_volume_scaling) {
            WarpX::GetInstance().ApplyInverseVolumeScalingToChargeDensity(rho[lev].get(), lev);
        }
#endif

        // Exchange guard cells
        if (!local) rho[lev]->SumBoundary( m_gdb->Geom(lev).periodicity() );
    }

    // Now that the charge has been deposited at each level,
    // we average down from fine to crse
    for (int lev = finest_level - 1; lev >= 0; --lev) {
        const DistributionMapping& fine_dm = rho[lev+1]->DistributionMap();
        BoxArray coarsened_fine_BA = rho[lev+1]->boxArray();
        coarsened_fine_BA.coarsen(m_gdb->refRatio(lev));
        MultiFab coarsened_fine_data(coarsened_fine_BA, fine_dm, rho[lev+1]->nComp(), 0);
        coarsened_fine_data.setVal(0.0);

        int const refinement_ratio = 2;

        interpolateDensityFineToCoarse( *rho[lev+1], coarsened_fine_data, refinement_ratio );
        rho[lev]->ParallelAdd( coarsened_fine_data, m_gdb->Geom(lev).periodicity() );
    }
}

std::unique_ptr<MultiFab>
WarpXParticleContainer::GetChargeDensity (int lev, bool local)
{
    const auto& gm = m_gdb->Geom(lev);
    const auto& ba = m_gdb->ParticleBoxArray(lev);
    const auto& dm = m_gdb->DistributionMap(lev);
    BoxArray nba = ba;
    nba.surroundingNodes();

    const int ng = WarpX::nox;

    auto rho = std::unique_ptr<MultiFab>(new MultiFab(nba,dm,WarpX::ncomps,ng));
    rho->setVal(0.0);

#ifdef _OPENMP
#pragma omp parallel
    {
#endif
#ifdef _OPENMP
        int thread_num = omp_get_thread_num();
#else
        int thread_num = 0;
#endif

        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            const long np = pti.numParticles();
            auto& wp = pti.GetAttribs(PIdx::w);

            pti.GetPosition(m_xp[thread_num], m_yp[thread_num], m_zp[thread_num]);

            int* AMREX_RESTRICT ion_lev;
            if (do_field_ionization){
                ion_lev = pti.GetiAttribs(particle_icomps["ionization_level"]).dataPtr();
            } else {
                ion_lev = nullptr;
            }

            DepositCharge(pti, wp, ion_lev, rho.get(), 0, 0, np,
                          thread_num, lev, lev);
        }
#ifdef _OPENMP
    }
#endif

#ifdef WARPX_DIM_RZ
    WarpX::GetInstance().ApplyInverseVolumeScalingToChargeDensity(rho.get(), lev);
#endif

    if (!local) rho->SumBoundary(gm.periodicity());

    return rho;
}

Real WarpXParticleContainer::sumParticleCharge(bool local) {

    amrex::Real total_charge = 0.0;

    const int nLevels = finestLevel();
    for (int lev = 0; lev < nLevels; ++lev)
    {

#ifdef _OPENMP
#pragma omp parallel reduction(+:total_charge)
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& wp = pti.GetAttribs(PIdx::w);
            for (unsigned long i = 0; i < wp.size(); i++) {
                total_charge += wp[i];
            }
        }
    }

    if (!local) ParallelDescriptor::ReduceRealSum(total_charge);
    total_charge *= this->charge;
    return total_charge;
}

std::array<Real, 3> WarpXParticleContainer::meanParticleVelocity(bool local) {

    amrex::Real vx_total = 0.0;
    amrex::Real vy_total = 0.0;
    amrex::Real vz_total = 0.0;

    long np_total = 0;

    amrex::Real inv_clight_sq = 1.0/PhysConst::c/PhysConst::c;

    const int nLevels = finestLevel();
    for (int lev = 0; lev <= nLevels; ++lev) {

#ifdef _OPENMP
#pragma omp parallel reduction(+:vx_total, vy_total, vz_total, np_total)
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& ux = pti.GetAttribs(PIdx::ux);
            auto& uy = pti.GetAttribs(PIdx::uy);
            auto& uz = pti.GetAttribs(PIdx::uz);

            np_total += pti.numParticles();

            for (unsigned long i = 0; i < ux.size(); i++) {
                Real usq = (ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i])*inv_clight_sq;
                Real gaminv = 1.0/std::sqrt(1.0 + usq);
                vx_total += ux[i]*gaminv;
                vy_total += uy[i]*gaminv;
                vz_total += uz[i]*gaminv;
            }
        }
    }

    if (!local) {
        ParallelDescriptor::ReduceRealSum(vx_total);
        ParallelDescriptor::ReduceRealSum(vy_total);
        ParallelDescriptor::ReduceRealSum(vz_total);
        ParallelDescriptor::ReduceLongSum(np_total);
    }

    std::array<Real, 3> mean_v;
    if (np_total > 0) {
        mean_v[0] = vx_total / np_total;
        mean_v[1] = vy_total / np_total;
        mean_v[2] = vz_total / np_total;
    }

    return mean_v;
}

Real WarpXParticleContainer::maxParticleVelocity(bool local) {

    amrex::ParticleReal max_v = 0.0;

    const int nLevels = finestLevel();
    for (int lev = 0; lev <= nLevels; ++lev)
    {

#ifdef _OPENMP
#pragma omp parallel reduction(max:max_v)
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& ux = pti.GetAttribs(PIdx::ux);
            auto& uy = pti.GetAttribs(PIdx::uy);
            auto& uz = pti.GetAttribs(PIdx::uz);
            for (unsigned long i = 0; i < ux.size(); i++) {
                max_v = std::max(max_v, std::sqrt(ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i]));
            }
        }
    }

    if (!local) ParallelAllReduce::Max(max_v, ParallelDescriptor::Communicator());
    return max_v;
}

void
WarpXParticleContainer::PushXES (Real dt)
{
    BL_PROFILE("WPC::PushXES()");

    const int num_levels = finestLevel() + 1;

    for (int lev = 0; lev < num_levels; ++lev) {
        const auto& gm = m_gdb->Geom(lev);
        const RealBox& prob_domain = gm.ProbDomain();
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti) {
            auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const long np  = pti.numParticles();

            auto& attribs = pti.GetAttribs();
            auto& uxp = attribs[PIdx::ux];
            auto& uyp = attribs[PIdx::uy];
            auto& uzp = attribs[PIdx::uz];

            WRPX_PUSH_LEAPFROG_POSITIONS(particles.dataPtr(), nstride, np,
                                         uxp.dataPtr(), uyp.dataPtr(),
#if AMREX_SPACEDIM == 3
                                         uzp.dataPtr(),
#endif
                                         &dt,
                                         prob_domain.lo(), prob_domain.hi());
        }
    }
}

void
WarpXParticleContainer::PushX (amrex::Real dt)
{
    const int nLevels = finestLevel();
    for (int lev = 0; lev <= nLevels; ++lev) {
        PushX(lev, dt);
    }
}

void
WarpXParticleContainer::PushX (int lev, amrex::Real dt)
{
    BL_PROFILE("WPC::PushX()");

    if (do_not_push) return;

    MultiFab* cost = WarpX::getCosts(lev);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {

        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            Real wt = amrex::second();

            //
            // Particle Push
            //
            // Extract pointers to particle position and momenta, for this particle tile
            // - positions are stored as an array of struct, in `ParticleType`
            ParticleType * AMREX_RESTRICT pstructs = &(pti.GetArrayOfStructs()[0]);
            // - momenta are stored as a struct of array, in `attribs`
            auto& attribs = pti.GetAttribs();
            ParticleReal* AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
            ParticleReal* AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
            ParticleReal* AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();
#ifdef WARPX_DIM_RZ
            ParticleReal* AMREX_RESTRICT theta = attribs[PIdx::theta].dataPtr();
#endif
            // Loop over the particles and update their position
            amrex::ParallelFor( pti.numParticles(),
                [=] AMREX_GPU_DEVICE (long i) {
                    ParticleType& p = pstructs[i]; // Particle object that gets updated
                    ParticleReal x, y, z; // Temporary variables
#ifndef WARPX_DIM_RZ
                    GetPosition( x, y, z, p ); // Initialize x, y, z
                    UpdatePosition( x, y, z, ux[i], uy[i], uz[i], dt);
                    SetPosition( p, x, y, z ); // Update the object p
#else
                    // For WARPX_DIM_RZ, the particles are still pushed in 3D Cartesian
                    GetCartesianPositionFromCylindrical( x, y, z, p, theta[i] );
                    UpdatePosition( x, y, z, ux[i], uy[i], uz[i], dt);
                    SetCylindricalPositionFromCartesian( p, theta[i], x, y, z );
#endif
                }
            );

            if (cost) {
                const Box& tbx = pti.tilebox();
                wt = (amrex::second() - wt) / tbx.d_numPts();
                Array4<Real> const& costarr = cost->array(pti);
                amrex::ParallelFor(tbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    costarr(i,j,k) += wt;
                });
            }
        }
    }
}

// This function is called in Redistribute, just after locate
void
WarpXParticleContainer::particlePostLocate(ParticleType& p,
                                           const ParticleLocData& pld,
                                           const int lev)
{
    // Tag particle if goes to higher level.
    // It will be split later in the loop
    if (pld.m_lev == lev+1
        and p.m_idata.id != NoSplitParticleID
        and p.m_idata.id >= 0)
    {
        p.m_idata.id = DoSplitParticleID;
    }

    if (pld.m_lev == lev-1){
        // For the moment, do not do anything if particles goes
        // to lower level.
    }
}
