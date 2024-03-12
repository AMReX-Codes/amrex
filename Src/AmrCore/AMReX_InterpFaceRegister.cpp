#include <AMReX_InterpFaceRegister.H>
#include <AMReX_InterpFaceReg_C.H>

namespace amrex {

InterpFaceRegister::InterpFaceRegister (BoxArray const& fba, DistributionMapping const& fdm,
                                        Geometry const& fgeom, IntVect const& ref_ratio)
{
    define(fba, fdm, fgeom, ref_ratio);
}

void InterpFaceRegister::define (BoxArray const& fba, DistributionMapping const& fdm,
                                 Geometry const& fgeom, IntVect const& ref_ratio)
{
    AMREX_ASSERT(fba.ixType().cellCentered());
    m_fine_ba = fba;
    m_fine_dm = fdm;
    m_fine_geom = fgeom;
    m_ref_ratio = ref_ratio;

    m_crse_geom = amrex::coarsen(m_fine_geom, m_ref_ratio);

    constexpr int crse_fine_face = 1;
    constexpr int fine_fine_face = 0;
    // First, set the face mask to 1 (i.e., coarse/fine boundary),
    // except that it's 0 for non-periodic phys boundary.  For periodic
    // boundary, it's set to 1 here, and it will be fixed later.
    for (OrientationIter it; it; ++it) {
        Orientation face = it();
        int idim = face.coordDir();
        IndexType t;
        t.set(idim);
        m_fine_face_ba[face] = BoxArray(m_fine_ba, BATransformer(face,t,0,0,0));
        m_crse_face_ba[face] = amrex::coarsen(m_fine_face_ba[face], m_ref_ratio);
        m_face_mask[face].define(m_fine_face_ba[face], m_fine_dm, 1, 0);

        auto const& domface = amrex::bdryNode(m_fine_geom.Domain(), face);

#ifdef AMREX_USE_GPU
        Vector<Array4BoxValTag<int> > tags;
#endif
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(m_face_mask[face]); mfi.isValid(); ++mfi) {
            auto& fab = m_face_mask[face][mfi];
            Box const& dbox = fab.box();
            int value = (m_fine_geom.isPeriodic(idim) || !dbox.intersects(domface))
                ? crse_fine_face : fine_fine_face;
#ifdef AMREX_USE_GPU
            if (Gpu::inLaunchRegion()) {
                tags.emplace_back(Array4BoxValTag<int>{fab.array(),dbox,value});
            } else
#endif
            {
                fab.template setVal<RunOn::Host>(value, dbox);
            }
        }

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            ParallelFor(tags, [=] AMREX_GPU_DEVICE (int i, int j, int k,
                                                    Array4BoxValTag<int> const& tag) noexcept
            {
                tag.dfab(i,j,k) = tag.val;
            });
        }
#endif
    }

    // This set fine/fine boundary, including at periodic boundary.
    for (OrientationIter it; it; ++it) {
        Orientation face = it();
        Orientation other_face = face.flip();
        auto const& cpc = m_face_mask[face].getCPC(IntVect(0), m_face_mask[other_face],
                                                   IntVect(0), m_fine_geom.periodicity());
        m_face_mask[face].setVal(fine_fine_face, cpc, 0, 1);
    }
}

iMultiFab const&
InterpFaceRegister::mask (Orientation face) const
{
    return m_face_mask[face];
}

#ifdef AMREX_USE_GPU
namespace {
    struct IFRTag {
        Array4<Real> fine;
        Array4<Real> slope;
        Array4<Real const> crse;
        Array4<int const> mask;
        Box domface;
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Box box() const noexcept { return Box(mask); }
    };
}
#endif

void
InterpFaceRegister::interp (Array<MultiFab*, AMREX_SPACEDIM> const& fine, // NOLINT(readability-convert-member-functions-to-static)
                            Array<MultiFab const*, AMREX_SPACEDIM> const& crse,
                            int scomp, int ncomp)
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(fine,crse,scomp,ncomp);
    amrex::Abort("InterpFaceRegister::interp: 1D not supported");
#else
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Orientation olo(idim, Orientation::low);
        Orientation ohi(idim, Orientation::high);
        IntVect ng(1);
        ng[idim] = 0;
        MultiFab clodata(m_crse_face_ba[olo], m_fine_dm, ncomp, ng);
        MultiFab chidata(m_crse_face_ba[ohi], m_fine_dm, ncomp, ng);
        clodata.ParallelCopy_nowait(*crse[idim], scomp, 0, ncomp, IntVect(0), ng,
                                    m_crse_geom.periodicity());
        chidata.ParallelCopy_nowait(*crse[idim], scomp, 0, ncomp, IntVect(0), ng,
                                    m_crse_geom.periodicity());
        clodata.ParallelCopy_finish();
        chidata.ParallelCopy_finish();

        IntVect rr = m_ref_ratio;

        Box const& domain = m_crse_geom.growPeriodicDomain(1);
        Box const& domlo = amrex::bdryLo(domain, idim);
        Box const& domhi = amrex::bdryHi(domain, idim);

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            auto const& mlo_mf = m_face_mask[Orientation(idim,Orientation::low)];
            auto const& mhi_mf = m_face_mask[Orientation(idim,Orientation::high)];
            MultiFab slope_lo(mlo_mf.boxArray(), mlo_mf.DistributionMap(), ncomp, 0);
            MultiFab slope_hi(mhi_mf.boxArray(), mhi_mf.DistributionMap(), ncomp, 0);

            Vector<IFRTag> tags;
            for (MFIter mfi(*fine[idim]); mfi.isValid(); ++mfi) {
                auto const& fine_arr = fine[idim]->array(mfi);
                auto const& slo_arr = slope_lo.array(mfi);
                auto const& shi_arr = slope_hi.array(mfi);
                auto const& clo_arr = clodata.const_array(mfi);
                auto const& chi_arr = chidata.const_array(mfi);
                auto const& mlo_arr = mlo_mf.const_array(mfi);
                auto const& mhi_arr = mhi_mf.const_array(mfi);
                tags.push_back(IFRTag{fine_arr, slo_arr, clo_arr, mlo_arr, domlo});
                tags.push_back(IFRTag{fine_arr, shi_arr, chi_arr, mhi_arr, domhi});
            }

            ParallelFor(tags, [=] AMREX_GPU_DEVICE (int i, int j, int k, IFRTag const& tag) noexcept
            {
                if (tag.mask(i,j,k)) {
                    interp_face_reg(AMREX_D_DECL(i,j,k), rr, tag.fine, scomp, tag.crse,
                                    tag.slope, ncomp, tag.domface, idim);
                }
            });
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
            {
                FArrayBox slope;
                for (MFIter mfi(*fine[idim]); mfi.isValid(); ++mfi) {
                    Box const& vbx = mfi.validbox();
                    Box const& lobx = amrex::bdryLo(vbx, idim);
                    Box const& hibx = amrex::bdryHi(vbx, idim);
                    auto const& fine_arr = fine[idim]->array(mfi);
                    auto const& clo_arr = clodata.const_array(mfi);
                    auto const& chi_arr = chidata.const_array(mfi);
                    auto const& mlo_arr =
                        m_face_mask[Orientation(idim,Orientation::low)].const_array(mfi);
                    auto const& mhi_arr =
                        m_face_mask[Orientation(idim,Orientation::high)].const_array(mfi);
                    slope.resize(lobx, ncomp);
                    Array4<Real> slope_arr = slope.array();
                    amrex::LoopOnCpu(lobx, [&] (int i, int j, int k)
                    {
                        if (mlo_arr(i,j,k)) {
                            interp_face_reg(AMREX_D_DECL(i,j,k), rr, fine_arr, scomp, clo_arr,
                                            slope_arr, ncomp, domlo, idim);
                        }
                    });
                    slope.resize(hibx, ncomp);
                    slope_arr = slope.array();
                    amrex::LoopOnCpu(hibx, [&] (int i, int j, int k)
                    {
                        if (mhi_arr(i,j,k)) {
                            interp_face_reg(AMREX_D_DECL(i,j,k), rr, fine_arr, scomp, chi_arr,
                                            slope_arr, ncomp, domhi, idim);
                        }
                    });
                }
            }
        }
    }
#endif
}

}
