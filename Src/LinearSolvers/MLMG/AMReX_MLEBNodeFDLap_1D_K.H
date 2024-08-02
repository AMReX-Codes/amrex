#ifndef AMREX_MLEBNODEFDLAP_1D_K_H_
#define AMREX_MLEBNODEFDLAP_1D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_adotx (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& /*y*/,
                        Array4<Real const> const& /*x*/, Array4<int const> const& /*dmsk*/,
                        Real /*bx*/) noexcept
{
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_gsrb (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& /*x*/,
                       Array4<Real const> const& /*rhs*/, Array4<int const> const& /*dmsk*/,
                       Real /*bx*/, int /*redblack*/) noexcept
{
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_adotx (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& /*y*/,
                            Array4<Real const> const& /*x*/, Array4<int const> const& /*dmsk*/,
                            Array4<Real const> const& /* sig */, Real /*bx*/) noexcept
{
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_sig_gsrb (int /*i*/, int /*j*/, int /*k*/, Array4<Real> const& /*x*/,
                           Array4<Real const> const& /*rhs*/, Array4<int const> const& /*dmsk*/,
                           Array4<Real const> const& /* sig */,
                           Real /*bx*/, int /*redblack*/) noexcept
{
}

}

#endif
