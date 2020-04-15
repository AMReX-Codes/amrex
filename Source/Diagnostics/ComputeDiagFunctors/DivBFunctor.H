#ifndef WARPX_DIVBFUNCTOR_H_
#define WARPX_DIVBFUNCTOR_H_

#include "ComputeDiagFunctor.H"

/**
 * \brief Functor to compute divB into mf_out.
 */
class
DivBFunctor final : public ComputeDiagFunctor
{
public:
    /** Constructor.
     * \param[in] arr_mf_src source multifabs (3 elements for x y z).
     * \param[in] lev level of multifab.
     * \param[in] crse_ratio for interpolating field values from simulation MultiFabs
                  to the output diagnostic MultiFab, mf_dst.
     * \param[in] ncomp Number of component of mf_src to cell-center in dst multifab.
     */
    DivBFunctor(const std::array<const amrex::MultiFab* const, 3> arr_mf_src, const int lev, amrex::IntVect crse_ratio, const int ncomp=1);

    /** \brief Compute DivB directly into mf_dst.
     *
     * \param[out] mf_dst output MultiFab where the result is written
     * \param[in] dcomp first component of mf_dst in which cell-centered
     *            data is stored
     */
    virtual void operator()(amrex::MultiFab& mf_dst, const int dcomp) const override;
private:
    /** Vector of pointer to source multifab Bx, By, Bz */
    std::array<const amrex::MultiFab * const, 3> m_arr_mf_src;
    int const m_lev; /**< level on which mf_src is defined (used in cylindrical) */
};

#endif // WARPX_DIVBFUNCTOR_H_
