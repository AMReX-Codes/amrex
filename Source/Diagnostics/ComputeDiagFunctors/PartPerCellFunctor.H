#ifndef WARPX_PARTPERCELLFUNCTOR_H_
#define WARPX_PARTPERCELLFUNCTOR_H_

#include "ComputeDiagFunctor.H"

/**
 * \brief Functor to cell-center MF and store result in mf_out.
 */
class
PartPerCellFunctor final : public ComputeDiagFunctor
{
public:
    /** Constructor.
     * \param[in] mf_src source multifab. Must be nullptr as not source MF is needed
     *            to compute the number of particles per cell, banane.
     * \param[in] lev level of multifab. Used for averaging in rz.
     * \param[in] ncomp Number of component of mf_src to cell-center in dst multifab.
     */
    PartPerCellFunctor(const amrex::MultiFab * const mf_src, const int lev, const int ncomp=1);

    /** \brief Compute the number of particles per cell directly into mf_dst.
     *
     * \param[out] mf_dst output MultiFab where the result is written
     * \param[in] dcomp first component of mf_dst in which cell-centered
     *            data is stored
     */
    virtual void operator()(amrex::MultiFab& mf_dst, const int dcomp) const override;
private:
    int const m_lev; /**< level on which mf_src is defined */
};

#endif // WARPX_PARTPERCELLFUNCTOR_H_
