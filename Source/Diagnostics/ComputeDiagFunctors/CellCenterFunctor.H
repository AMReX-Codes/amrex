#ifndef WARPX_CELLCENTERFUNCTOR_H_
#define WARPX_CELLCENTERFUNCTOR_H_

#include "ComputeDiagFunctor.H"

/**
 * \brief Functor to cell-center MF and store result in mf_out.
 */
class
CellCenterFunctor final : public ComputeDiagFunctor
{
public:
    /** Constructor.
     * \param[in] mf_src source multifab.
     * \param[in] lev level of multifab. Used for averaging in rz.
     * \param[in] convert_RZmodes2cartesian (in cylindrical) whether to
     *            sum all modes in mf_src before cell-centering into dst multifab.
     * \param[in] ncomp Number of component of mf_src to cell-center in dst multifab.
     */
    CellCenterFunctor(amrex::MultiFab const * mf_src, int lev,
                      bool convertRZmodes2cartesian=true, int ncomp=1);

    /** \brief Cell-center m_mf_src and write the result in mf_dst.
     *
     * In cylindrical geometry, by default this functor average all components
     * of a MultiFab and writes into one single component.
     *
     * \param[out] mf_dst output MultiFab where the result is written
     * \param[in] dcomp first component of mf_dst in which cell-centered
     *            data is stored
     */
    virtual void operator()(amrex::MultiFab& mf_dst, int dcomp) const override;
private:
    /** pointer to source multifab (can be multi-component) */
    amrex::MultiFab const * const m_mf_src;
    int m_lev; /**< level on which mf_src is defined (used in cylindrical) */
    /**< (for cylindrical) whether to average all modes into 1 comp */
    bool m_convertRZmodes2cartesian;
};

#endif // WARPX_CELLCENTERFUNCTOR_H_
