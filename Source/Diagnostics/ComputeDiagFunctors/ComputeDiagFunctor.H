#ifndef WARPX_COMPUTEDIAGFUNCTOR_H_
#define WARPX_COMPUTEDIAGFUNCTOR_H_

#include <AMReX_MultiFab.H>

/**
 * \brief Functor to compute a diagnostic and store the result in existing
 * MultiFab
 */
class
ComputeDiagFunctor
{
public:
    ComputeDiagFunctor( int ncomp, amrex::IntVect crse_ratio) :
                        m_ncomp(ncomp), m_crse_ratio(crse_ratio) {};
    /** Compute a field and store the result in mf_dst
     * \param[out] mf_dst output MultiFab where the result is written
     * \param[in] dcomp first component of mf_dst in which the result is written
     * \param[in] ncomp number of components of mf_dst to update
     * \param[in] crse_ratio for interpolating field values from simulation MultiFab, mf_src,
                  to the output diagnostic MultiFab, mf_dst.
     */
    virtual void operator() (amrex::MultiFab& mf_dst, int dcomp) const = 0;
    /** Number of component from the input multifab to write to the output
     * multifab */
    int nComp () const { return m_ncomp; };
private:
    /** Number of components of mf_dst that this functor updates. */
    int m_ncomp;
protected:
    /** Coarsening ratio used to interpolate fields from simulation MultiFabs to output MultiFab. */
    amrex::IntVect m_crse_ratio;
};

#endif // WARPX_COMPUTEDIAGFUNCTOR_H_
