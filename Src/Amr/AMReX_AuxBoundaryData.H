
#ifndef AMREX_AuxBoundaryData_H_
#define AMREX_AuxBoundaryData_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace amrex {
// \cond CODEGEN
class AuxBoundaryData
{
public:

    AuxBoundaryData () noexcept;

    AuxBoundaryData (const BoxArray& grids,
                     int             n_grow,
                     int             n_comp,
                     const Geometry& geom);

    AuxBoundaryData (const AuxBoundaryData& rhs);

    AuxBoundaryData& operator= (const AuxBoundaryData& rhs) = delete;

    void copyTo (MultiFab& destmf,
                 int       src_comp,
                 int       dst_comp,
                 int       num_comp) const;

    void copyFrom (const MultiFab& srcmf,
                   int             src_comp,
                   int             dst_comp,
                   int             num_comp,
		   int             src_ng = 0);

    size_t size () const noexcept
    {
        BL_ASSERT(!m_empty); BL_ASSERT(m_initialized); return m_fabs.size();
    }

    void copy (const AuxBoundaryData& src,
               int                    src_comp,
               int                    dst_comp,
               int                    num_comp);

    void initialize (const BoxArray& grids,
		     int             n_grow,
                     int             n_comp,
                     const Geometry& geom);

    const BoxArray& equivBoxArray () const noexcept
    {
        BL_ASSERT(!m_empty); BL_ASSERT(m_initialized); return m_fabs.boxArray();
    }

    void setVal (Real r) { BL_ASSERT(m_initialized); if (!m_empty) m_fabs.setVal(r); }

    const DistributionMapping& DistributionMap () const noexcept
    {
        BL_ASSERT(!m_empty); BL_ASSERT(m_initialized); return m_fabs.DistributionMap();
    }

    FArrayBox&       operator[] (const MFIter& mfi) noexcept
    {
        BL_ASSERT(!m_empty); BL_ASSERT(m_initialized); return m_fabs[mfi];
    }
    const FArrayBox& operator[] (const MFIter& mfi) const noexcept
    {
        BL_ASSERT(!m_empty); BL_ASSERT(m_initialized); return m_fabs[mfi];
    }

    FArrayBox      * fabPtr (const MFIter& mfi) noexcept
    {
        BL_ASSERT(!m_empty); BL_ASSERT(m_initialized); return m_fabs.fabPtr(mfi);
    }
    FArrayBox const* fabPtr (const MFIter& mfi) const noexcept
    {
        BL_ASSERT(!m_empty); BL_ASSERT(m_initialized); return m_fabs.fabPtr(mfi);
    }

    int nGrow () const noexcept { BL_ASSERT(m_initialized); return m_ngrow; }

    int nComp () const noexcept
    {
        BL_ASSERT(!m_empty); BL_ASSERT(m_initialized); return m_fabs.nComp();
    }

    bool isEmpty () const noexcept { return m_empty; }

private:

    MultiFab m_fabs;
    int      m_ngrow;
    bool     m_empty;
    bool     m_initialized;
};
// \endcond

}

#endif
