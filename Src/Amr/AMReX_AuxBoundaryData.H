
#ifndef AMREX_AuxBoundaryData_H_
#define AMREX_AuxBoundaryData_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace amrex {
// \cond CODEGEN
class AuxBoundaryData
{
public:

    AuxBoundaryData () noexcept = default;

    AuxBoundaryData (const BoxArray& ba,
                     int             n_grow,
                     int             n_comp,
                     const Geometry& geom);

    ~AuxBoundaryData () = default;

    AuxBoundaryData (AuxBoundaryData&& rhs) = default;
    AuxBoundaryData& operator= (AuxBoundaryData&& rhs) = default;

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

    void initialize (const BoxArray& ba,
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

    int nGrow () const noexcept { BL_ASSERT(m_initialized); return m_ngrow; }

    int nComp () const noexcept
    {
        BL_ASSERT(!m_empty); BL_ASSERT(m_initialized); return m_fabs.nComp();
    }

    bool isEmpty () const noexcept { return m_empty; }

private:

    MultiFab m_fabs;
    int      m_ngrow{0};
    bool     m_empty{false};
    bool     m_initialized{false};
};
// \endcond

}

#endif
