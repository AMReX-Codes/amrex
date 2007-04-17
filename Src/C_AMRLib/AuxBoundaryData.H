
#ifndef _AuxBoundaryData_H_
#define _AuxBoundaryData_H_

#include <Geometry.H>
#include <MultiFab.H>

class AuxBoundaryData
{
public:

    AuxBoundaryData ();

    AuxBoundaryData (const BoxArray& grids,
                     int             n_grow,
                     int             n_comp,
                     const Geometry& geom);

    AuxBoundaryData (const AuxBoundaryData& rhs);

    void copyTo (MultiFab& destmf,
                 int       src_comp,
                 int       dst_comp,
                 int       num_comp) const;

    void copyFrom (const MultiFab& srcmf,
                   int             src_comp,
                   int             dst_comp,
                   int             num_comp);

    size_t size () const;

    void copy (const AuxBoundaryData& src,
               int                    src_comp,
               int                    dst_comp,
               int                    num_comp);

    void initialize (const BoxArray& grids,
		     int             n_grow,
                     int             n_comp,
                     const Geometry& geom);
    
    const BoxArray& equivBoxArray () const;

    void setVal (Real r);

    const DistributionMapping& DistributionMap () const;

    FArrayBox&       operator[] (const MFIter& mfi);
    const FArrayBox& operator[] (const MFIter& mfi) const;

    int nGrow () const;

    int nComp () const;

protected:

    MultiFab m_fabs;
    int      m_ngrow;
    bool     m_initialized;

private:
    //
    // Not defined -- but could be if it's needed.
    // 
    AuxBoundaryData& operator= (const AuxBoundaryData& rhs);
};

inline
void
AuxBoundaryData::setVal (Real r)
{
    BL_ASSERT(m_initialized);
    m_fabs.setVal(r);
}

inline
const DistributionMapping&
AuxBoundaryData::DistributionMap () const
{
    BL_ASSERT(m_initialized);
    return m_fabs.DistributionMap();
}

inline
int
AuxBoundaryData::nComp () const
{
    BL_ASSERT(m_initialized);
    return m_fabs.nComp();
}

inline
size_t
AuxBoundaryData::size () const
{
    BL_ASSERT(m_initialized);
    return m_fabs.size();
}

inline
int
AuxBoundaryData::nGrow () const
{
    BL_ASSERT(m_initialized);
    return m_ngrow;
}

inline
const BoxArray&
AuxBoundaryData::equivBoxArray () const
{
    BL_ASSERT(m_initialized);
    return m_fabs.boxArray();
}

inline
FArrayBox&
AuxBoundaryData::operator[] (const MFIter& mfi)
{
    BL_ASSERT(m_initialized);
    return m_fabs[mfi];
}

inline
const FArrayBox&
AuxBoundaryData::operator[] (const MFIter& mfi) const
{
    BL_ASSERT(m_initialized);
    return m_fabs[mfi];
}

#endif
