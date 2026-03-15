
#ifndef AMREX_AuxBoundaryData_H_
#define AMREX_AuxBoundaryData_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

/**
 * \file AMReX_AuxBoundaryData.H
 * \brief Lightweight container for cached boundary data adjacent to AMR levels.
 */

namespace amrex {
class AuxBoundaryData
{
public:

    /**
     * \brief Construct an empty container; call initialize() before use.
     */
    AuxBoundaryData () noexcept = default;

    /**
     * \brief Construct and allocate storage for boundary data.
     *
     * \param ba    BoxArray describing the owning level.
     * \param n_grow Number of grow cells to capture around each box.
     * \param n_comp Number of components per FArrayBox.
     * \param geom  Geometry describing the level (used for periodicity).
     */
    AuxBoundaryData (const BoxArray& ba,
                     int             n_grow,
                     int             n_comp,
                     const Geometry& geom);

    //! Trivial destructor; MultiFab members handle their own cleanup.
    ~AuxBoundaryData () = default;

    //! Move constructor keeps the bookkeeping simple when containers relocate.
    AuxBoundaryData (AuxBoundaryData&& rhs) = default;
    //! Move assignment likewise transfers ownership of the cached MultiFab.
    AuxBoundaryData& operator= (AuxBoundaryData&& rhs) = default;

    /**
     * \brief Deep-copy constructor.
     *
     * \param rhs Source AuxBoundaryData whose storage is duplicated.
     */
    AuxBoundaryData (const AuxBoundaryData& rhs);
    //! Copy assignment is disabled; construct a new container instead.
    AuxBoundaryData& operator= (const AuxBoundaryData& rhs) = delete;

    /**
     * \brief Copy data from the auxiliary storage into \p destmf.
     *
     * \param destmf   Destination MultiFab.
     * \param src_comp Starting component inside this container.
     * \param dst_comp Destination component index.
     * \param num_comp Number of components to copy.
     */
    void copyTo (MultiFab& destmf,
                 int       src_comp,
                 int       dst_comp,
                 int       num_comp) const;

    /**
     * \brief Populate the auxiliary storage from \p srcmf.
     *
     * \param srcmf    Source MultiFab supplying data.
     * \param src_comp Starting component inside \p srcmf.
     * \param dst_comp Destination component inside this container.
     * \param num_comp Number of components to copy.
     * \param src_ng   Number of valid grow cells present in \p srcmf.
     */
    void copyFrom (const MultiFab& srcmf,
                   int             src_comp,
                   int             dst_comp,
                   int             num_comp,
                   int             src_ng = 0);

    //! Number of FArrayBoxes currently cached (0 if empty/uninitialized).
    size_t size () const noexcept
    {
        BL_ASSERT(!m_empty); BL_ASSERT(m_initialized); return m_fabs.size();
    }

    /**
     * \brief Copy data from another AuxBoundaryData \p src (component-wise).
     *
     * \param src      Source container.
     * \param src_comp Starting component inside \p src.
     * \param dst_comp Destination component inside this container.
     * \param num_comp Number of components to copy.
     */
    void copy (const AuxBoundaryData& src,
               int                    src_comp,
               int                    dst_comp,
               int                    num_comp);

    /**
     * \brief Allocate and configure storage.
     *
     * \param ba     BoxArray describing the owning level.
     * \param n_grow Number of grow cells to capture.
     * \param n_comp Number of components per face.
     * \param geom   Geometry describing the level.
     */
    void initialize (const BoxArray& ba,
                     int             n_grow,
                     int             n_comp,
                     const Geometry& geom);

    //! \return BoxArray describing the cached boundary regions.
    const BoxArray& equivBoxArray () const noexcept
    {
        BL_ASSERT(!m_empty); BL_ASSERT(m_initialized); return m_fabs.boxArray();
    }

    //! Fill all stored data with scalar value \p r (if initialized).
    void setVal (Real r) { BL_ASSERT(m_initialized); if (!m_empty) { m_fabs.setVal(r); } }

    //! Distribution mapping associated with the cached grids.
    const DistributionMapping& DistributionMap () const noexcept
    {
        BL_ASSERT(!m_empty); BL_ASSERT(m_initialized); return m_fabs.DistributionMap();
    }

    //! Mutable access to the boundary data at given \p mfi.
    FArrayBox&       operator[] (const MFIter& mfi) noexcept
    {
        BL_ASSERT(!m_empty); BL_ASSERT(m_initialized); return m_fabs[mfi];
    }
    //! Read-only access to the boundary data at given \p mfi.
    const FArrayBox& operator[] (const MFIter& mfi) const noexcept
    {
        BL_ASSERT(!m_empty); BL_ASSERT(m_initialized); return m_fabs[mfi];
    }

    //! Number of ghost cells captured when the container was initialized.
    int nGrow () const noexcept { BL_ASSERT(m_initialized); return m_ngrow; }

    //! Number of components stored per entry.
    int nComp () const noexcept
    {
        BL_ASSERT(!m_empty); BL_ASSERT(m_initialized); return m_fabs.nComp();
    }

    //! True if the container manages no data (either default constructed or after reset()).
    bool isEmpty () const noexcept { return m_empty; }

private:

    MultiFab m_fabs;
    int      m_ngrow{0};
    bool     m_empty{false};
    bool     m_initialized{false};
};

}

#endif
