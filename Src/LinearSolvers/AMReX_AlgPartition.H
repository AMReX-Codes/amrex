#ifndef AMREX_ALG_PARTITION_H_
#define AMREX_ALG_PARTITION_H_
#include <AMReX_Config.H>

#include <AMReX_INT.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Vector.H>

#include <memory>

namespace amrex {

/**
 * \file AMReX_AlgPartition.H
 *
 * Lightweight helpers that describe how global sparse rows are distributed
 * across MPI ranks for algebra kernels.
 */

class AlgPartition
{
public:
    //! Construct an empty partition that must be defined later.
    AlgPartition ();
    //! Evenly partition \p global_size rows across all ranks.
    explicit AlgPartition (Long global_size);
    //! Adopt an explicit prefix array \p rows describing local row extents.
    explicit AlgPartition (Vector<Long> const& rows);
    //! Adopt an explicit prefix array \p rows, transferring ownership.
    explicit AlgPartition (Vector<Long>&& rows) noexcept;

    //! Recompute an even partition spanning \p global_size rows.
    void define (Long global_size);
    //! Replace the partition with the explicit prefix array \p rows.
    void define (Vector<Long> const& rows);
    //! Replace the partition, taking ownership of \p rows.
    void define (Vector<Long>&& rows);

    //! True if the partition contains no rows.
    [[nodiscard]] bool empty () const { return m_ref->m_row.empty(); }

    //! Starting global row index owned by rank \p i.
    [[nodiscard]] Long operator[] (int i) const { return m_ref->m_row[i]; }
    /**
     * \brief Total number of rows covered by the partition.
     */
    [[nodiscard]] Long numGlobalRows () const { return m_ref->m_row.back(); }
    //! Number of MPI ranks that own at least one row.
    [[nodiscard]] int numActiveProcs () const { return m_ref->m_n_active_procs; }

    //! Underlying prefix array describing row offsets (size nproc+1).
    [[nodiscard]] Vector<Long> const& dataVector () const { return m_ref->m_row; }

    //! Compare partitions for identical layouts.
    [[nodiscard]] bool operator== (AlgPartition const& rhs) const noexcept;
    //! Negation of operator==.
    [[nodiscard]] bool operator!= (AlgPartition const& rhs) const noexcept;

private:
    struct Ref
    {
        friend class AlgPartition;
        Ref () = default;
        explicit Ref (Long global_size);
        explicit Ref (Vector<Long> const& rows);
        explicit Ref (Vector<Long>&& rows);
        void define (Long global_size);
        void define (Vector<Long> const& rows);
        void define (Vector<Long>&& rows);
        void update_n_active_procs ();

        Vector<Long> m_row; // size: nprocs + 1
        int m_n_active_procs = 0;
    };

    std::shared_ptr<Ref> m_ref;
};

}

#endif
