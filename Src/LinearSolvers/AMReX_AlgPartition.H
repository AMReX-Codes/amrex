#ifndef AMREX_ALG_PARTITION_H_
#define AMREX_ALG_PARTITION_H_
#include <AMReX_Config.H>

#include <AMReX_INT.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Vector.H>

#include <memory>

namespace amrex {

class AlgPartition
{
public:
    AlgPartition ();
    explicit AlgPartition (Long global_size);
    explicit AlgPartition (Vector<Long> const& rows);
    explicit AlgPartition (Vector<Long>&& rows) noexcept;

    void define (Long global_size);
    void define (Vector<Long> const& rows);
    void define (Vector<Long>&& rows);

    [[nodiscard]] bool empty () const { return m_ref->m_row.empty(); }

    [[nodiscard]] Long operator[] (int i) const { return m_ref->m_row[i]; }
    [[nodiscard]] Long numGlobalRows () const { return m_ref->m_row.back(); }
    [[nodiscard]] int numActiveProcs () const { return m_ref->m_n_active_procs; }

    [[nodiscard]] Vector<Long> const& dataVector () const { return m_ref->m_row; }

    [[nodiscard]] bool operator== (AlgPartition const& rhs) const noexcept;
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
