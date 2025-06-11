#include <AMReX_AlgPartition.H>

namespace amrex {

AlgPartition::AlgPartition ()
    : m_ref(std::make_shared<Ref>())
{}

AlgPartition::AlgPartition (Long global_size)
    : m_ref(std::make_shared<Ref>(global_size))
{}

AlgPartition::AlgPartition (Vector<Long> const& rows)
    : m_ref(std::make_shared<Ref>(rows))
{}

AlgPartition::AlgPartition (Vector<Long>&& rows) noexcept
    : m_ref(std::make_shared<Ref>(std::move(rows)))
{}

void AlgPartition::define (Long global_size)
{
    m_ref->define(global_size);
}

void AlgPartition::define (Vector<Long> const& rows)
{
    m_ref->define(rows);
}

void AlgPartition::define (Vector<Long>&& rows)
{
    m_ref->define(std::move(rows));
}

bool AlgPartition::operator== (AlgPartition const& rhs) const noexcept
{
    return m_ref == rhs.m_ref || m_ref->m_row == rhs.m_ref->m_row;
}

bool AlgPartition::operator!= (AlgPartition const& rhs) const noexcept
{
    return !operator==(rhs);
}

AlgPartition::Ref::Ref (Long global_size)
{
    define(global_size);
}

AlgPartition::Ref::Ref (Vector<Long> const& rows)
    : m_row(rows)
{
    update_n_active_procs();
}

AlgPartition::Ref::Ref (Vector<Long>&& rows)
    : m_row(std::move(rows))
{
    update_n_active_procs();
}

void AlgPartition::Ref::define (Long global_size)
{
    auto nprocs = Long(ParallelDescriptor::NProcs());
    Long sz = global_size / nprocs;
    Long extra = global_size - sz*nprocs;
    m_row.resize(nprocs+1);
    for (Long i = 0; i < nprocs; ++i) {
        if (i < extra) {
            m_row[i] = i*(sz+1);
        } else {
            m_row[i] = i*sz + extra;
        }
    }
    m_row[nprocs] = global_size;

    update_n_active_procs();
}

void AlgPartition::Ref::define (Vector<Long> const& rows)
{
    m_row = rows;
    update_n_active_procs();
}

void AlgPartition::Ref::define (Vector<Long>&& rows)
{
    m_row = std::move(rows);
    update_n_active_procs();
}

void AlgPartition::Ref::update_n_active_procs ()
{
    AMREX_ASSERT(m_row.size() == ParallelDescriptor::NProcs()+1);
    m_n_active_procs = 0;
    for (int i = 0, N = int(m_row.size())-1; i < N; ++i) {
        if (m_row[i] < m_row[i+1]) { ++m_n_active_procs; }
    }
}

}
