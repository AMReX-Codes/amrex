#include <AMReX_EBCellFlag.H>
#include <AMReX_Reduce.H>
#include <iostream>
#include <mutex>

namespace amrex {

namespace {
    /**
     * Guards EBCellFlagFab::m_typemap, the lazily filled cell-count cache that
     * getType() and the getNum*Cells() queries share.
     *
     * A real mutex rather than `omp critical`, which excludes only threads in
     * the same OpenMP contention group and is compiled out entirely in a build
     * without OpenMP -- these are const queries, so any host thread can be the
     * one that fills the cache.
     */
    std::mutex ebcellflag_typemap_mutex;
}

EBCellFlagFab::EBCellFlagFab (Arena* ar) noexcept
    : BaseFab<EBCellFlag>(ar)
{}

EBCellFlagFab::EBCellFlagFab (const Box& b, int n, Arena* ar)
    : BaseFab<EBCellFlag>(b,n,ar)
{}

EBCellFlagFab::EBCellFlagFab (const Box& b, int n, bool alloc, bool shared, Arena* ar)
    : BaseFab<EBCellFlag>(b,n,alloc,shared,ar)
{}

EBCellFlagFab::EBCellFlagFab (const EBCellFlagFab& rhs,
                              MakeType make_type, int scomp, int ncomp)
    : BaseFab<EBCellFlag>(rhs,make_type,scomp,ncomp)
{}

namespace {
EBCellFlagFab::NumCells countCells (Array4<EBCellFlag const> const& flag, const Box& bx) noexcept
{
    int nregular=0, nsingle=0, nmulti=0;
    int ncells = static_cast<int>(bx.numPts());
    AMREX_ASSERT(bx.numPts() <= static_cast<Long>(std::numeric_limits<int>::max()));

    if (Gpu::inLaunchRegion())
    {
        ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpSum> reduce_op;
        ReduceData<int,int,int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            int nr=0, ns=0, nm=0;
            auto f = flag(i,j,k);
            if (f.isRegular()) {
                ++nr;
            } else if (f.isSingleValued()) {
                ++ns;
            } else if (f.isMultiValued()) {
                ++nm;
            }
            return {nr, ns, nm};
        });
        ReduceTuple hv = reduce_data.value(reduce_op);
        nregular = amrex::get<0>(hv);
        nsingle  = amrex::get<1>(hv);
        nmulti   = amrex::get<2>(hv);
    }
    else
    {
        amrex::LoopOnCpu(bx,
        [=,&nregular,&nsingle,&nmulti] (int i, int j, int k) noexcept
        {
            auto f = flag(i,j,k);
            if (f.isRegular()) {
                ++nregular;
            } else if (f.isSingleValued()) {
                ++nsingle;
            } else if (f.isMultiValued()) {
                ++nmulti;
            }
        });
    }

    int ncovered = ncells - nregular - nsingle - nmulti;

    EBCellFlagFab::NumCells r;
    if (nregular == ncells) {
        r.type = FabType::regular;
    } else if (ncovered == ncells) {
        r.type = FabType::covered;
    } else if (nmulti > 0) {
        r.type = FabType::multivalued;
    } else {
        r.type = FabType::singlevalued;
    }
    r.nregular = nregular;
    r.nsingle = nsingle;
    r.nmulti = nmulti;
    r.ncovered = ncovered;
    return r;
}
}

FabType
EBCellFlagFab::getType (const Box& bx_in) const noexcept
{
    FabType thistype = getType();

    if (thistype == FabType::regular)
    {
        return FabType::regular;
    }
    else if (thistype == FabType::covered)
    {
        return FabType::covered;
    }
    else
    {
        const Box& bx = amrex::enclosedCells(bx_in);
        {
            std::scoped_lock lock(ebcellflag_typemap_mutex);
            auto it = m_typemap.find(bx);
            if (it != m_typemap.end())
            {
                return it->second.type;
            }
        }

        auto const& flag = this->const_array();
        auto const& t = countCells(flag, bx);

        {
            std::scoped_lock lock(ebcellflag_typemap_mutex);
            m_typemap.insert({bx,t});
        }

        return t.type;
    }
}

void
EBCellFlagFab::resetType (int ng)
{
    {
        // getType() below takes this lock too, so it has to be released first
        // -- ebcellflag_typemap_mutex is a plain std::mutex.
        std::scoped_lock lock(ebcellflag_typemap_mutex);
        this->setType(FabType::undefined);
        m_typemap.clear();
    }

    Box const& bx = this->box();
    auto typ = this->getType(bx);
    this->setType(typ);
    for (int nshrink = 1; nshrink < ng; ++nshrink) {
        const Box& b = amrex::grow(bx,-nshrink);
        this->getType(b);
    }
}

int
EBCellFlagFab::getNumRegularCells (const Box& bx_in) const noexcept
{
    FabType thistype = getType();

    const Box& bx = amrex::enclosedCells(bx_in);

    if (thistype == FabType::regular)
    {
        return static_cast<int>(bx.numPts());
    }
    else if (thistype == FabType::covered)
    {
        return 0;
    }
    else
    {
        {
            std::scoped_lock lock(ebcellflag_typemap_mutex);
            auto it = m_typemap.find(bx);
            if (it != m_typemap.end())
            {
                return it->second.nregular;
            }
        }

        auto const& flag = this->const_array();
        auto const& t = countCells(flag, bx);

        {
            std::scoped_lock lock(ebcellflag_typemap_mutex);
            m_typemap.insert({bx,t});
        }

        return t.nregular;
    }
}

int
EBCellFlagFab::getNumCutCells (const Box& bx_in) const noexcept
{
    FabType thistype = getType();

    const Box& bx = amrex::enclosedCells(bx_in);

    if (thistype == FabType::regular ||
        thistype == FabType::covered)
    {
        return 0;
    }
    else
    {
        {
            std::scoped_lock lock(ebcellflag_typemap_mutex);
            auto it = m_typemap.find(bx);
            if (it != m_typemap.end())
            {
                return it->second.nsingle;
            }
        }

        auto const& flag = this->const_array();
        auto const& t = countCells(flag, bx);

        {
            std::scoped_lock lock(ebcellflag_typemap_mutex);
            m_typemap.insert({bx,t});
        }

        return t.nsingle;
    }
}

int
EBCellFlagFab::getNumCoveredCells (const Box& bx_in) const noexcept
{
    FabType thistype = getType();

    const Box& bx = amrex::enclosedCells(bx_in);

    if (thistype == FabType::regular)
    {
        return 0;
    }
    else if (thistype == FabType::covered)
    {
        return static_cast<int>(bx.numPts());
    }
    else
    {
        {
            std::scoped_lock lock(ebcellflag_typemap_mutex);
            auto it = m_typemap.find(bx);
            if (it != m_typemap.end())
            {
                return it->second.ncovered;
            }
        }

        auto const& flag = this->const_array();
        auto const& t = countCells(flag, bx);

        {
            std::scoped_lock lock(ebcellflag_typemap_mutex);
            m_typemap.insert({bx,t});
        }

        return t.ncovered;
    }
}

std::ostream&
operator<< (std::ostream& os, const EBCellFlag& flag)
{
    std::ios_base::fmtflags old_fmt = os.flags();
    os << std::hex << flag.getValue() << ":" << std::dec;

    if (flag.isRegular()) {
        os << "R";
    } else if (flag.isSingleValued()) {
        os << "S";
    } else if (flag.isCovered()) {
        os << "C";
    } else {
        os << "M";
    }

#if (AMREX_SPACEDIM == 3)
    for (int k = -1; k <= 1; ++k) {
#endif
        for (int j = -1; j <= 1; ++j) {
            for (int i = -1; i <= 1; ++i) {
                os << static_cast<int>(flag.isConnected(IntVect{AMREX_D_DECL(i,j,k)}));
            }
        }
#if (AMREX_SPACEDIM == 3)
    }
#endif

    os.flags(old_fmt);

    return os;
}

}
