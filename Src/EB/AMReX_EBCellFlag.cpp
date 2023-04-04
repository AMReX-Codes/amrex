#include <AMReX_EBCellFlag.H>
#include <AMReX_Reduce.H>
#include <iostream>

namespace amrex {

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
        std::map<Box,NumCells>::iterator it;
#ifdef AMREX_USE_OMP
#pragma omp critical (amrex_ebcellflagfab_gettype)
#endif
        it = m_typemap.find(bx);
        if (it != m_typemap.end())
        {
            return it->second.type;
        }
        else
        {
            auto const& flag = this->const_array();
            auto const& t = countCells(flag, bx);

#ifdef AMREX_USE_OMP
#pragma omp critical (amrex_ebcellflagfab_gettype)
#endif
            m_typemap.insert({bx,t});

            return t.type;
        }
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
        std::map<Box,NumCells>::iterator it;
#ifdef AMREX_USE_OMP
#pragma omp critical (amrex_ebcellflagfab_gettype)
#endif
        it = m_typemap.find(bx);
        if (it != m_typemap.end())
        {
            return it->second.nregular;
        }
        else
        {
            auto const& flag = this->const_array();
            auto const& t = countCells(flag, bx);

#ifdef AMREX_USE_OMP
#pragma omp critical (amrex_ebcellflagfab_gettype)
#endif
            m_typemap.insert({bx,t});

            return t.nregular;
        }
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
        std::map<Box,NumCells>::iterator it;
#ifdef AMREX_USE_OMP
#pragma omp critical (amrex_ebcellflagfab_gettype)
#endif
        it = m_typemap.find(bx);
        if (it != m_typemap.end())
        {
            return it->second.nsingle;
        }
        else
        {
            auto const& flag = this->const_array();
            auto const& t = countCells(flag, bx);

#ifdef AMREX_USE_OMP
#pragma omp critical (amrex_ebcellflagfab_gettype)
#endif
            m_typemap.insert({bx,t});

            return t.nsingle;
        }
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
        std::map<Box,NumCells>::iterator it;
#ifdef AMREX_USE_OMP
#pragma omp critical (amrex_ebcellflagfab_gettype)
#endif
        it = m_typemap.find(bx);
        if (it != m_typemap.end())
        {
            return it->second.ncovered;
        }
        else
        {
            auto const& flag = this->const_array();
            auto const& t = countCells(flag, bx);

#ifdef AMREX_USE_OMP
#pragma omp critical (amrex_ebcellflagfab_gettype)
#endif
            m_typemap.insert({bx,t});

            return t.ncovered;
        }
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
