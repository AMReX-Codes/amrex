#include <AMReX_EBCellFlag.H>
#include <AMReX_Reduce.H>

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
        std::map<Box,FabType>::iterator it;
#ifdef _OPENMP
#pragma omp critical (amrex_ebcellflagfab_gettype)
#endif
        it = m_typemap.find(bx);
        if (it != m_typemap.end())
        {
            return it->second;
        }
        else
        {
            auto const& flag = this->const_array();
            int nregular=0, nsingle=0, nmulti=0, ncovered=0;

            if (Gpu::inLaunchRegion())
            {
                ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum> reduce_op;
                ReduceData<int,int,int,int> reduce_data(reduce_op);
                using ReduceTuple = typename decltype(reduce_data)::Type;
                reduce_op.eval(bx, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                    int nr=0, ns=0, nm=0, nc=0;
                    auto f = flag(i,j,k);
                    if (f.isRegular()) {
                        ++nr;
                    } else if (f.isSingleValued()) {
                        ++ns;
                    } else if (f.isMultiValued()) {
                        ++nm;
                    } else {
                        ++nc;
                    }
                    return {nr, ns, nm, nc};
                });
                ReduceTuple hv = reduce_data.value();
                nregular = amrex::get<0>(hv);
                nsingle  = amrex::get<1>(hv);
                nmulti   = amrex::get<2>(hv);
                ncovered = amrex::get<3>(hv);
            }
            else
            {
                amrex::LoopOnCpu(bx,
                [=,&nregular,&nsingle,&nmulti,&ncovered] (int i, int j, int k) noexcept
                {
                    auto f = flag(i,j,k);
                    if (f.isRegular()) {
                        ++nregular;
                    } else if (f.isSingleValued()) {
                        ++nsingle;
                    } else if (f.isMultiValued()) {
                        ++nmulti;
                    } else {
                        ++ncovered;
                    }
                });
            }

            int ncells = bx.numPts();

            FabType t = FabType::undefined;
            if (nregular == ncells) {
                t = FabType::regular;
            } else if (ncovered == ncells) {
                t = FabType::covered;
            } else if (nmulti > 0) {
                t = FabType::multivalued;
            } else {
                t = FabType::singlevalued;
            }

#ifdef _OPENMP
#pragma omp critical (amrex_ebcellflagfab_gettype)
#endif
            m_typemap.insert({bx,t});

            return t;
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
