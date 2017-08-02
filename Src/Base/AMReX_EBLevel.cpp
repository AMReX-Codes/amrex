
#include <AMReX_EBLevel.H>
#include <type_traits>

namespace amrex {

EBLevel::EBLevel ()
{
}

EBLevel::~EBLevel ()
{
}

EBLevel::EBLevel (const BoxArray& ba, const DistributionMapping& dm, const Box& domain, const int ng)
    : EBLevelGrid(ba,dm,domain,ng),
      m_flags(ba, dm, 1, ng)
{
    static_assert(sizeof(EBCellFlag) == 4, "sizeof EBCellFlag != 4");
    static_assert(std::is_standard_layout<EBCellFlag>::value == true, "EBCellFlag is not pod");

    
}

}
