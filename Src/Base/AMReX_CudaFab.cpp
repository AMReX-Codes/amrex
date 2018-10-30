
#include <AMReX_CudaFab.H>
#include <AMReX_FArrayBox.H>

namespace amrex {
namespace Cuda {

DeviceFab::DeviceFab ()
{
    static_assert(AMREX_IS_TRIVIALLY_COPYABLE(BaseFabData<Real>),
                  "BaseFabData must be trivially copyable");
//    m_fab = new FArrayBox();
}

DeviceFab::DeviceFab (Box const& bx, int ncomp)
{
}

DeviceFab::~DeviceFab ()
{

}

}
}
