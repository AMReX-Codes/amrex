
#include <AMReX_EB2.H>

namespace amrex { namespace EB2 {

IndexSpace* IndexSpace::m_pinstance = nullptr;
std::function<void()> IndexSpace::m_finalizer;

}}
