#include <AMReX_AnyDimUtil.H>

namespace amrex {

std::ostream&
operator<< (std::ostream& os, const BoxND& b)
{
    if (b.m_dim == 1) {
        os << "("
           << "(" << b.m_box.smallEnd(0) << ")" << " "
           << "(" << b.m_box.bigEnd(0) << ")" << " "
           << "(" << b.m_box.type(0) << ")"
           << ")";
    } else if (b.m_dim == 2) {
        os << "("
           << "(" << b.m_box.smallEnd(0) << "," << b.m_box.smallEnd(1) << ")" << " "
           << "(" << b.m_box.bigEnd(0) << "," << b.m_box.bigEnd(1) << ")" << " "
           << "(" << b.m_box.type(0) << "," << b.m_box.type(1) << ")"
           << ")";
    } else {
        os << b.m_box;
    }
    return os;
}

}


