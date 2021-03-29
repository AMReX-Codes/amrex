
#include <AMReX_Dim3.H>
#include <iostream>

namespace amrex {

std::ostream& operator<< (std::ostream& os, const Dim3& d)
{
    os << '(' << d.x << ',' << d.y << ',' << d.z << ')';
    return os;
}

}
