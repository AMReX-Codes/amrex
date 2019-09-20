#include <AMReX_LO_BCTYPES.H>
#include <AMReX.H>
#include <iostream>

namespace amrex {

std::ostream& operator<< (std::ostream& os, const LinOpBCType& t)
{
    switch (t) {
        case LinOpBCType::interior:
        {
            os << "interior";
            break;
        }
        case LinOpBCType::Dirichlet:
        {
            os << "Dirichlet";
            break;
        }
        case LinOpBCType::Neumann:
        {
            os << "Neumann";
            break;
        }
        case LinOpBCType::reflect_odd:
        {
            os << "reflect_odd";
            break;
        }
        case LinOpBCType::Marshak:
        {
            os << "Marshak";
            break;
        }
        case LinOpBCType::SanchezPomraning:
        {
            os << "SanchezPomraning";
            break;
        }
        case LinOpBCType::inflow:
        {
            os << "inflow";
            break;
        }
        case LinOpBCType::inhomogNeumann:
        {
            os << "inhomogeneous Neumann";
            break;
        }
        case LinOpBCType::Periodic:
        {
            os << "Periodic";
            break;
        }
        default:
        {
            os << "bogus";
        }
    };
    return os;
}

}
