
#include <limits>
#include <AMReX_Periodicity.H>

namespace amrex {

std::vector<IntVect>
Periodicity::shiftIntVect () const
{
    std::vector<IntVect> r;

    int per[3] = {0,0,0};
    int jmp[3] = {1,1,1};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
	if (isPeriodic(i)) {
	    per[i] = jmp[i] = period[i];
	}
    }

    for (int i = -per[0]; i <= per[0]; i += jmp[0]) {
    for (int j = -per[1]; j <= per[1]; j += jmp[1]) {
    for (int k = -per[2]; k <= per[2]; k += jmp[2]) {
	r.push_back(IntVect(AMREX_D_DECL(i,j,k)));
    }
    }
    }

    return r;
}

Box
Periodicity::Domain () const noexcept
{
    Box pdomain;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
	if (isPeriodic(i)) {
	    pdomain.setSmall(i,0);
	    pdomain.setBig  (i,period[i]-1);
	} else {
	    pdomain.setSmall(i, std::numeric_limits<int>::min());
	    pdomain.setBig(i, std::numeric_limits<int>::max()-1); // so that it can be nodalized.
	}
    }
    return pdomain;
}

const Periodicity&
Periodicity::NonPeriodic () noexcept
{
    static const Periodicity np(IntVect(AMREX_D_DECL(0,0,0)));
    return np;
}

}
