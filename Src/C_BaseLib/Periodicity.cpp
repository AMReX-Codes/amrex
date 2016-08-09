
#include <Periodicity.H>

std::vector<IntVect>
Periodicity::shiftIntVect () const
{
    std::vector<IntVect> r;

    int per[3] = {0,0,0};
    int jmp[3] = {1,1,1};

    for (int i = 0; i < BL_SPACEDIM; ++i) {
	if (isPeriodic(i)) {
	    per[i] = jmp[i] = period[i];
	}
    }

    for (int i = -per[0]; i <= per[0]; i += jmp[0]) {
    for (int j = -per[1]; j <= per[1]; j += jmp[1]) {
    for (int k = -per[2]; k <= per[2]; k += jmp[2]) {
	r.push_back(IntVect(D_DECL(i,j,k)));
    }
    }
    }

    return r;
}
