
#include <Periodicity.H>

std::vector<IntVect>
Periodicity::shiftIntVect () const
{
    std::vector<IntVect> r;

    int p[3] = { D_DECL(period[0], period[1], period[2]) };

    for (int i = -p[0]; i <= p[0]; i += p[0]) {
    for (int j = -p[1]; j <= p[1]; j += p[1]) {
    for (int k = -p[2]; k <= p[2]; k += p[2]) {
	r.push_back(IntVect(D_DECL(i,j,k)));
    }
    }
    }

    return r;
}
