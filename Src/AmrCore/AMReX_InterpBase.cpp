#include <AMReX_InterpBase.H>

namespace amrex {

InterpolaterBoxCoarsener::InterpolaterBoxCoarsener (InterpBase* mapper_, const IntVect& ratio_)
    : mapper(mapper_), ratio(ratio_)
{}

Box
InterpolaterBoxCoarsener::doit (const Box& fine) const
{
    return mapper->CoarseBox(fine, ratio);
}

BoxConverter*
InterpolaterBoxCoarsener::clone () const
{
    return new InterpolaterBoxCoarsener(mapper, ratio);
}

InterpolaterBoxCoarsener
InterpBase::BoxCoarsener (const IntVect& ratio)
{
    return InterpolaterBoxCoarsener(this, ratio);
}

Vector<int>
InterpBase::GetBCArray (const Vector<BCRec>& bcr)
{
    Vector<int> bc(2*AMREX_SPACEDIM*bcr.size());

    for (int n = 0; n < bcr.size(); n++)
    {
        const int* b_rec = bcr[n].vect();

        for (int m = 0; m < 2*AMREX_SPACEDIM; m++)
        {
            bc[2*AMREX_SPACEDIM*n + m] = b_rec[m];
        }
    }

    return bc;
}

}
