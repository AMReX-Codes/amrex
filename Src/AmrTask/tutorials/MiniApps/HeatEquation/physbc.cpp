
#include <AMReX_BCUtil.H>

using namespace amrex;

void fill_physbc (MultiFab& phi, const Geometry& geom)
{
    if (Geometry::isAllPeriodic()) return;

    // Set up BC; see Src/Base/AMReX_BC_TYPES.H for supported types
    Vector<BCRec> bc(phi.nComp());
    for (int n = 0; n < phi.nComp(); ++n)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            if (Geometry::isPeriodic(idim))
            {
                bc[n].setLo(idim, BCType::int_dir); // interior
                bc[n].setHi(idim, BCType::int_dir);
            }
            else
            {
                bc[n].setLo(idim, BCType::foextrap); // first-order extrapolation.
                bc[n].setHi(idim, BCType::foextrap);
            }
        }
    }

    amrex::FillDomainBoundary(phi, geom, bc);
}
