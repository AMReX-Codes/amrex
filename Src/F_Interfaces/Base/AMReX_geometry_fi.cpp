#include <AMReX_Geometry.H>

using namespace amrex;

extern "C"
{
    void amrex_fi_new_geometry (Geometry*& geom, int lo[3], int hi[3])
    {
	const Box domain{IntVect{AMREX_D_DECL(lo[0],lo[1],lo[2])},
                         IntVect{AMREX_D_DECL(hi[0],hi[1],hi[2])}};
	geom = new Geometry(domain);
    }

    void amrex_fi_delete_geometry (Geometry* geom)
    {
	delete geom;
    }

    void amrex_fi_geometry_get_pmask (int is_per[3])
    {
        Geometry* gg = AMReX::top()->getDefaultGeometry();
	for (int i = 0; i < BL_SPACEDIM; ++i)
	    is_per[i] = gg->isPeriodic(i);
    }

    void amrex_fi_geometry_get_probdomain (Real problo[3], Real probhi[3])
    {
        Geometry* gg = AMReX::top()->getDefaultGeometry();
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    problo[i] = gg->ProbLo(i);
	    probhi[i] = gg->ProbHi(i);
	}
    }

    void amrex_fi_geometry_get_intdomain (const Geometry* geom, int lo[3], int hi[3])
    {
	const Box& bx = geom->Domain();
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    lo[i] = bx.smallEnd(i);
	    hi[i] = bx.bigEnd(i);
	}
    }
}
