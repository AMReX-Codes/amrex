#include <Geometry.H>

extern "C"
{
    void fi_new_geometry (Geometry*& geom, int lo[3], int hi[3])
    {
	Box domain(IntVect(D_DECL(lo[0],lo[1],lo[2])),
		   IntVect(D_DECL(hi[0],hi[1],hi[2])));
	geom = new Geometry(domain);
    }

    void fi_delete_geometry (Geometry* geom)
    {
	delete geom;
    }

    void fi_geometry_get_pmask (Geometry* geom, int is_per[3])
    {
	for (int i = 0; i < BL_SPACEDIM; ++i)
	    is_per[i] = geom->isPeriodic(i);
    }

    void fi_geometry_get_probdomain (Geometry* geom, double problo[3], double probhi[3])
    {
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    problo[i] = geom->ProbLo(i);
	    probhi[i] = geom->ProbHi(i);
	}
    }
}
