
#include <AMReX_DistributionMapping.H>

using namespace amrex;

extern "C" {

    void amrex_fi_new_distromap (DistributionMapping*& dm, const BoxArray* ba)
    {
	dm = new DistributionMapping(*ba);
    }

    void amrex_fi_delete_distromap (DistributionMapping* dm)
    {
	delete dm;
    }

    void amrex_fi_clone_distromap (DistributionMapping*& dmo, const DistributionMapping* dmi)
    {
	delete dmo;
	dmo = new DistributionMapping(*dmi);
    }

}
