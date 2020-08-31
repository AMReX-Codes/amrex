
#include <AMReX_DistributionMapping.H>
#include <AMReX_Print.H>

using namespace amrex;

extern "C" {

    void amrex_fi_new_distromap (DistributionMapping*& dm, const BoxArray* ba)
    {
	dm = new DistributionMapping(*ba);
    }

    void amrex_fi_new_distromap_from_pmap (DistributionMapping*& dm, const int* pmap, const int plen)
    {
        Vector<int> PMap(pmap,pmap+plen);
        dm = new DistributionMapping(std::move(PMap));
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

    void amrex_fi_distromap_get_pmap (const DistributionMapping* dm, int* pmap, const int plen)
    {
	Long dmsize = dm->size();
	AMREX_ASSERT(plen >= dmsize);
	for (int i = 0; i < dmsize && i < plen; ++i)
            pmap[i] = (*dm)[i];
    }

    void amrex_fi_print_distromap (const DistributionMapping* dm)
    {
	AllPrint() << *dm;
    }
}
