
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>

using namespace amrex;

extern "C"
{
    void amrex_fi_write_multifab (const MultiFab* mf, const char* name)
    {
	VisMF::Write(*mf, std::string(name));
    }

    void amrex_fi_read_multifab (MultiFab* mf, const char* name)
    {
	BL_ASSERT(mf != nullptr);
	VisMF::Read(*mf, std::string(name));
    }
}
