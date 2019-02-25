
#include <AMReX.H>

extern "C"
{
    void amrex_fi_error (const char* message)
    {
	amrex::Error(message);
    }

    void amrex_fi_abort (const char* message)
    {
	amrex::Abort(message);
    }

    void amrex_fi_warning (const char* message)
    {
	amrex::Warning(message);
    }
}
