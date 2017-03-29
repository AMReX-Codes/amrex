
#include <WarpXWrappers.h>

extern void inject_particles ();

int main(int argc, char* argv[])
{
    amrex_init(argc, argv);

    static_assert(BL_SPACEDIM == 3, "3D only");

    warpx_init();

    for (int i = 1; i <= 10; ++i)
    {
	warpx_evolve(i);
	inject_particles();
    }

    warpx_finalize();

    amrex_finalize(1);
}
