#include <AMReX_Sundials.H>

namespace amrex {
namespace sundials {

namespace {
    Vector<int> initialized;
    Vector<::sundials::Context*> the_sundials_context;
}

void Initialize(int nthreads)
{
    amrex::Print() << "Initializing SUNDIALS with " << nthreads << " threads...\n";

    // Initalize the sundials context
    if (initialized.size() == 0) {
        initialized.resize(nthreads);
        std::fill(initialized.begin(), initialized.end(), 0);
        the_sundials_context.resize(nthreads);
        std::fill(the_sundials_context.begin(), the_sundials_context.end(), nullptr);
    }
    for (int i = 0; i < nthreads; i++) {
        if (initialized[i]) continue;
        initialized[i] = 1;
        BL_ASSERT(the_sundials_context[i] == nullptr);
        the_sundials_context[i] = new ::sundials::Context();
    }

    // Initialize the memory helper
    MemoryHelper::Initialize(nthreads);

    amrex::Print() << "SUNDIALS initialized.\n";
}

void Finalize()
{
    // Cleanup the memory helpers
    MemoryHelper::Finalize();

    // Clean up the sundials contexts
    for (int i = 0; i < initialized.size(); i++) {
        initialized[i] = 0;
        delete the_sundials_context[i];
        the_sundials_context[i] = nullptr;
    }
}

::sundials::Context* The_Sundials_Context(int i)
{
    BL_ASSERT(the_sundials_context[i] != nullptr);
    return the_sundials_context[i];
}

}//sundials
}//amrex
