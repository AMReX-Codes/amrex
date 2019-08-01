
#include <AMReX.H>
#include <AMReX_Scan.H>

using namespace amrex;

void main_main();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    main_main();
    amrex::Finalize();
}

void main_main ()
{
    const int N = 512*1024*1024;
    Vector<int> h_in(N);
    for (auto& x: h_in) {
        x = static_cast<int>((Random()-0.5)*10.);
    }
    Vector<int> h_inclusive_thrust(N);
    Vector<int> h_inclusive_amrex(N);
    Vector<int> h_exclusive_thrust(N);
    Vector<int> h_exclusive_amrex(N);

    std::size_t nbytes = h_in.size()*sizeof(int);
    int* d_in = (int*)The_Device_Arena()->alloc(nbytes);
    int* d_out = (int*)The_Device_Arena()->alloc(nbytes);

    amrex::Gpu::htod_memcpy(d_in, h_in.data(), nbytes);

    // warm up
    Scan::InclusiveSum(N, d_in, d_out);
    Scan::ExclusiveSum(N, d_in, d_out);
    Gpu::synchronize();
    thrust::inclusive_scan(thrust::cuda::par(Cuda::The_ThrustCachedAllocator()),
                           d_in, d_in+N, d_out);
    thrust::exclusive_scan(thrust::cuda::par(Cuda::The_ThrustCachedAllocator()),
                           d_in, d_in+N, d_out);

    {
        BL_PROFILE("amrex::inclusive_scan");
        Scan::InclusiveSum(N, d_in, d_out);
        Gpu::synchronize();
    }
    amrex::Gpu::dtoh_memcpy(h_inclusive_amrex.data(), d_out, nbytes);

    {
        BL_PROFILE("amrex::exclusive_scan");
        Scan::ExclusiveSum(N, d_in, d_out);
        Gpu::synchronize();
    }
    amrex::Gpu::dtoh_memcpy(h_exclusive_amrex.data(), d_out, nbytes);

    {
        BL_PROFILE("thrust::inclusive_scan");
        thrust::inclusive_scan(thrust::cuda::par(Cuda::The_ThrustCachedAllocator()),
                               d_in, d_in+N, d_out);
    }
    amrex::Gpu::dtoh_memcpy(h_inclusive_thrust.data(), d_out, nbytes);

    {
        BL_PROFILE("thrust::exclusive_scan");
        thrust::exclusive_scan(thrust::cuda::par(Cuda::The_ThrustCachedAllocator()),
                               d_in, d_in+N, d_out);
    }
    amrex::Gpu::dtoh_memcpy(h_exclusive_thrust.data(), d_out, nbytes);

    for (int i = 0; i < N; ++i) {
        if (h_inclusive_thrust[i] != h_inclusive_amrex[i]) {
            amrex::Print() << "i = " << i << ", thrust = " << h_inclusive_thrust[i]
                           << ", amrex = " << h_inclusive_amrex[i] << "\n"; 
            amrex::Abort("inclusive scan failed");
        }
        if (h_exclusive_thrust[i] != h_exclusive_amrex[i]) {
            amrex::Print() << "i = " << i << ", thrust = " << h_exclusive_thrust[i]
                           << ", amrex = " << h_exclusive_amrex[i] << "\n"; 
            amrex::Abort("exclusive scan failed");
        }
    }

    The_Device_Arena()->free(d_in);
    The_Device_Arena()->free(d_out);
}
