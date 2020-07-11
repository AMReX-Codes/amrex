#include <AMReX_Arena.H>
#include <AMReX_Vector.H>
#include <AMReX_Box.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunchGlobal.H>
#include <cstring>

namespace amrex {
namespace Gpu {

#ifdef AMREX_USE_CUDA

typedef void (*Lambda3DLauncher)(char*,int,int,int);

struct Fuse3DHelper {
    Lambda3DLauncher m_fp = nullptr; // Device function pointer to wrapper around Lambda's operator()
    std::size_t m_offset; // This tells us where the lambda object is.
    Box m_bx;
};

template <typename Lambda>
AMREX_GPU_DEVICE
void
LaunchLambda (char* buf, int i, int j, int k)
{
    (*reinterpret_cast<Lambda*>(buf))(i,j,k);
}

template <typename Lambda>
std::size_t
PutLambda (Lambda const& f, char* buf) {
    new (buf) Lambda(f);
    return Arena::align(sizeof(Lambda));
}

template <typename Lambda>
void
PutLambda3DLauncher (Fuse3DHelper* helper) {
    static Lambda3DLauncher fp = nullptr;
    if (fp == nullptr) {
        launch_global<<<1,1>>>([=] AMREX_GPU_DEVICE
        {
            helper->m_fp = &LaunchLambda<Lambda>;
        });
        Gpu::synchronize();
        fp = helper->m_fp;
    } else {
        helper->m_fp = fp;
    }
}

class Fuser3D
{
public:

    Fuser3D ();
    ~Fuser3D ();

    template <typename F>
    void
    Register (Box const& bx, F&& f)
    {
        if (bx.isEmpty()) return;

        if (m_nlambdas >= m_nhelpers_buf) {
            m_nhelpers_buf += m_nhelpers_buf/2;
            auto p = (Fuse3DHelper*)The_Pinned_Arena()->alloc
                (m_nhelpers_buf*sizeof(Fuse3DHelper));
            for (int i = 0; i < m_nlambdas; ++i) {
                new (p+i) Fuse3DHelper (m_helper_buf[i]);
                (m_helper_buf+i)->~Fuse3DHelper();
            }
            The_Pinned_Arena()->free(m_helper_buf);
            m_helper_buf = p;
        }

        int ikernel = m_nlambdas++;
        new (m_helper_buf+ikernel) Fuse3DHelper;
        m_helper_buf[ikernel].m_offset = m_nbytes_used_lambda_buf;
        m_helper_buf[ikernel].m_bx = bx;

        using Lambda = typename std::decay<F>::type;
        static_assert(std::alignment_of<Lambda>::value <= Arena::align_size,
                      "GpuFuse: alignment error");
        if (m_nbytes_used_lambda_buf + Arena::align(sizeof(Lambda)) > m_nbytes_lambda_buf) {
            m_nbytes_lambda_buf += m_nbytes_lambda_buf/2;
            auto p = (char*)The_Pinned_Arena()->alloc(m_nbytes_lambda_buf);
            std::memcpy(p, m_lambda_buf, m_nbytes_used_lambda_buf);
            The_Pinned_Arena()->free(m_lambda_buf);
            m_lambda_buf = p;
        }
        char* p_lambda = m_lambda_buf + m_nbytes_used_lambda_buf;
        std::size_t sz_lambda =  PutLambda<Lambda>(f, p_lambda);
        m_nbytes_used_lambda_buf += sz_lambda;

        // how about ~Lambda?????

        PutLambda3DLauncher<Lambda>(m_helper_buf+ikernel);
    }

    void Launch ();

private:

    char* m_lambda_buf = nullptr;
    std::size_t m_nbytes_used_lambda_buf = 0;
    std::size_t m_nbytes_lambda_buf = 2*1024*1024;

    Fuse3DHelper* m_helper_buf = nullptr;
    std::size_t m_nhelpers_buf = 1024;

    int m_nlambdas = 0;
};

#endif

}}
