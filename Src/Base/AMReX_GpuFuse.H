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
typedef void (*Lambda4DLauncher)(char*,int,int,int,int);
typedef void (*LambdaDtor)(char*);

struct FuseHelper {
    union LambdaLauncher {
        Lambda3DLauncher L3D;
        Lambda4DLauncher L4D;
    };
    // Device function pointer to wrapper around Lambda's operator()
    LambdaLauncher m_fp;
    std::size_t m_offset; // This tells us where the lambda object is.
    Box m_bx;
    int m_ncomp;
};

template <typename Lambda>
AMREX_GPU_DEVICE
void
LaunchLambda3D (char* buf, int i, int j, int k)
{
    (*reinterpret_cast<Lambda*>(buf))(i,j,k);
}

template <typename Lambda>
AMREX_GPU_DEVICE
void
LaunchLambda4D (char* buf, int i, int j, int k, int n)
{
    (*reinterpret_cast<Lambda*>(buf))(i,j,k,n);
}

template <typename Lambda>
void
DestroyLambda (char* buf)
{
    reinterpret_cast<Lambda*>(buf)->~Lambda();
}

template <typename Lambda>
std::size_t
PutLambda (Lambda const& f, char* buf)
{
    new (buf) Lambda(f);
    return Arena::align(sizeof(Lambda));
}

template <typename Lambda>
void
PutLambda3DLauncher (FuseHelper* helper)
{
    static Lambda3DLauncher fp = nullptr;
    if (fp == nullptr) {
        launch_global<<<1,1>>>([=] AMREX_GPU_DEVICE
        {
            helper->m_fp.L3D = &LaunchLambda3D<Lambda>;
        });
        Gpu::synchronize();
        fp = helper->m_fp.L3D;
    } else {
        helper->m_fp.L3D = fp;
    }
}

template <typename Lambda>
void
PutLambda4DLauncher (FuseHelper* helper)
{
    static Lambda4DLauncher fp = nullptr;
    if (fp == nullptr) {
        launch_global<<<1,1>>>([=] AMREX_GPU_DEVICE
        {
            helper->m_fp.L4D = &LaunchLambda4D<Lambda>;
        });
        Gpu::synchronize();
        fp = helper->m_fp.L4D;
    } else {
        helper->m_fp.L4D = fp;
    }
}

template <typename Lambda>
void
PutLambdaDtor (LambdaDtor& lambda_dtor)
{
    lambda_dtor = &DestroyLambda<Lambda>;
}

class Fuser
{
public:

    Fuser ();
    ~Fuser ();

    template <typename F>
    void
    Register (Box const& bx, F&& f)
    {
        if (bx.isEmpty()) return;

        if (m_nlambdas >= m_nhelpers_buf) {
            resize_helper_buf();
        }

        int ikernel = m_nlambdas++;
        new (m_helper_buf+ikernel) FuseHelper;
        m_helper_buf[ikernel].m_offset = m_nbytes_used_lambda_buf;
        m_helper_buf[ikernel].m_bx = bx;
        m_helper_buf[ikernel].m_ncomp = 0;

        using Lambda = typename std::decay<F>::type;
        static_assert(std::alignment_of<Lambda>::value <= Arena::align_size,
                      "GpuFuse: alignment error");
        if (m_nbytes_used_lambda_buf + Arena::align(sizeof(Lambda)) > m_nbytes_lambda_buf) {
            resize_lambda_buf();
        }
        char* p_lambda = m_lambda_buf + m_nbytes_used_lambda_buf;
        std::size_t sz_lambda =  PutLambda<Lambda>(f, p_lambda);
        m_nbytes_used_lambda_buf += sz_lambda;
        AMREX_ASSERT(m_nbytes_used_lambda_buf <= m_nbytes_lambda_buf);

        PutLambda3DLauncher<Lambda>(m_helper_buf+ikernel);
        m_dtor_buf.push_back(nullptr);
        PutLambdaDtor<Lambda>(m_dtor_buf.back());
    }

    template <typename F>
    void
    Register (Box const& bx, int ncomp, F&& f)
    {
        if (bx.isEmpty()) return;

        if (m_nlambdas >= m_nhelpers_buf) {
            resize_helper_buf();
        }

        int ikernel = m_nlambdas++;
        new (m_helper_buf+ikernel) FuseHelper;
        m_helper_buf[ikernel].m_offset = m_nbytes_used_lambda_buf;
        m_helper_buf[ikernel].m_bx = bx;
        m_helper_buf[ikernel].m_ncomp = ncomp;

        using Lambda = typename std::decay<F>::type;
        static_assert(std::alignment_of<Lambda>::value <= Arena::align_size,
                      "GpuFuse: alignment error");
        if (m_nbytes_used_lambda_buf + Arena::align(sizeof(Lambda)) > m_nbytes_lambda_buf) {
            resize_lambda_buf();
        }
        char* p_lambda = m_lambda_buf + m_nbytes_used_lambda_buf;
        std::size_t sz_lambda =  PutLambda<Lambda>(f, p_lambda);
        m_nbytes_used_lambda_buf += sz_lambda;
        AMREX_ASSERT(m_nbytes_used_lambda_buf <= m_nbytes_lambda_buf);

        PutLambda4DLauncher<Lambda>(m_helper_buf+ikernel);
        m_dtor_buf.push_back(nullptr);
        PutLambdaDtor<Lambda>(m_dtor_buf.back());
    }

    void Launch ();

private:

    char* m_lambda_buf = nullptr;
    std::size_t m_nbytes_used_lambda_buf = 0;
    std::size_t m_nbytes_lambda_buf = 2*1024*1024;

    FuseHelper* m_helper_buf = nullptr;
    std::size_t m_nhelpers_buf = 1024;

    Vector<LambdaDtor> m_dtor_buf;

    int m_nlambdas = 0;

    void resize_lambda_buf ();
    void resize_helper_buf ();
};

#endif

}}
