#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_CudaDevice.H>

namespace amrex {

struct CopyMemory
{
    Array4<void const> src{nullptr,{},{}};  // think of int as void
    Array4<void> dst{nullptr,{},{}};
    int scomp = 0;
    int ncomp = 0;

    template <class T>
    AMREX_GPU_HOST_DEVICE
    Array4<T const> getSrc () { return ToArray4<T const>(src); }

    template <class T>
    AMREX_GPU_HOST_DEVICE
    Array4<T> getDst () { return ToArray4<T>(dst); }
};

template <typename T, typename U>
CopyMemory
makeCopyMemory (Array4<T> const& src, Array4<U> const& dst, int scomp, int ncomp)
{
    return { ToArray4<void const>(src), ToArray4<void>(dst), scomp, ncomp };
}

template <typename T>
struct CudaGraph
{
    cudaGraphExec_t m_graph;
    Vector<T> m_parms;
    T* m_parms_d = nullptr;

    CudaGraph(int num)
        : m_parms(num)
    {
        cudaMalloc(&m_parms_d, sizeof(T)*m_parms.size());
    }
    ~CudaGraph()
    {
        cudaFree(m_parms_d);
    }
    void setGraph(cudaGraphExec_t const& graph)
    { 
        m_graph = graph;
    }
    void setParams(int idx, T const& a_parm)
    {
        m_parms[idx] = a_parm;
    }
    T* getDevicePtr (int idx) {
        return m_parms_d + idx;
    }
    void executeGraph()
    {
        static_assert(std::is_trivially_copyable<T>::value, "CudaGraph's T must be trivially copyable");
        cudaMemcpy(m_parms_d, m_parms.data(), std::size_t(sizeof(T)*m_parms.size()),
                   cudaMemcpyHostToDevice);

        Cuda::Device::setStreamIndex(0);

        cudaGraphLaunch(m_graph, amrex::Cuda::Device::cudaStream());

        Cuda::Device::resetStreamIndex();
        Cuda::Device::synchronize();
    }
};

}
