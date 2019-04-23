#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_CudaDevice.H>

namespace amrex {

struct CopyMemory
{
    void* src = nullptr;
    void* dst = nullptr;
    Dim3 src_begin = {0,0,0};
    Dim3 src_end = {0,0,0};
    Dim3 dst_begin = {0,0,0};
    Dim3 dst_end = {0,0,0};
    int scomp = 0;
    int ncomp = 0;

    template <class T>
    AMREX_GPU_HOST_DEVICE
    Array4<T const> getSrc () { return Array4<T const>(static_cast<T const*>(src), src_begin, src_end); }

    template <class T>
    AMREX_GPU_HOST_DEVICE
    Array4<T> getDst () { return Array4<T>(static_cast<T*>(dst), dst_begin, dst_end); }
};

template <typename T, typename U>
CopyMemory
makeCopyMemory (Array4<T> const& src, Array4<U> const& dst, int scomp, int ncomp)
{
    return { src.p, dst.p, src.begin, src.end, dst.begin, dst.end, scomp, ncomp };
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
        m_parms_d = static_cast<T*>( The_Device_Arena()->alloc(sizeof(T)*m_parms.size()) );
    }
    ~CudaGraph()
    {
        The_Device_Arena()->free(m_parms_d);
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
