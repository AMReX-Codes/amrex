#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_CudaDevice.H>

namespace amrex {

template <typename T>
struct CopyMemory
{
    Array4<T const> src{nullptr,{},{}};
    Array4<T> dst{nullptr,{},{}};
    Dim3 offset = {0,0,0};
    int scomp = 0;
    int ncomp = 0;
};

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
//        static_assert(std::is_trivially_copyable<T>::value, "CudaGraph's T must be trivially copyable");
//        cudaMemcpy(m_parms_d, m_parms.data(), std::size_t(sizeof(T)*m_parms.size()),
//                   cudaMemcpyHostToDevice);
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
