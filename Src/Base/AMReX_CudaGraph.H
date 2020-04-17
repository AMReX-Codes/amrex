#ifndef AMREX_CUDA_GRAPH_H_
#define AMREX_CUDA_GRAPH_H_

#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10))

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_GpuDevice.H>

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
    Array4<T const> getSrc () { return Array4<T const>(static_cast<T const*>(src), src_begin, src_end, scomp+ncomp); }

    template <class T>
    AMREX_GPU_HOST_DEVICE
    Array4<T> getDst () { return Array4<T>(static_cast<T*>(dst), dst_begin, dst_end, scomp+ncomp); }
};

template <typename T, typename U>
CopyMemory
makeCopyMemory (Array4<T> const& src, Array4<U> const& dst, int scomp, int ncomp)
{
#if __cplusplus < 201402L
    CopyMemory mem;
    mem.src = (void*)(src.p);
    mem.dst = (void*)(dst.p);
    mem.src_begin = src.begin;
    mem.src_end = src.end;
    mem.dst_begin = dst.begin;
    mem.dst_end = dst.end;
    mem.scomp = scomp;
    mem.ncomp = ncomp;
    return mem;

#else

    return CopyMemory{ (void*)(src.p), (void*)(dst.p), src.begin, src.end, dst.begin, dst.end, scomp, ncomp };
#endif
}

// ======================================================================================

template <typename T>
struct CudaGraph
{
    cudaGraphExec_t m_graph;
    Vector<T> m_parms;
    T* m_parms_d = nullptr;
    bool graph_is_ready = false;

    CudaGraph()
        : m_parms(0)
    {
        static_assert(AMREX_IS_TRIVIALLY_COPYABLE(T), "CudaGraph's T must be trivially copyable");
    }
    CudaGraph(int num)
        : m_parms(num)
    {
        static_assert(AMREX_IS_TRIVIALLY_COPYABLE(T), "CudaGraph's T must be trivially copyable");
        m_parms_d = static_cast<T*>( The_Device_Arena()->alloc(sizeof(T)*m_parms.size()) );
    }
    ~CudaGraph() {
        The_Device_Arena()->free(m_parms_d);

        if (graph_is_ready)
        {
            AMREX_CUDA_SAFE_CALL(cudaGraphExecDestroy(m_graph));
        }
    }

    void resize(Long num) {
        m_parms.resize(num);
        if (m_parms_d != nullptr)
        {
            The_Device_Arena()->free(m_parms_d);
        }
        m_parms_d = static_cast<T*>( The_Device_Arena()->alloc(sizeof(T)*m_parms.size()) );
    }
    void setGraph(cudaGraphExec_t const& graph) { 
        m_graph = graph;
        graph_is_ready = true;
    }
    void setParams(int idx, T const& a_parm) { m_parms[idx] = a_parm; }

    T* getHostPtr (int idx) { return (m_parms.data() + idx); }
    T* getDevicePtr (int idx) const { return (m_parms_d + idx); }
    bool ready() const { return graph_is_ready; }

    void executeGraph (bool synch = true) const {
        Gpu::Device::executeGraph(m_graph, synch);
    }
};

}

#endif
#endif
