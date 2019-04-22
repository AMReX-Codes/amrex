#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_CudaDevice.H>

namespace amrex {

struct CopyMemory
{
    Array4<Real> src{nullptr,{},{}};
    Array4<Real> dst{nullptr,{},{}};
    Dim3 offset = {0,0,0};
    int scomp = 0;
    int ncomp = 0;
};

struct CudaGraph
{
    cudaGraphExec_t copyGraph;
    CopyMemory* mem = nullptr;
    CopyMemory* mem_d = nullptr;
    int size;

    CudaGraph(int num)
        : size(num)
    {
        // Placement new to allow graph creation before host memory definition. 
        mem = static_cast<CopyMemory*>(std::malloc(sizeof(CopyMemory)*size)); 
        for(int i=0; i<size; ++i)
        {
           new (mem+i) CopyMemory();
        }

        cudaMalloc(&mem_d, sizeof(CopyMemory)*size);
        cudaMemcpy(mem_d, mem, size_t(sizeof(CopyMemory)*size), cudaMemcpyHostToDevice);
    }
    ~CudaGraph()
    {
        for(int i=0; i<size; ++i)
        {
            (mem+i)->~CopyMemory();
        }
        std::free(mem);

        cudaFree(mem_d);         
    }
    void setGraph(cudaGraphExec_t &graph)
    { 
        copyGraph = graph;
    }
    void setParams(int idx, CopyMemory &set)
    {
        std::memcpy((mem+idx), &set, sizeof(CopyMemory)); 
    }
    void executeGraph()
    {
        cudaMemcpy(mem_d, mem, size_t(sizeof(CopyMemory)*size), cudaMemcpyHostToDevice);

        Cuda::Device::setStreamIndex(0);

        cudaGraphLaunch(copyGraph, amrex::Cuda::Device::cudaStream());

        Cuda::Device::resetStreamIndex();
        Cuda::Device::synchronize();
    }
};

}
