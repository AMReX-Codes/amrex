#include <cuda_runtime.h>
#include <iostream>

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

__global__ void func (double *ptr)
{
    *ptr = 3;
}

template<class L>
__global__ void lambda (L f0) { f0(); }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

int main (int argc, char* argv[])
{
    std::cout << "**********************************\n";
    // ===================================
    // Simple cuda action to make sure all tests have cuda.
    // Allows nvprof to return data.
    int devices = 0;
    cudaGetDeviceCount(&devices);

    std::cout << "GPU devices: " << devices << "\n";
    std::cout << "**********************************\n"; 
    // ===================================

    {
        // Malloc value for double testing
        double *p, *q;
        cudaMallocManaged(&p, sizeof(double));
        cudaMallocManaged(&q, sizeof(double));
        *p = 0.0;
        *q = 0.0;

        cudaStream_t    stream;
        cudaStreamCreate(&stream);

        cudaGraph_t     graphP, graphQ;
        cudaGraphExec_t graphExec;

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Create P graph, instantiate and execute. 

        {
#if (__CUDACC_VER_MAJOR__ == 10) && (__CUDACC_VER_MINOR__ == 0)
            cudaStreamBeginCapture(stream);
#else
            cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
#endif
            func<<<1,1,0,stream>>>(p);
            cudaStreamEndCapture(stream, &(graphP));

            cudaGraphInstantiate(&graphExec, graphP, NULL, NULL, 0);
            cudaGraphLaunch(graphExec, stream);  

            cudaDeviceSynchronize();

            std::cout << "1: p, q = " << *p << ", " << *q << std::endl;

// ------------------------------------------------------------------------------- 
//          Make Q graph

#if (__CUDACC_VER_MAJOR__ == 10) && (__CUDACC_VER_MINOR__ == 0)
            cudaStreamBeginCapture(stream);
#else
            cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
#endif

            func<<<1,1,0,stream>>>(q);
            cudaStreamEndCapture(stream, &(graphQ));

// ------------------------------------------------------------------------------- 
//          Get nodes from graphs

            size_t numP, numQ;
            cudaGraphGetNodes(graphP, NULL, &numP);
            cudaGraphGetNodes(graphQ, NULL, &numQ); 

            std::cout << "numP, numQ = " << numP << ", " << numQ << std::endl;

            cudaGraphNode_t nodesP[numP], nodesQ[numQ];
            cudaGraphGetNodes(graphP, nodesP, &numP);
            cudaGraphGetNodes(graphQ, nodesQ, &numQ);

// ------------------------------------------------------------------------------- 
//          Swap params in CudaGraphExec and re-execute without instantiation.

            cudaKernelNodeParams knp;
            cudaGraphKernelNodeGetParams(nodesQ[0], &knp);
#if (__CUDACC_VER_MAJOR__ == 10) && (__CUDACC_VER_MINOR__ == 0)
            cudaGraphExecKernelNodeSetParams(graphExec, nodesP[0], &knp);
#endif
            cudaGraphLaunch(graphExec, stream);  
            cudaDeviceSynchronize();

            // "2: p, q = 3, 3"?
            std::cout << "2: p, q = " << *p << ", " << *q << std::endl;
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        std::cout << "Double Test Completed." << std::endl << std::endl; 
    }

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    {
        std::cout << "**********************************\n"; 
        // ===================================

        // Malloc values for lambda testing.
        double *p, *q;
        cudaMallocManaged(&p, sizeof(double));
        cudaMallocManaged(&q, sizeof(double));
        *p = 0.0;
        *q = 0.0;

        cudaStream_t    stream;
        cudaStreamCreate(&stream);

        cudaGraph_t     graphP, graphQ;
        cudaGraphExec_t graphExec;

        auto plus  = [=] __device__ (double* p) { *p = *p + 1; };
        auto minus = [=] __device__ (double* p) { *p = *p - 1; };

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Create P graph, instantiate and execute. 

        {
#if (__CUDACC_VER_MAJOR__ == 10) && (__CUDACC_VER_MINOR__ == 0)
            cudaStreamBeginCapture(stream);
#else
            cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
#endif

            lambda<<<1,1,0,stream>>>( [=] __device__ () { plus(p); } );
            cudaStreamEndCapture(stream, &(graphP));

            cudaGraphInstantiate(&graphExec, graphP, NULL, NULL, 0);
            cudaGraphLaunch(graphExec, stream);  

            cudaDeviceSynchronize();

            std::cout << "1: p = " << *p << std::endl;

// ------------------------------------------------------------------------------- 
//          Make Q graph

#if (__CUDACC_VER_MAJOR__ == 10) && (__CUDACC_VER_MINOR__ == 0)
            cudaStreamBeginCapture(stream);
#else
            cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
#endif

            lambda<<<1,1,0,stream>>>( [=] __device__ () { minus(p); } );
            cudaStreamEndCapture(stream, &(graphQ));

// ------------------------------------------------------------------------------- 
//          Get nodes from graphs

            size_t numP, numQ;
            cudaGraphGetNodes(graphP, NULL, &numP);
            cudaGraphGetNodes(graphQ, NULL, &numQ); 

            std::cout << "numP, numQ = " << numP << ", " << numQ << std::endl;

            cudaGraphNode_t nodesP[numP], nodesQ[numQ];
            cudaGraphGetNodes(graphP, nodesP, &numP);
            cudaGraphGetNodes(graphQ, nodesQ, &numQ);

// ------------------------------------------------------------------------------- 
//          Swap params in CudaGraphExec and re-execute without instantiation.

            cudaKernelNodeParams knp, knpB;
            cudaGraphKernelNodeGetParams(nodesQ[0], &knp);

#if (__CUDACC_VER_MAJOR__ == 10) && (__CUDACC_VER_MINOR__ == 0)
            cudaGraphExecKernelNodeSetParams(graphExec, nodesP[0], &knp);
#endif
            cudaGraphLaunch(graphExec, stream);  
            cudaDeviceSynchronize();

            std::cout << "2: p = " << *p << std::endl;
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        std::cout << "Lambda Test Completed." << std::endl << std::endl;
    }

}
