#include "MyKernel.H"

#ifdef USE_CUDA
__global__
#endif
void plusone (double *data)
{
    *data += 1.0;
}
