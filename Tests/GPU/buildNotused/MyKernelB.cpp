#include "MyKernelB.H"

#ifdef USE_CUDA
__global__
#endif
void minusone (double *data)
{
    *data -= 1.0;
}
