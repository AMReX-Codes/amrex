#include "MyKernel.H"

__global__
void plusone (double *data)
{
    *data += 1.0;
}
