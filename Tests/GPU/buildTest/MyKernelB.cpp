#include "MyKernelB.H"

__global__
void minusone (double *data)
{
    *data -= 1.0;
}
