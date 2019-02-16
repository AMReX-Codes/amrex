#ifndef _MYKERNEL_H_
#define _MYKERNEL_H_

#ifdef USE_CUDA
__global__
#endif
void plusone (double *data);

#endif

