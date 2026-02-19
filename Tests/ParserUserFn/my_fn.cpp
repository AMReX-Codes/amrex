#include "my_fn.H"

double f1_h (double a)
{
    return 3*a;
}

AMREX_GPU_DEVICE double f1_d (double a)
{
    return 3*a;
}

AMREX_GPU_HOST_DEVICE double f1_hd (double a)
{
    return 3*a;
}

double f2_h (double a, double b)
{
    return a-3*b;
}

AMREX_GPU_DEVICE double f2_d (double a, double b)
{
    return a-3*b;
}

AMREX_GPU_HOST_DEVICE double f2_hd (double a, double b)
{
    return a-3*b;
}

double f3_h (double a, double b, double c)
{
    return 3*a+b*c;
}

AMREX_GPU_DEVICE double f3_d (double a, double b, double c)
{
    return 3*a+b*c;
}

AMREX_GPU_HOST_DEVICE double f3_hd (double a, double b, double c)
{
    return 3*a+b*c;
}

double f4_h (double a, double b, double c, double d)
{
    return 3*a+b*c-d;
}

AMREX_GPU_DEVICE double f4_d (double a, double b, double c, double d)
{
    return 3*a+b*c-d;
}

AMREX_GPU_HOST_DEVICE double f4_hd (double a, double b, double c, double d)
{
    return 3*a+b*c-d;
}
