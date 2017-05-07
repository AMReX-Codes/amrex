#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    int err_c_ret_int_f_ret_real ();

    int pass_c_ret_int_f_ret_int ();

    void err_c_2_f_4 (int* i, int* j);

    void err_c_value_f_pointer (double x);

    void err_c_poiner_f_value (amrex::Real* x);

    void err_c_reference_f_value (const amrex::Real& x);

    void pass_c_bl_fort_fab_arg (BL_FORT_FAB_ARG(fab));

    void pass_c_bl_fort_fab_arg_3d (BL_FORT_FAB_ARG_3D(fab));

    void err_c_float_f_double (BL_FORT_IFAB_ARG_3D(ifab), float*);

    void pass_c_real_f_real (float x, double y, amrex_real z);

    void pass_c_pointer_f_c_ptr (int * p);

    void err_c_pointer_f_c_ptr (amrex::Real * p);

    void pass_c_void_pointer_f_any (void* x, void* y, void* z);

#ifdef __cplusplus
}
#endif
