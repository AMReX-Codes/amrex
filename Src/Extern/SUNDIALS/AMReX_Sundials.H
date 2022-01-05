#ifndef AMREX_SUNDIALS_H_
#define AMREX_SUNDIALS_H_

#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_NVector_MultiFab.H>
static_assert(std::is_same<amrex::Real,realtype>::value, "amrex::Real must be the same as SUNDIALS realtype");
#include <AMReX_SUNMemory.H>

namespace amrex {
namespace sundials {

void Initialize(int nthreads);
void Finalize();
::sundials::Context* The_Sundials_Context(int i = amrex::OpenMP::get_thread_num());

}//sundials
}//amrex

#endif
