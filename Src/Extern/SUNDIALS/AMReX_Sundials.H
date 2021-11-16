#ifndef AMREX_SUNDIALS_H_
#define AMREX_SUNDIALS_H_

#include <AMReX_SUNMemory.H>
#include <AMReX_NVector_MultiFab.H>
static_assert(std::is_same<amrex::Real,realtype>::value, "amrex::Real must be the same as SUNDIALS realtype");
#include <sundials/sundials_context.h>

namespace amrex {
namespace sundials {

void Initialize();
void Finalize();
::sundials::Context* The_Sundials_Context();

}//sundials
}//amrex

#endif
