#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

void main_main ();

void doDeviceSleep (amrex::MultiFab& mf, int& n);

#endif
