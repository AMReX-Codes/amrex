#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

void main_main ();

void doRandomSleep (amrex::MultiFab& mf);

void init_mf (amrex::MultiFab& mf_new, amrex::Geometry const& geom);
#endif
