#ifndef IO_H_
#define IO_H_

#include <AMReX_REAL.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

void
WritePlotFile (const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
               const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
               const amrex::MultiFab& jx, const amrex::MultiFab& jy, const amrex::MultiFab& jz,
               const amrex::Geometry& geom, amrex::Real time, int step);

#endif
