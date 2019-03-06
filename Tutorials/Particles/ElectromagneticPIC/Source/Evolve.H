#ifndef EVOLVE_H_
#define EVOLVE_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

amrex::Real compute_dt(const amrex::Geometry& geom);

void evolve_electric_field(amrex::MultiFab& Ex, amrex::MultiFab& Ey, amrex::MultiFab& Ez,
                           const amrex::MultiFab& Bx,
                           const amrex::MultiFab& By,
                           const amrex::MultiFab& Bz,
                           const amrex::MultiFab& jx,
                           const amrex::MultiFab& jy,
                           const amrex::MultiFab& jz,
                           const amrex::Geometry& geom, amrex::Real dt);

void evolve_magnetic_field(const amrex::MultiFab& Ex,
                           const amrex::MultiFab& Ey,
                           const amrex::MultiFab& Ez,
                           amrex::MultiFab& Bx, amrex::MultiFab& By, amrex::MultiFab& Bz,
                           const amrex::Geometry& geom, amrex::Real dt);

void fill_boundary_electric_field(amrex::MultiFab& Ex,
                                  amrex::MultiFab& Ey,
                                  amrex::MultiFab& Ez,
                                  const amrex::Geometry& geom);

void fill_boundary_magnetic_field(amrex::MultiFab& Bx,
                                  amrex::MultiFab& By,
                                  amrex::MultiFab& Bz,
                                  const amrex::Geometry& geom);

#endif
