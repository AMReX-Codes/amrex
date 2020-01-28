/* Copyright 2019 Axel Huebl, David Grote, Igor Andriyash
 * Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_FielIO_H_
#define WARPX_FielIO_H_

#include <WarpX.H>
#ifdef WARPX_USE_OPENPMD
#   include <openPMD/openPMD.hpp>
#endif

void
PackPlotDataPtrs (amrex::Vector<const amrex::MultiFab*>& pmf,
          const std::array<std::unique_ptr<amrex::MultiFab>,3>& data);

void
AverageAndPackVectorField( amrex::MultiFab& mf_avg,
                         const std::array< std::unique_ptr<amrex::MultiFab>, 3 >& vector_field,
                         const amrex::DistributionMapping& dm,
                         const int dcomp, const int ngrow );

void
AverageAndPackScalarField( amrex::MultiFab& mf_avg,
                         const amrex::MultiFab & scalar_field,
                         const int dcomp, const int ngrow );

void
WriteRawField( const amrex::MultiFab& F,
               const amrex::DistributionMapping& dm,
               const std::string& filename,
               const std::string& level_prefix,
               const std::string& field_name,
               const int lev, const bool plot_guards );

void
WriteZeroRawField( const amrex::MultiFab& F,
             const amrex::DistributionMapping& dm,
             const std::string& filename,
             const std::string& level_prefix,
             const std::string& field_name,
             const int lev, const int ng );

void
WriteCoarseScalar( const std::string field_name,
    const std::unique_ptr<amrex::MultiFab>& F_cp,
    const std::unique_ptr<amrex::MultiFab>& F_fp,
    const amrex::DistributionMapping& dm,
    const std::string& filename,
    const std::string& level_prefix,
    const int lev, const bool plot_guards,
    const int icomp=0 );

void
WriteCoarseVector( const std::string field_name,
    const std::unique_ptr<amrex::MultiFab>& Fx_cp,
    const std::unique_ptr<amrex::MultiFab>& Fy_cp,
    const std::unique_ptr<amrex::MultiFab>& Fz_cp,
    const std::unique_ptr<amrex::MultiFab>& Fx_fp,
    const std::unique_ptr<amrex::MultiFab>& Fy_fp,
    const std::unique_ptr<amrex::MultiFab>& Fz_fp,
    const amrex::DistributionMapping& dm,
    const std::string& filename,
    const std::string& level_prefix,
    const int lev, const bool plot_guards );

std::unique_ptr<amrex::MultiFab>
getInterpolatedScalar(
    const amrex::MultiFab& F_cp, const amrex::MultiFab& F_fp,
    const amrex::DistributionMapping& dm, const int r_ratio,
    const amrex::Real* dx, const int ngrow );

std::array<std::unique_ptr<amrex::MultiFab>, 3>
getInterpolatedVector(
    const std::unique_ptr<amrex::MultiFab>& Fx_cp,
    const std::unique_ptr<amrex::MultiFab>& Fy_cp,
    const std::unique_ptr<amrex::MultiFab>& Fz_cp,
    const std::unique_ptr<amrex::MultiFab>& Fx_fp,
    const std::unique_ptr<amrex::MultiFab>& Fy_fp,
    const std::unique_ptr<amrex::MultiFab>& Fz_fp,
    const amrex::DistributionMapping& dm,
    const int r_ratio, const amrex::Real* dx,
    const int ngrow );

void
coarsenCellCenteredFields(
    amrex::Vector<amrex::MultiFab>& coarse_mf,
    amrex::Vector<amrex::Geometry>& coarse_geom,
    const amrex::Vector<amrex::MultiFab>& source_mf,
    const amrex::Vector<amrex::Geometry>& source_geom,
    int coarse_ratio, int finest_level );

#ifdef WARPX_USE_OPENPMD
void
setOpenPMDUnit( openPMD::Mesh mesh, const std::string field_name );

std::vector<std::uint64_t>
getReversedVec( const amrex::IntVect& v );

std::vector<double>
getReversedVec( const amrex::Real* v );

void
WriteOpenPMDFields( const std::string& filename,
                    const std::vector<std::string>& varnames,
                    const amrex::MultiFab& mf, const amrex::Geometry& geom,
                    const int iteration, const double time );
#endif // WARPX_USE_OPENPMD

#endif // WARPX_FielIO_H_
