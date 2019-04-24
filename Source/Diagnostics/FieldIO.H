#ifndef WARPX_FielIO_H_
#define WARPX_FielIO_H_

#include <WarpX.H>
#ifdef WARPX_USE_OPENPMD
#include <openPMD/openPMD.hpp>
#endif

using namespace amrex;

void
PackPlotDataPtrs (Vector<const MultiFab*>& pmf,
          const std::array<std::unique_ptr<MultiFab>,3>& data);

void
AverageAndPackVectorField( MultiFab& mf_avg,
                         const std::array< std::unique_ptr<MultiFab>, 3 >& vector_field,
                         const int dcomp, const int ngrow );

void
AverageAndPackScalarField( MultiFab& mf_avg,
                         const MultiFab & scalar_field,
                         const int dcomp, const int ngrow );

void
WriteRawField( const MultiFab& F, const DistributionMapping& dm,
             const std::string& filename,
             const std::string& level_prefix,
             const std::string& field_name,
             const int lev, const bool plot_guards );

void
WriteZeroRawField( const MultiFab& F, const DistributionMapping& dm,
             const std::string& filename,
             const std::string& level_prefix,
             const std::string& field_name,
             const int lev, const int ng );

void
WriteCoarseScalar( const std::string field_name,
    const std::unique_ptr<MultiFab>& F_cp,
    const std::unique_ptr<MultiFab>& F_fp,
    const DistributionMapping& dm,
    const std::string& filename,
    const std::string& level_prefix,
    const int lev, const bool plot_guards,
    const int icomp=0 );

void
WriteCoarseVector( const std::string field_name,
    const std::unique_ptr<MultiFab>& Fx_cp,
    const std::unique_ptr<MultiFab>& Fy_cp,
    const std::unique_ptr<MultiFab>& Fz_cp,
    const std::unique_ptr<MultiFab>& Fx_fp,
    const std::unique_ptr<MultiFab>& Fy_fp,
    const std::unique_ptr<MultiFab>& Fz_fp,
    const DistributionMapping& dm,
    const std::string& filename,
    const std::string& level_prefix,
    const int lev, const bool plot_guards );

std::unique_ptr<MultiFab>
getInterpolatedScalar(
    const MultiFab& F_cp, const MultiFab& F_fp,
    const DistributionMapping& dm, const int r_ratio,
    const Real* dx, const int ngrow );

std::array<std::unique_ptr<MultiFab>, 3>
getInterpolatedVector(
    const std::unique_ptr<MultiFab>& Fx_cp,
    const std::unique_ptr<MultiFab>& Fy_cp,
    const std::unique_ptr<MultiFab>& Fz_cp,
    const std::unique_ptr<MultiFab>& Fx_fp,
    const std::unique_ptr<MultiFab>& Fy_fp,
    const std::unique_ptr<MultiFab>& Fz_fp,
    const DistributionMapping& dm,
    const int r_ratio, const Real* dx,
    const int ngrow );

void
coarsenCellCenteredFields(
    Vector<MultiFab>& coarse_mf, Vector<Geometry>& coarse_geom,
    const Vector<MultiFab>& source_mf, const Vector<Geometry>& source_geom,
    int coarse_ratio, int finest_level );

#ifdef WARPX_USE_OPENPMD
void
setOpenPMDUnit( openPMD::Mesh mesh, const std::string field_name );

std::vector<std::uint64_t>
getReversedVec( const IntVect& v );

std::vector<double>
getReversedVec( const Real* v );

void
WriteOpenPMDFields( const std::string& filename,
                    const std::vector<std::string>& varnames,
                    const MultiFab& mf, const Geometry& geom,
                    const int iteration, const double time );
#endif // WARPX_USE_OPENPMD

#endif // WARPX_FielIO_H_
