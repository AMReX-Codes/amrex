#ifndef WARPX_FielIO_H_
#define WARPX_FielIO_H_

#include <WarpX.H>

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
             const int lev, const bool plot_guards );

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

#endif
