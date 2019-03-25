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

#endif
