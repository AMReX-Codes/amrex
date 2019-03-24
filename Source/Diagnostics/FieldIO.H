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
coarsenCellCenteredFields(
	Vector<MultiFab>& coarse_mf, Vector<Geometry>& coarse_geom,
	const Vector<MultiFab>& source_mf, const Vector<Geometry>& source_geom,
	int coarse_ratio, int finest_level );

void
WriteOpenPMDFields( const std::string& filename,
                    const std::vector<std::string>& varnames,
                    const MultiFab& mf, const Geometry& geom,
                    const int iteration, const double time );

#endif
