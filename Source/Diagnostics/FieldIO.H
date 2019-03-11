#ifndef WARPX_FielIO_H_
#define WARPX_FielIO_H_

#include <WarpX.H>

using namespace amrex;

void
PackPlotDataPtrs (Vector<const MultiFab*>& pmf,
		  const std::array<std::unique_ptr<MultiFab>,3>& data);

#endif
