#include "NodalFlags.H"

using namespace amrex;

IntVect YeeGrid::Bx_nodal_flag(1,0,0);
IntVect YeeGrid::By_nodal_flag(0,1,0);
IntVect YeeGrid::Bz_nodal_flag(0,0,1);

IntVect YeeGrid::Ex_nodal_flag(0,1,1);
IntVect YeeGrid::Ey_nodal_flag(1,0,1);
IntVect YeeGrid::Ez_nodal_flag(1,1,0);

IntVect YeeGrid::jx_nodal_flag(0,1,1);
IntVect YeeGrid::jy_nodal_flag(1,0,1);
IntVect YeeGrid::jz_nodal_flag(1,1,0);
