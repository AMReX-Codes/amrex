/*--------------------------------------------------------------------
  standard includes
  --------------------------------------------------------------------*/
#include <AMReX_BCRec.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <array>

/*--------------------------------------------------------------------
  public free function declarations
  --------------------------------------------------------------------*/
amrex::Vector<amrex::Geometry> DefineGeometry ( int nx, int ny, int nz, double dx );

amrex::Vector<amrex::BoxArray> DefineBoxArrays ( int nx, int ny, int nz );

amrex::Vector<amrex::DistributionMapping> DefineDMs (
    const amrex::Vector<amrex::BoxArray>& bas );

void DefineFABs (  //
    amrex::Vector<std::array<amrex::MultiFab, 3>>& vels,
    amrex::Vector<amrex::MultiFab>& temps,
    const amrex::Vector<amrex::BoxArray>& bas,
    const amrex::Vector<amrex::DistributionMapping>& dms );

void CalculateVelocities (  //
    amrex::Vector<std::array<amrex::MultiFab, 3>>& vels );

void CalculateTemperatures (  //
    amrex::Vector<amrex::MultiFab>& temps,
    const amrex::Vector<amrex::Geometry>& geoms,
    const amrex::BCRec& bcrec );

void CalculateBuoyancy (  //
    amrex::Vector<std::array<amrex::MultiFab, 3>>& vels,
    const amrex::Vector<amrex::MultiFab>& temps,
    amrex::Real dt,
    const std::array<amrex::Real, 3>& gs,
    amrex::Real beta,
    amrex::Real T_ref );

void CheckResults (  //
    const amrex::Vector<std::array<amrex::MultiFab, 3>>& vels,
    const amrex::Vector<amrex::Geometry>& geoms );

/*--------------------------------------------------------------------
  End of file
  --------------------------------------------------------------------*/
