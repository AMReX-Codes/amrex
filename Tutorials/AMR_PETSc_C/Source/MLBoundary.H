#ifndef _MLBoundary_H_
#define _MLBoundary_H_

#include <MFTower.H>
#include <BCRec.H>
#include <Layout.H>

#include <petsc.h>

class MLBoundary
{
public:

  MLBoundary(Layout&      layout,
	     const BCRec& bc);

  ~MLBoundary();

  const BCRec& BC() const {return bc;}

  virtual void DefineDirichletValues(FArrayBox&         bcfab, 
                                     const Orientation& face, 
                                     int                dComp, 
                                     int                nComp);

  virtual void DefineNeumannValues(FArrayBox&         bcfab, 
                                   const Orientation& face, 
                                   int                dComp, 
                                   int                nComp);
  
  void SetDirichletValues(MFTower& pressure,
			  int      pComp);

  void SetInflowFlux(PArray<MFTower>& flux,
                     int              fComp);
  
  void SetInflowFlux(MultiFab& fmf,
                     int       fComp,
                     int       lev,
                     int       d);

  void SetDirichletValues(MultiFab& pmf,
                          int       pComp,
                          int       lev); 

  Array<std::map<int, FArrayBox* > > bc_pressure_values; // val[lev][gidx] is a Fab of values
  Array< Array<std::map<int, FArrayBox* > > > bc_flux_values; // flux[d][lev][gidx] is a Fab of fluxes
  Array<Orientation> dirichlet_faces;
  Array<Orientation> neumann_faces;

protected:
  BCRec bc;
  int nLevs;
  Real bcval, influx;
};

#endif
