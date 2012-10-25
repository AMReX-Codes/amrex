#include <MLBoundary.H>
#include <ParmParse.H>

static int  pressure_comp_bc = 0; // Component for pressure values in internal boundary value fab
static int  flux_comp_bc = 0; // Component for flux values in internal boundary value fab

static int nGrow = 1;
static int nComp = 1;

static Real bcval_DEF = 0;
static Real influx_DEF = 1.e-7;

MLBoundary::MLBoundary(Layout&      layout,
		       const BCRec& _bc)
  : bc(_bc)
{
  nLevs = layout.NumLevels();

  ParmParse pp("mlb");
  influx = influx_DEF; pp.query("influx",influx);
  bcval = bcval_DEF; pp.query("bcval",bcval);

  for (int d=0; d<BL_SPACEDIM; ++d) {
    if (bc.lo(d) == EXT_DIR) {
      dirichlet_faces.push_back(Orientation(d,Orientation::low));
    }
    else if (bc.lo(d) == FOEXTRAP) {
      neumann_faces.push_back(Orientation(d,Orientation::low));
    }

    if (bc.hi(d) == EXT_DIR) {
      dirichlet_faces.push_back(Orientation(d,Orientation::high));
    }
    else if (bc.hi(d) == FOEXTRAP) {
      neumann_faces.push_back(Orientation(d,Orientation::high));
    }
  }

  bc_pressure_values.resize(nLevs);

  std::vector< std::pair<int,Box> > isects;
  for (int i=0; i<dirichlet_faces.size(); ++i) {
    const Orientation& dface = dirichlet_faces[i];
    for (int lev=0; lev<nLevs; ++lev) {
      Box gbox = BoxLib::adjCell(layout.GeomArray()[lev].Domain(), dface, nGrow);
      BoxArray gba = BoxArray(layout.GridArray()[lev]).grow(nGrow);
      isects = gba.intersections(gbox);
      if (isects.size()>0) {
        const DistributionMapping& dm = layout.DistributionMap(lev);
        for (int j=0; j<isects.size(); ++j) {
          int idx = isects[j].first;
          const Box& bbox = isects[j].second;
          if (dm[idx] == ParallelDescriptor::MyProc()) {
            FArrayBox* fptr = new FArrayBox(bbox,nComp);
            bc_pressure_values[lev][idx].push_back(fptr);
            int bsize = bc_pressure_values[lev][idx].size();
            DefineDirichletValues(*bc_pressure_values[lev][idx][bsize-1],dface,0,nComp);
          }
        }
      }
    }
  }

  bc_flux_values.resize(BL_SPACEDIM,Array<std::map<int, Array<FArrayBox*> > >(nLevs));
  for (int i=0; i<neumann_faces.size(); ++i) {
    const Orientation& nface = neumann_faces[i];
    int dir = nface.coordDir();
    for (int lev=0; lev<nLevs; ++lev) {
      Box ebox = BoxLib::bdryNode(layout.GeomArray()[lev].Domain(),nface,nGrow);
      BoxArray eba = BoxArray(layout.GridArray()[lev]).surroundingNodes(dir);
      isects = eba.intersections(ebox);
      if (isects.size()>0) {
        const DistributionMapping& dm = layout.DistributionMap(lev);
        for (int j=0; j<isects.size(); ++j) {
          int idx = isects[j].first;
          const Box& bbox = isects[j].second;
          if (dm[idx] == ParallelDescriptor::MyProc()) {
            FArrayBox* fptr = new FArrayBox(bbox,nComp);
            bc_flux_values[dir][lev][idx].push_back(fptr);
            int bsize = bc_flux_values[dir][lev][idx].size();
            DefineNeumannValues(*(bc_flux_values[dir][lev][idx][bsize-1]),nface,0,nComp);
          }
        }
      }
    }
  }
}

static Orientation XLO=Orientation(0,Orientation::low);
static Orientation XHI=Orientation(0,Orientation::high);
static Orientation YLO=Orientation(1,Orientation::low);
static Orientation YHI=Orientation(1,Orientation::high);

void 
MLBoundary::DefineDirichletValues(FArrayBox&         bcfab, 
                                  const Orientation& face, 
                                  int                dComp, 
                                  int                nComp)
{
  bcfab.setVal(bcval,bcfab.box(),dComp,nComp);    
}

void 
MLBoundary::DefineNeumannValues(FArrayBox&         bcfab, 
                                const Orientation& face, 
                                int                dComp, 
                                int                nComp)
{
  int sgn = (face.isLow() ? +1 : -1);
  bcfab.setVal(sgn * influx,bcfab.box(),dComp,nComp);
}

void 
MLBoundary::SetInflowFlux(PArray<MFTower>& flux,
			  int              fComp)
{
  for (int d=0; d<BL_SPACEDIM; ++d) {
    for (int lev=0; lev<nLevs; ++lev) {
      MultiFab& fmf = flux[d][lev];
      for (std::map<int, Array<FArrayBox*> >::const_iterator it=bc_flux_values[d][lev].begin(), 
             End=bc_flux_values[d][lev].end(); it!=End; ++it) {
        FArrayBox& ffab = fmf[it->first];
        const Array<FArrayBox*>& fluxes = it->second;
        for (int j=0; j<fluxes.size(); ++j) {
          ffab.copy(*fluxes[j],flux_comp_bc,fComp+flux[d].BaseComp());
        }
      }
    }
  }
}

void
MLBoundary::SetDirichletValues(MFTower& pressure,
			       int      pComp)
{
  for (int lev=0; lev<nLevs; ++lev) {
    MultiFab& pmf = pressure[lev];
    for (std::map<int, Array<FArrayBox*> >::const_iterator it=bc_pressure_values[lev].begin(), 
           End=bc_pressure_values[lev].end(); it!=End; ++it) {
      FArrayBox& pfab = pmf[it->first];
      const Array<FArrayBox*>& vals = it->second;
      for (int j=0; j<vals.size(); ++j) {
        pfab.copy(*vals[j],pressure_comp_bc,pComp+pressure.BaseComp());
      }
    }
  }
}
