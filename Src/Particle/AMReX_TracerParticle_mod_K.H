#ifndef AMREX_TRACERPARTICLE_MOD_K_H
#define AMREX_TRACERPARTICLE_MOD_K_H
#include <AMReX_TracerParticles.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>
#include <cmath>
#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>

//template <typename P> 
//AMREX_GPU_HOST_DEVICE AMREX_INLINE
//void amrex_interpolate_CIC_2(const int p, 
//			     const GeometryData& geomdata, 
//			     array4<Real> const& uccar, 
//			     int idx, 
//			     int val, 
//			     int comp){

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void amrex_interpolate_CIC_2 (const P& p, 
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& idx,
			      const amrex::Array4<const double> &  uccarr,
			      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> &  val,
			      const int comp)  
{


  //BL_Profile("Particle<Nreal, Nint>::Interp()");
  //BL_ASSERT(idx != 0);
  BL_ASSERT(val != 0);
  
  const int   M   = AMREX_D_TERM(2,+2,+4);
  //const amrex::FArrayBox fab = *p_fab;

  amrex::Real    fracs[M];
  amrex::IntVect cells[M];
  
  const amrex::Real len[AMREX_SPACEDIM] = { AMREX_D_DECL(p.m_rdata.pos[0] - plo[0]*idx[0] + 
							 amrex::Real(0.5), 
							 p.m_rdata.pos[1] - plo[1]*idx[1] + amrex::Real(0.5), 
							 p.m_rdata.pos[2] - plo[2]*idx[2] + amrex::Real(0.5))};

  amrex::IntVect cell(AMREX_D_DECL(std::floor(len[0]), std::floor(len[1]), std::floor(len[2])));

  const amrex::Real frac[AMREX_SPACEDIM] = { AMREX_D_DECL(len[0]-cell[0], len[1]-cell[1], len[2]-cell[2])};
  

 
#if (AMREX_SPACEDIM == 1)
  // high
  fracs[0] = frac[0];
  cells[0] = cell;
  // low
  fracs[1] = 1-frac[0];
  cell[0]  = cell[0] - 1;
  cells[1] = cell;

#elif (AMREX_SPACEDIM == 2)
  // HH
  fracs[0] = frac[0] * frac[1];
  cells[0] = cell;
  // LH
  fracs[1] = (1-frac[0])*frac[1];
  cell[0]  = cell[0] - 1;
  cells[1] = cell;
  // LL
  fracs[2] = (1-frac[0]) * (1-frac[1]);
  cell[1]  = cell[1] - 1;
  cells[2] = cell;
  // HL
  fracs[3] = frac[0] * (1-frac[1]);
  cell[0]  = cell[0] + 1;
  cells[3] = cell;


#elif (AMREX_SPACEDIM == 3)

  //! HHH                                                                                            
  fracs[0] = frac[0] * frac[1] * frac[2];
  cells[0] = cell;
  // LHH                                                                                             
  fracs[1] = (1-frac[0]) * frac[1] * frac[2];
  cell[0]  = cell[0] - 1;
  cells[1] = cell;
  // LLH                                                                                             
  fracs[2] = (1-frac[0]) * (1-frac[1]) * frac[2];
  cell[1]  = cell[1] - 1;
  cells[2] = cell;
  // HLH                                                                                             
  fracs[3] = frac[0] * (1-frac[1]) * frac[2];
  cell[0]  = cell[0] + 1;
  cells[3] = cell;
  // HHL                                                                                             
  fracs[4] = frac[0] * frac[1] * (1-frac[2]);
  cell[2]  = cell[2] - 1;
  cells[4] = cell;
  // LHL                                                                                             
  fracs[5] = (1-frac[0]) * frac[1] * (1-frac[2]);
  cell[0]  = cell[0] - 1;
  cells[5] = cell;
  // LLL                                                                                             
  fracs[6] = (1-frac[0]) * (1-frac[1]) * (1-frac[2]);
  cell[1]  = cell[1] - 1;
  cells[6] = cell;
  // HLL                                                                                             
  fracs[7] = frac[0] * (1-frac[1]) * (1-frac[2]);
  cell[0]  = cell[0] + 1;
  cells[7] = cell;
                                                  
#else
  amrex::abort("Invalid AMREX_SPACEDIM")
#endif

    for (int i = 0; i < AMREX_SPACEDIM; i++) 
    {
      //BL_ASSERT( (idx[i] >= 0) && (idx[i] < uccarr.nComp()) );
      val[i] = 0; 
	for (int j = 0; j<M; j++)
	{ 
          val[j] += uccarr(i,j, comp)*fracs[j];
	  //amrex::print(val[j]);
        } 
   
    }
}
#endif
