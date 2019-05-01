#ifndef AMREX_FUNCTORS_H_
#define AMREX_FUNCTORS_H_

#include <AMReX_Geometry.H>
#include <AMReX_BaseFab.H>
#include <AMReX_Particle.H>
#include <AMReX_ParticleUtil.H>

#ifdef AMREX_USE_CUDA

#include <thrust/tuple.h>


namespace amrex
{

namespace functors
{

    struct tupleAdd
    {
        template <typename T>
        AMREX_GPU_HOST_DEVICE
        T operator()(T v1, T v2)  {
            thrust::get<0>(v1) = thrust::get<0>(v2) + thrust::get<0>(v1);
            thrust::get<1>(v1) = thrust::get<1>(v2) + thrust::get<1>(v1);
            return v1 ;
        }
    };
    
    struct applyParticlePeriodicBCs
    {
        GeometryData domain;

        AMREX_GPU_HOST_DEVICE
        applyParticlePeriodicBCs(const GeometryData& a_domain)
            : domain(a_domain) {}

        template <typename Particle>
        AMREX_GPU_HOST_DEVICE
        Particle operator()(const Particle& in) const {
            Particle p = in;
            
            for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
	      if (not domain.isPeriodic(dim)) continue;
	      if (p.m_rdata.pos[dim] > domain.ProbHi()[dim]) { 
		p.m_rdata.pos[dim] -= (domain.ProbHi()[dim]-domain.ProbLo()[dim]); 
	      } else if (p.m_rdata.pos[dim] < domain.ProbLo()[dim]) {
		p.m_rdata.pos[dim] += (domain.ProbHi()[dim]-domain.ProbLo()[dim]);
	      }
            }
            
            return p;
        }
    };

    struct assignParticleGrid
    {        
        GeometryData domain;
        Box box;
        BaseFab<int>* mask_ptr;
        
        AMREX_GPU_HOST_DEVICE
        assignParticleGrid(const GeometryData& a_domain,
                           const Box&          a_box,
                           BaseFab<int>*       a_mask_ptr) 
            : domain(a_domain), box(a_box), mask_ptr(a_mask_ptr) {}
        
        template <typename Particle>
        AMREX_GPU_HOST_DEVICE
        int operator()(const Particle& p) const {
	    IntVect iv = IntVect(
	AMREX_D_DECL(floor((p.m_rdata.pos[0]-domain.ProbLo()[0])/domain.CellSize()[0]),
		     floor((p.m_rdata.pos[1]-domain.ProbLo()[1])/domain.CellSize()[1]),
		     floor((p.m_rdata.pos[2]-domain.ProbLo()[2])/domain.CellSize()[2]))
				 );

	    iv += domain.Domain().smallEnd();

	    return (*mask_ptr)(iv); 
        }
    };

    struct assignParticleTile
    {        
        GeometryData domain;
        Box box;
        BaseFab<int>* mask_ptr;
        
        AMREX_GPU_HOST_DEVICE
        assignParticleTile(const GeometryData& a_domain,
			   const Box&          a_box,
			   BaseFab<int>*       a_mask_ptr) 
            : domain(a_domain), box(a_box), mask_ptr(a_mask_ptr) {}
        
        template <typename Particle>
        AMREX_GPU_HOST_DEVICE
        int operator()(const Particle& p) const {
	    IntVect iv = IntVect(
	AMREX_D_DECL(floor((p.m_rdata.pos[0]-domain.ProbLo()[0])/domain.CellSize()[0]),
		     floor((p.m_rdata.pos[1]-domain.ProbLo()[1])/domain.CellSize()[1]),
		     floor((p.m_rdata.pos[2]-domain.ProbLo()[2])/domain.CellSize()[2]))
				 );

	    iv += domain.Domain().smallEnd();

	    return (*mask_ptr)(iv, 1);
        }
    };

    struct periodicTransformation
    {
        GeometryData domain;
        IntVect pshift;

        AMREX_GPU_HOST_DEVICE
        periodicTransformation(const GeometryData& a_domain,
                               const IntVect&      a_pshift)
            : domain(a_domain), pshift(a_pshift) {}
        
        template <typename Particle>
        AMREX_GPU_HOST_DEVICE
        Particle operator()(const Particle& in) const {
            Particle p = in;
            
            for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
                if (not domain.isPeriodic(dim)) continue;
                if (pshift[dim] > 0) 
                    p.m_rdata.pos[dim] += (domain.ProbHi()[dim]-domain.ProbLo()[dim]); 
                else if (pshift[dim] < 0)
                    p.m_rdata.pos[dim] -= (domain.ProbHi()[dim]-domain.ProbLo()[dim]);
            }
            
            return p;
        }
    };

    struct assignParticleBin
    {
        GeometryData domain;
        Box box;
        IntVect bin_size;

        AMREX_GPU_HOST_DEVICE
        assignParticleBin(const GeometryData& a_domain,
			  const Box&          a_box,
			  const IntVect&      a_bin_size) 
            : domain(a_domain), box(a_box), bin_size(a_bin_size) {}
      
        template <typename Particle>
        AMREX_GPU_HOST_DEVICE
        int operator()(const Particle& p) const {
	    IntVect iv = IntVect(
	AMREX_D_DECL(floor((p.m_rdata.pos[0]-domain.ProbLo()[0])/domain.CellSize()[0]),
		     floor((p.m_rdata.pos[1]-domain.ProbLo()[1])/domain.CellSize()[1]),
		     floor((p.m_rdata.pos[2]-domain.ProbLo()[2])/domain.CellSize()[2]))
				 );

	    iv += domain.Domain().smallEnd();
	    
	    Box tbx;
	    return getTileIndex(iv, box, true, bin_size, tbx);	    
	}
    };

    struct assignParticleGridAndTile
    {        
        GeometryData domain;
        Box box;
        BaseFab<int>* mask_ptr;
        
        AMREX_GPU_HOST_DEVICE
        assignParticleGridAndTile(const GeometryData& a_domain,
				  const Box&          a_box,
				  BaseFab<int>*       a_mask_ptr) 
            : domain(a_domain), box(a_box), mask_ptr(a_mask_ptr) {}
        
        template <typename Particle>
        AMREX_GPU_HOST_DEVICE
        thrust::tuple<int, int> operator()(const Particle& p) const {
	    IntVect iv = IntVect(
	AMREX_D_DECL(floor((p.m_rdata.pos[0]-domain.ProbLo()[0])/domain.CellSize()[0]),
		     floor((p.m_rdata.pos[1]-domain.ProbLo()[1])/domain.CellSize()[1]),
		     floor((p.m_rdata.pos[2]-domain.ProbLo()[2])/domain.CellSize()[2]))
				 );

	    iv += domain.Domain().smallEnd();

	    return thrust::make_tuple((*mask_ptr)(iv, 0), (*mask_ptr)(iv, 1));
        }
    };

    struct assignParticleCell
    {        
        GeometryData domain;
        
        AMREX_GPU_HOST_DEVICE
        assignParticleCell(const GeometryData& a_domain)
            : domain(a_domain) {}
        
        template <typename Particle>
        AMREX_GPU_HOST_DEVICE
        IntVect operator()(const Particle& p) const {
	    IntVect iv = IntVect(
	AMREX_D_DECL(floor((p.m_rdata.pos[0]-domain.ProbLo()[0])/domain.CellSize()[0]),
		     floor((p.m_rdata.pos[1]-domain.ProbLo()[1])/domain.CellSize()[1]),
		     floor((p.m_rdata.pos[2]-domain.ProbLo()[2])/domain.CellSize()[2]))
				 );

	    iv += domain.Domain().smallEnd();

	    return iv;
        }
    };
    
    struct grid_is
    {
        int grid_id;
        
        AMREX_GPU_HOST_DEVICE
        grid_is(int a_grid_id) : grid_id(a_grid_id) {}
        
        AMREX_GPU_HOST_DEVICE
        bool operator()(int grid) const
        {
            return grid_id == grid;
        }
    };
    
    struct grid_is_not
    {
        int grid_id;
        
        AMREX_GPU_HOST_DEVICE
        grid_is_not(int a_grid_id) : grid_id(a_grid_id) {}
    
        AMREX_GPU_HOST_DEVICE
        bool operator()(int grid) const
        {
            return grid_id != grid;
        }
    };

    struct getPosition
    {        
        AMREX_GPU_HOST_DEVICE
        getPosition() {}
        
        template <typename Particle>
        AMREX_GPU_HOST_DEVICE
        thrust::tuple<AMREX_D_DECL(Real, Real, Real)> operator()(const Particle& p) const
        {
            return thrust::make_tuple(AMREX_D_DECL(p.m_rdata.pos[0],
                                                   p.m_rdata.pos[1], 
                                                   p.m_rdata.pos[2]));
        }
    };

    template <typename Particle>
    struct setPosition
    {
        AMREX_GPU_HOST_DEVICE
        setPosition() {}

        template <typename Tuple>
        AMREX_GPU_HOST_DEVICE
        Particle operator()(const Tuple& pos) const
        {
            Particle p;
            p.m_rdata.pos[0] = thrust::get<0>(pos);
#if AMREX_SPACEDIM > 1
            p.m_rdata.pos[1] = thrust::get<1>(pos);
#endif
#if AMREX_SPACEDIM > 2
            p.m_rdata.pos[2] = thrust::get<2>(pos);
#endif
            return p;
        }
    };
}

}

#endif 

#endif
