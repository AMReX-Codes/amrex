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
}

}

#endif 

#endif
