/* Copyright 2016-2020 Maxence Thevenet, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef FILTERFUNCTORS_H
#define FILTERFUNCTORS_H

#include "WarpX.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Gpu.H>
#include <AMReX_Random.H>

using SuperParticleType = typename WarpXParticleContainer::SuperParticleType;

/**
 * \brief Functor that returns 0 or 1 depending on a random draw per particle
 */
struct RandomFilter
{
    /** constructor
     * \param a_is_active whether the test is active
     * \param a_fraction fraction of particles to select
     */
    RandomFilter(bool a_is_active, amrex::Real a_fraction)
        : m_is_active(a_is_active), m_fraction(a_fraction) {};

    /**
     * \brief draw random number, return 1 if number < m_fraction, 1 otherwise
     * \param p one particle
     * \return whether or not the particle is selected
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator () (const SuperParticleType& p) const noexcept
    {
        if ( !m_is_active ) return 1;
        else if ( amrex::Random() < m_fraction ) return 1;
        else return 0;
    }
private:
    const bool m_is_active = false; //! select all particles if false
    const amrex::Real m_fraction = 1.0; //! range: [0.0:1.0] where 0 is no & 1 is all particles
};

/**
 * \brief Functor that returns 1 if stride divide particle_id, 0 otherwise
 */
struct UniformFilter
{
    /** constructor
     * \param a_is_active whether the test is active
     * \param a_stride one particle every a_stride is written to file
     */
    UniformFilter(bool a_is_active, int a_stride)
        : m_is_active(a_is_active), m_stride(a_stride) {};

    /**
     * \brief return 1 if stride divide particle_id, 0 otherwise
     * \param p one particle
     * \return whether or not the particle is selected
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator () (const SuperParticleType& p) const noexcept
    {
        if ( !m_is_active ) return 1;
        else if ( p.id()%m_stride == 0 ) return 1;
        else return 0;
    }
private:
    const bool m_is_active = false; //! select all particles if false
    const int m_stride = 0; //! selection of every n-th particle
};

/**
 * \brief Functor that returns 0 or 1 depending on a parser selection
 */
struct ParserFilter
{
    /** constructor
     * \param a_is_active whether the test is active
     */
    ParserFilter(bool a_is_active, int a_species_index)
    : m_is_active(a_is_active), m_species_index(a_species_index)
    {
        m_t = WarpX::GetInstance().gett_new(0);
        auto & mypc   = WarpX::GetInstance().GetPartContainer();
        auto & myspc  = mypc.GetParticleContainer(m_species_index);
        m_function_partparser = myspc.m_particle_filter_parser.get();
    };

    /**
     * \brief return 1 if the particle is selected by the parser
     * \param p one particle
     * \return whether or not the particle is selected
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator () (const SuperParticleType& p) const noexcept
    {
        if ( !m_is_active ) return 1;
        amrex::Real const x  = p.pos(0);
        amrex::Real const y  = p.pos(1);
        amrex::Real const z  = p.pos(2);
        amrex::Real const ux = p.rdata(PIdx::ux)/PhysConst::c;
        amrex::Real const uy = p.rdata(PIdx::uy)/PhysConst::c;
        amrex::Real const uz = p.rdata(PIdx::uz)/PhysConst::c;
        if ( (*m_function_partparser)(m_t,x,y,z,ux,uy,uz) > 0.5 ) return 1;
        else return 0;
    }
public:
    amrex::Real m_t;
    ParserWrapper<7> * m_function_partparser;
private:
    const bool m_is_active = false; //! select all particles if false
    const int  m_species_index = 0; //! species index
};

#endif // FILTERFUNCTORS_H
