#ifndef WARPX_PARSER_WRAPPER_H_
#define WARPX_PARSER_WRAPPER_H_

#include <WarpXParser.H>
#include <AMReX_Gpu.H>
#include <GpuParser.H>
/**
 * \brief
 *  The ParserWrapper struct is constructed to safely use the GpuParser class
 *  that can essentially be though of as a raw pointer. The GpuParser does
 *  not have an explicit destructor and the AddPlasma subroutines handle the GpuParser
 *  in a safe way. The ParserWrapper struct is used to avoid memory leak
 *  in the EB parser functions.
 */
struct ParserWrapper
    : public amrex::Gpu::Managed
{
    ParserWrapper (WarpXParser const& a_parser) noexcept
       : m_parser(a_parser) {}

    ~ParserWrapper() {
       m_parser.clear();
    }

    AMREX_GPU_HOST_DEVICE
    amrex::Real
    getField (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
       return m_parser(x,y,z);
    }

    GpuParser m_parser;
};

#endif
