/* Copyright 2020 Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARSER_WRAPPER_H_
#define WARPX_PARSER_WRAPPER_H_

#include "Parser/WarpXParser.H"
#include "Parser/GpuParser.H"

#include <AMReX_Gpu.H>

/**
 * \brief
 *  The ParserWrapper struct is constructed to safely use the GpuParser class
 *  that can essentially be though of as a raw pointer. The GpuParser does
 *  not have an explicit destructor and the AddPlasma subroutines handle the GpuParser
 *  in a safe way. The ParserWrapper struct is used to avoid memory leak
 *  in the EB parser functions.
 */
template <int N>
struct ParserWrapper
    : public amrex::Gpu::Managed, public GpuParser<N>
{
    using GpuParser<N>::GpuParser;

    ParserWrapper (ParserWrapper<N> const&) = delete;
    void operator= (ParserWrapper<N> const&) = delete;

    ~ParserWrapper() { GpuParser<N>::clear(); }
};

#endif
