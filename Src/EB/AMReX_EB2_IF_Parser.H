#ifndef AMREX_EB2_IF_PARSER_H_
#define AMREX_EB2_IF_PARSER_H_
#include <AMReX_Config.H>

#include <AMReX_EB2_IF_Base.H>
#include <AMReX_Parser.H>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex::EB2 {

class ParserIF
    : public amrex::GPUable
{
public:
    ParserIF (const ParserExecutor<3>& a_parser)
        : m_parser(a_parser)
        {}

    AMREX_GPU_HOST_DEVICE inline
    amrex::Real operator() (AMREX_D_DECL(amrex::Real x, amrex::Real y,
                                         amrex::Real z)) const noexcept {
#if (AMREX_SPACEDIM == 2)
        return static_cast<Real>(m_parser({x,y,Real(0.0)}));
#else
        return static_cast<Real>(m_parser({x,y,z}));
#endif
    }

    inline amrex::Real operator() (const amrex::RealArray& p) const noexcept {
        return this->operator()(AMREX_D_DECL(p[0],p[1],p[2]));
    }

private:
    ParserExecutor<3> m_parser;
};

}

#endif
