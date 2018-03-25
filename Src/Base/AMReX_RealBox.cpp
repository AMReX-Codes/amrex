
#include <iostream>
#include <string>

#include <AMReX_RealBox.H>
#include <AMReX_CArena.H>

namespace amrex {

#ifdef AMREX_USE_CUDA
int RealBox_init::m_cnt = 0;

namespace
{
    Arena* the_realbox_arena = 0;
}

RealBox_init::RealBox_init ()
{
    if (m_cnt++ == 0)
    {
        BL_ASSERT(the_realbox_arena == 0);

        const std::size_t hunk_size = 64 * 1024;

        the_realbox_arena = new CArena(hunk_size);

	the_realbox_arena->SetHostAlloc();
    }
}

RealBox_init::~RealBox_init ()
{
    if (--m_cnt == 0)
        delete the_realbox_arena;
}

Arena*
The_RealBox_Arena ()
{
    BL_ASSERT(the_realbox_arena != 0);

    return the_realbox_arena;
}
#endif

RealBox::RealBox (const Box&  bx,
                  const Real* dx,
                  const Real* base)
{
    const int* blo = bx.loVect();
    const int* bhi = bx.hiVect();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        xlo[i] = base[i] + dx[i]*blo[i];
        int shft = (bx.type(i) == IndexType::CELL ? 1 : 0);
        xhi[i] = base[i] + dx[i]*(bhi[i]+ shft);
    }
    nullify_device_memory();
}

RealBox::RealBox ()
{
    AMREX_D_TERM(xlo[0] , = xlo[1] , = xlo[2] ) = 0.;
    AMREX_D_TERM(xhi[0] , = xhi[1] , = xhi[2] ) = -1.;
    nullify_device_memory();
}

RealBox::RealBox (const Real* a_lo,
                  const Real* a_hi)
{
    AMREX_D_EXPR(xlo[0] = a_lo[0] , xlo[1] = a_lo[1] , xlo[2] = a_lo[2]);
    AMREX_D_EXPR(xhi[0] = a_hi[0] , xhi[1] = a_hi[1] , xhi[2] = a_hi[2]);
    nullify_device_memory();
}

RealBox::RealBox (const std::array<Real,AMREX_SPACEDIM>& a_lo,
                  const std::array<Real,AMREX_SPACEDIM>& a_hi)
{
    AMREX_D_EXPR(xlo[0] = a_lo[0] , xlo[1] = a_lo[1] , xlo[2] = a_lo[2]);
    AMREX_D_EXPR(xhi[0] = a_hi[0] , xhi[1] = a_hi[1] , xhi[2] = a_hi[2]);
    nullify_device_memory();
}

RealBox::RealBox (AMREX_D_DECL(Real x0, Real y0, Real z0),
                  AMREX_D_DECL(Real x1, Real y1, Real z1))
{
    AMREX_D_EXPR(xlo[0] = x0 , xlo[1] = y0 , xlo[2] = z0);
    AMREX_D_EXPR(xhi[0] = x1 , xhi[1] = y1 , xhi[2] = z1);
    nullify_device_memory();
}

void
RealBox::nullify_device_memory() const
{
#ifdef AMREX_USE_CUDA
    xlo_d = nullptr;
    xhi_d = nullptr;
#endif
}

void
RealBox::initialize_device_memory() const
{
#ifdef AMREX_USE_CUDA
    initialize_lo();
    initialize_hi();
#endif
}

void
RealBox::initialize_lo() const
{
#ifdef AMREX_USE_CUDA
    const size_t sz = 3 * sizeof(Real);

    Real* xlo_temp = static_cast<Real*>(amrex::The_RealBox_Arena()->alloc(sz));
    xlo_d.reset(xlo_temp, [](Real* ptr) { amrex::The_RealBox_Arena()->free(ptr); });
    copy_xlo();
#endif
}

void
RealBox::initialize_hi() const
{
#ifdef AMREX_USE_CUDA
    const size_t sz = 3 * sizeof(Real);

    Real* xhi_temp = static_cast<Real*>(amrex::The_RealBox_Arena()->alloc(sz));
    xhi_d.reset(xhi_temp, [](Real* ptr) { amrex::The_RealBox_Arena()->free(ptr); });
    copy_xhi();
#endif
}

void
RealBox::copy_device_memory() const
{
  copy_xlo();
  copy_xhi();
}

void
RealBox::copy_xlo() const
{
#ifdef AMREX_USE_CUDA
    for (int i = 0; i < BL_SPACEDIM; ++i)
	xlo_d.get()[i] = xlo[i];
    for (int i = BL_SPACEDIM; i < 3; ++i)
	xlo_d.get()[i] = 0;
#endif
}

void
RealBox::copy_xhi() const
{
#ifdef AMREX_USE_CUDA
    for (int i = 0; i < BL_SPACEDIM; ++i)
	xhi_d.get()[i] = xhi[i];
    for (int i = BL_SPACEDIM; i < 3; ++i)
	xhi_d.get()[i] = 0;
#endif
}

bool
RealBox::contains (const RealBox& rb, Real eps) const
{
    return contains(rb.xlo, eps) && contains(rb.xhi, eps);
}

bool
RealBox::ok () const
{
    return (length(0) >= 0.0)
#if (AMREX_SPACEDIM > 1)
        && (length(1) >= 0.0)
#endif   
#if (AMREX_SPACEDIM > 2)
        && (length(2) >= 0.0)
#endif
   ;
}

Real
RealBox::volume () const
{
    if (ok()) return AMREX_D_TERM(length(0), *length(1), *length(2));
    return 0.0;
}

bool
RealBox::contains (const Real* point, Real eps) const
{
    return  AMREX_D_TERM((xlo[0]-eps < point[0]) && (point[0] < xhi[0]+eps),
                   && (xlo[1]-eps < point[1]) && (point[1] < xhi[1]+eps),
                   && (xlo[2]-eps < point[2]) && (point[2] < xhi[2]+eps));
}

std::ostream&
operator << (std::ostream &os, const RealBox& b)
{
    os << "(RealBox ";
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        os << b.lo(i) << ' ' << b.hi(i) << ' ';
    os << ')';
    return os;
}

//
// Copied from <Utility.H>
//
#define BL_IGNORE_MAX 100000

std::istream&
operator >> (std::istream &is, RealBox& b)
{
    is.ignore(BL_IGNORE_MAX,'(');

    std::string s;

    is >> s;

    if (s != "RealBox")
    {
        std::cerr << "unexpected token in RealBox: " << s << '\n';
        amrex::Abort();
    }

    Real lo[AMREX_SPACEDIM];
    Real hi[AMREX_SPACEDIM];
#ifdef BL_USE_FLOAT
    double dlotemp, dhitemp;
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        is >> dlotemp >> dhitemp;
        lo[i] = dlotemp;
        hi[i] = dhitemp;
    }
#else
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        is >> lo[i] >> hi[i];
#endif

    is.ignore(BL_IGNORE_MAX, ')');

    b = RealBox(lo,hi);

    return is;
}

}
