
#include <AMReX_BLassert.H>
#include <AMReX_ErrorList.H>
#include <AMReX_SPACE.H>

#include <iostream>

namespace amrex {

ErrorRec::ErrorFunc::ErrorFunc ()
    :
    m_func(nullptr),
    m_func3D(nullptr)
{}

ErrorRec::ErrorFunc::ErrorFunc (ErrorFuncDefault inFunc)
    :
    m_func(inFunc),
    m_func3D(nullptr)
{}

ErrorRec::ErrorFunc::ErrorFunc (ErrorFunc3DDefault inFunc)
    :
    m_func(nullptr),
    m_func3D(inFunc)
{}

ErrorRec::ErrorFunc*
ErrorRec::ErrorFunc::clone () const
{
    return new ErrorFunc(*this);
}

// \cond CODEGEN
void
ErrorRec::ErrorFunc::operator () (int* tag, AMREX_D_DECL(const int&tlo0,const int&tlo1,const int&tlo2),
                                  AMREX_D_DECL(const int&thi0,const int&thi1,const int&thi2),
                                  const int* tagval, const int* clearval,
                                  Real* data, AMREX_D_DECL(const int&dlo0,const int&dlo1,const int&dlo2),
                                  AMREX_D_DECL(const int&dhi0,const int&dhi1,const int&dhi2),
                                  const int* lo, const int * hi, const int* nvar,
                                  const int* domain_lo, const int* domain_hi,
                                  const Real* dx, const Real* xlo,
                                  const Real* prob_lo, const Real* time,
                                  const int* level) const
{
    BL_ASSERT(m_func != nullptr);

    m_func(tag,AMREX_D_DECL(tlo0,tlo1,tlo2),AMREX_D_DECL(thi0,thi1,thi2),
           tagval,clearval,data,AMREX_D_DECL(dlo0,dlo1,dlo2),AMREX_D_DECL(dhi0,dhi1,dhi2),lo,hi,nvar,
           domain_lo,domain_hi,dx,xlo,prob_lo,time,level);
}

void
ErrorRec::ErrorFunc::operator () (int* tag, const int* tlo, const int* thi,
                                  const int* tagval, const int* clearval,
                                  Real* data, const int* dlo, const int* dhi,
                                  const int* lo, const int * hi, const int* nvar,
                                  const int* domain_lo, const int* domain_hi,
                                  const Real* dx, const Real* xlo,
                                  const Real* prob_lo, const Real* time,
                                  const int* level) const
{
    BL_ASSERT(m_func3D != nullptr);

    m_func3D(tag,AMREX_ARLIM_3D(tlo),AMREX_ARLIM_3D(thi),
             tagval,clearval,data,AMREX_ARLIM_3D(dlo),AMREX_ARLIM_3D(dhi),
             AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),nvar,
             AMREX_ARLIM_3D(domain_lo),AMREX_ARLIM_3D(domain_hi),
             AMREX_ZFILL(dx),AMREX_ZFILL(xlo),AMREX_ZFILL(prob_lo),time,level);
}
// \endcond

ErrorRec::ErrorFunc2::ErrorFunc2 ()
    :
    m_func(nullptr)
{}

ErrorRec::ErrorFunc2::ErrorFunc2 (ErrorFunc2Default inFunc)
    :
    m_func(inFunc)
{}

ErrorRec::ErrorFunc2*
ErrorRec::ErrorFunc2::clone () const
{
    return new ErrorFunc2(*this);
}

void
ErrorRec::ErrorFunc2::operator () (int* tag, AMREX_D_DECL(const int&tlo0,const int&tlo1,const int&tlo2),
                                   AMREX_D_DECL(const int&thi0,const int&thi1,const int&thi2),
                                   const int* tagval, const int* clearval,
                                   Real* data, AMREX_D_DECL(const int&dlo0,const int&dlo1,const int&dlo2),
                                   AMREX_D_DECL(const int&dhi0,const int&dhi1,const int&dhi2),
                                   const int* lo, const int * hi, const int* nvar,
                                   const int* domain_lo, const int* domain_hi,
                                   const Real* dx, const int* level, const Real* avg) const
{
    BL_ASSERT(m_func != nullptr);

    m_func(tag,AMREX_D_DECL(tlo0,tlo1,tlo2),AMREX_D_DECL(thi0,thi1,thi2),
           tagval,clearval,data,AMREX_D_DECL(dlo0,dlo1,dlo2),AMREX_D_DECL(dhi0,dhi1,dhi2),lo,hi,nvar,
           domain_lo,domain_hi,dx,level,avg);
}


ErrorRec::ErrorRec (std::string nm, int ng, ErrorRec::ErrorType etyp,
                    const ErrorRec::ErrorFunc2& f2)
    :
    derive_name(std::move(nm)),
    ngrow(ng),
    err_type(etyp),
    err_func(nullptr),
    err_func2(f2.clone())
{}

ErrorRec::ErrorRec (std::string nm, int ng, ErrorRec::ErrorType etyp,
                    const ErrorRec::ErrorFunc& f)
    :
    derive_name(std::move(nm)),
    ngrow(ng),
    err_type(etyp),
    err_func(f.clone()),
    err_func2(nullptr)
{}

const std::string&
ErrorRec::name () const noexcept
{
    return derive_name;
}

int
ErrorRec::nGrow () const noexcept
{
    return ngrow;
}

ErrorRec::ErrorType
ErrorRec::errType () const noexcept
{
    return err_type;
}

const ErrorRec::ErrorFunc&
ErrorRec::errFunc () const
{
    return *err_func;
}

const ErrorRec::ErrorFunc2&
ErrorRec::errFunc2() const
{
    return *err_func2;
}

ErrorRec::~ErrorRec()
{
    delete err_func;
    delete err_func2;
}

int
ErrorList::size () const noexcept
{
    return static_cast<int>(vec.size());
}

void
ErrorList::add (const std::string&         name,
                int                        nextra,
                ErrorRec::ErrorType        typ,
                const ErrorRec::ErrorFunc& func)
{
    //
    // Keep list in order of definition, append().
    //
    auto n = vec.size();
    vec.resize(n+1);
    vec[n] = std::make_unique<ErrorRec>(name, nextra, typ, func);
}

void
ErrorList::add (const std::string&          name,
                int                         nextra,
                ErrorRec::ErrorType         typ,
                const ErrorRec::ErrorFunc2& func2)
{
    //
    // Keep list in order of definition, append().
    //
    auto n = vec.size();
    vec.resize(n+1);
    vec[n] = std::make_unique<ErrorRec>(name, nextra, typ, func2);
}

const ErrorRec&
ErrorList::operator[] (int k) const noexcept
{
    BL_ASSERT(k < size());

    return *vec[k];
}

static const char* err_name[] = { "Special", "Standard", "UseAverage" };

std::ostream&
operator << (std::ostream&    os,
             const ErrorList& elst)
{
    for (int i = 0; i < elst.size(); i++)
    {
        os << elst[i].name()
           << ' '
           << elst[i].nGrow()
           << ' '
           << err_name[elst[i].errType()]
           << '\n';
    }
    return os;
}

int
AMRErrorTag::SetNGrow () const noexcept
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_test != USER, "Do not call SetNGrow with USER test");
    static std::map<TEST,int> ng = { {GRAD,1}, {RELGRAD,1}, {LESS,0}, {GREATER,0}, {VORT,0}, {BOX,0} };
    return ng[m_test];
}

void
AMRErrorTag::operator() (TagBoxArray&    tba,
                         const MultiFab* mf,
                         char            clearval,
                         char            tagval,
                         Real            time,
                         int             level,
                         const Geometry& geom) const noexcept
{
    BL_PROFILE("AMRErrorTag::operator()");

    if (m_test == USER)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_userfunc!=nullptr,
                                         "UserFunc not properly set in AMRErrorTag");

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(tba,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const auto& bx    = mfi.tilebox();
            auto const& dat   = mf->array(mfi);
            auto tag          = tba.array(mfi);
            (*m_userfunc)(bx,dat,tag,time,level,tagval,clearval);
        }
    }
    else
    {
        if ((level <  m_info.m_max_level) &&
            (time  >= m_info.m_min_time ) &&
            (time  <= m_info.m_max_time ) )
        {
            auto const& tagma = tba.arrays();
            if (m_test == BOX)
            {
                const auto plo = geom.ProbLoArray();
                const auto dx  = geom.CellSizeArray();
                const auto tag_rb = m_info.m_realbox;
                ParallelFor(tba, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
                {
                    GpuArray<Real,AMREX_SPACEDIM> pt
                        {AMREX_D_DECL(plo[0]+(Real(i)+Real(0.5))*dx[0],
                                      plo[1]+(Real(j)+Real(0.5))*dx[1],
                                      plo[2]+(Real(k)+Real(0.5))*dx[2])};
                    if (tag_rb.contains(pt.data())) {
                        tagma[bi](i,j,k) = tagval;
                    }
                });
            }
            else
            {
                auto const& datma   = mf->const_arrays();
                auto threshold = m_value[level];
                auto const volume_weighting = m_info.m_volume_weighting;
                auto geomdata = geom.data();
                auto tag_update = tagval;
                if (m_info.m_derefine) {
                    tag_update = clearval;
                }

                if (m_test == GRAD)
                {
#ifdef AMREX_USE_EB
                    if (mf->hasEBFabFactory()) {
                        auto const& ebfact =
                            dynamic_cast<amrex::EBFArrayBoxFactory const&>(mf->Factory());
                        auto const& flags = ebfact.getMultiEBCellFlagFab().arrays();
                        ParallelFor(tba, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
                        {
                            auto const& dat  = datma[bi];
                            auto const& flag = flags[bi];

                            Real ax = 0.; Real ay = 0.;
                            if (flag(i,j,k).isConnected(1,0,0)) {
                                ax = amrex::max(ax,std::abs(dat(i+1,j,k) - dat(i,j,k)));
                            }
                            if (flag(i,j,k).isConnected(-1,0,0)) {
                                ax = amrex::max(ax,std::abs(dat(i,j,k) - dat(i-1,j,k)));
                            }
                            if (flag(i,j,k).isConnected(0,1,0)) {
                                ay = amrex::max(ay,std::abs(dat(i,j+1,k) - dat(i,j,k)));
                            }
                            if (flag(i,j,k).isConnected(0,-1,0)) {
                                ay = amrex::max(ay,std::abs(dat(i,j,k) - dat(i,j-1,k)));
                            }
#if AMREX_SPACEDIM > 2
                            Real az = 0.;
                            if (flag(i,j,k).isConnected(0,0,1)) {
                                az = amrex::max(az,std::abs(dat(i,j,k+1) - dat(i,j,k)));
                            }
                            if (flag(i,j,k).isConnected(0,0,-1)) {
                                az = amrex::max(az,std::abs(dat(i,j,k) - dat(i,j,k-1)));
                            }
#endif
                            if (amrex::max(AMREX_D_DECL(ax,ay,az)) >= threshold) {
                                tagma[bi](i,j,k) = tag_update;
                            }
                        });
                    } else
#endif
                    {
                        ParallelFor(tba, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
                        {
                            auto const& dat = datma[bi];

                            Real ax = 0.;
                            ax = std::abs(dat(i+1,j,k) - dat(i,j,k));
                            ax = amrex::max(ax,std::abs(dat(i,j,k) - dat(i-1,j,k)));
#if AMREX_SPACEDIM == 1
                            if (ax >= threshold) { tagma[bi](i,j,k) = tag_update;}
#else
                            Real ay = 0.;
                            ay = std::abs(dat(i,j+1,k) - dat(i,j,k));
                            ay = amrex::max(ay,std::abs(dat(i,j,k) - dat(i,j-1,k)));
#if AMREX_SPACEDIM > 2
                            Real az = 0.;
                            az = std::abs(dat(i,j,k+1) - dat(i,j,k));
                            az = amrex::max(az,std::abs(dat(i,j,k) - dat(i,j,k-1)));
#endif // DIM > 2
                            if (amrex::max(AMREX_D_DECL(ax,ay,az)) >= threshold) {
                                tagma[bi](i,j,k) = tag_update;
                            }
#endif // DIM > 1
                       });
                    }
                }
                else if (m_test == RELGRAD)
                {
#ifdef AMREX_USE_EB
                    if (mf->hasEBFabFactory()) {
                        auto const& ebfact =
                            dynamic_cast<amrex::EBFArrayBoxFactory const&>(mf->Factory());
                        auto const& flags = ebfact.getMultiEBCellFlagFab().arrays();
                        ParallelFor(tba, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
                        {
                            auto const& dat  = datma[bi];
                            auto const& flag = flags[bi];

                            Real ax = 0.; Real ay = 0.;

                            if (flag(i,j,k).isConnected(1,0,0)) {
                                ax = amrex::max(ax,std::abs(dat(i+1,j,k) - dat(i,j,k)));
                            }
                            if (flag(i,j,k).isConnected(-1,0,0)) {
                                ax = amrex::max(ax,std::abs(dat(i,j,k) - dat(i-1,j,k)));
                            }
                            if (flag(i,j,k).isConnected(0,1,0)) {
                                ay = amrex::max(ay,std::abs(dat(i,j+1,k) - dat(i,j,k)));
                            }
                            if (flag(i,j,k).isConnected(0,-1,0)) {
                                ay = amrex::max(ay,std::abs(dat(i,j,k) - dat(i,j-1,k)));
                            }
#if AMREX_SPACEDIM > 2
                            Real az = 0.;
                            if (flag(i,j,k).isConnected(0,0,1)) {
                                az = amrex::max(az,std::abs(dat(i,j,k+1) - dat(i,j,k)));
                            }
                            if (flag(i,j,k).isConnected(0,0,-1)) {
                                az = amrex::max(az,std::abs(dat(i,j,k) - dat(i,j,k-1)));
                            }
#endif // DIM > 2
                            if (amrex::max(AMREX_D_DECL(ax,ay,az))
                                >= threshold * std::abs(dat(i,j,k))) {
                                tagma[bi](i,j,k) = tag_update;
                            }
                        });
                    } else
#endif
                    {
                        ParallelFor(tba, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
                        {
                            auto const& dat = datma[bi];

                            Real ax = std::abs(dat(i+1,j,k) - dat(i,j,k));
                            ax = amrex::max(ax,std::abs(dat(i,j,k) - dat(i-1,j,k)));
#if AMREX_SPACEDIM == 1
                            if (ax >= threshold * std::abs(dat(i,j,k))) { tagma[bi](i,j,k) = tag_update;}
#else
                            Real ay = std::abs(dat(i,j+1,k) - dat(i,j,k));
                            ay = amrex::max(ay,std::abs(dat(i,j,k) - dat(i,j-1,k)));
#if AMREX_SPACEDIM > 2
                            Real az = std::abs(dat(i,j,k+1) - dat(i,j,k));
                            az = amrex::max(az,std::abs(dat(i,j,k) - dat(i,j,k-1)));
#endif // DIM > 2
                            if (amrex::max(AMREX_D_DECL(ax,ay,az))
                                >= threshold * std::abs(dat(i,j,k))) {
                                tagma[bi](i,j,k) = tag_update;
                            }
#endif // DIM > 1
                        });
                    }
                }
                else if (m_test == LESS)
                {
#ifdef AMREX_USE_EB
                    if (mf->hasEBFabFactory()) {
                        auto const& ebfact =
                            dynamic_cast<amrex::EBFArrayBoxFactory const&>(mf->Factory());
                        auto const& flags = ebfact.getMultiEBCellFlagFab().arrays();
                        ParallelFor(tba, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
                        {
                            Real vol = volume_weighting ? Geometry::Volume(IntVect{AMREX_D_DECL(i,j,k)}, geomdata) : 1.0_rt;
                            auto const& flag = flags[bi];
                            if (!flag(i,j,k).isCovered()) {
                                if (datma[bi](i,j,k) * vol <= threshold) {
                                    tagma[bi](i,j,k) = tag_update;
                                }
                            }
                        });
                    } else
#endif
                    {
                    ParallelFor(tba, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept {
                        Real vol = volume_weighting ? Geometry::Volume(IntVect{AMREX_D_DECL(i,j,k)}, geomdata) : 1.0_rt;
                        if (datma[bi](i,j,k) * vol <= threshold) {
                            tagma[bi](i,j,k) = tag_update;
                        }
                    });
                    }
                }
                else if (m_test == GREATER)
                {
#ifdef AMREX_USE_EB
                    if (mf->hasEBFabFactory()) {
                        auto const& ebfact =
                            dynamic_cast<amrex::EBFArrayBoxFactory const&>(mf->Factory());
                        auto const& flags = ebfact.getMultiEBCellFlagFab().arrays();
                        ParallelFor(tba, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
                        {
                            Real vol = volume_weighting ? Geometry::Volume(IntVect{AMREX_D_DECL(i,j,k)}, geomdata) : 1.0_rt;
                            auto const& flag = flags[bi];
                            if (!flag(i,j,k).isCovered()) {
                                if (datma[bi](i,j,k) * vol >= threshold) {
                                    tagma[bi](i,j,k) = tag_update;
                                }
                            }
                        });
                    } else
#endif
                    ParallelFor(tba, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
                    {
                        Real vol = volume_weighting ? Geometry::Volume(IntVect{AMREX_D_DECL(i,j,k)}, geomdata) : 1.0_rt;
                            if (datma[bi](i,j,k) * vol >= threshold) {
                                tagma[bi](i,j,k) = tag_update;
                            }
                    });
                }
                else if (m_test == VORT)
                {
                    const Real fac = threshold * Real(std::pow(2,level));
#ifdef AMREX_USE_EB
                    if (mf->hasEBFabFactory()) {
                        auto const& ebfact =
                            dynamic_cast<amrex::EBFArrayBoxFactory const&>(mf->Factory());
                        auto const& flags = ebfact.getMultiEBCellFlagFab().arrays();
                        ParallelFor(tba, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
                        {
                            auto const& flag = flags[bi];
                            if (!flag(i,j,k).isCovered()) {
                                if (datma[bi](i,j,k) >= fac) {
                                    tagma[bi](i,j,k) = tag_update;
                                }
                            }
                        });
                    } else
#endif
                    {
                        ParallelFor(tba, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
                        {
                            if (datma[bi](i,j,k) >= fac) {
                                tagma[bi](i,j,k) = tag_update;
                            }
                        });
                    }
                }
                else
                {
                    Abort("Bad AMRErrorTag test flag");
                }
            }
            Gpu::streamSynchronize();
        }
    }
}

}
