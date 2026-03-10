
#include <AMReX_StateDescriptor.H>
#include <AMReX_BCRec.H>

#include <algorithm>
#include <string>
#include <iostream>

namespace amrex {

int StateDescriptor::bf_ext_dir_threadsafe = 0;

bool
StateDescriptor::bf_thread_safety (const int* /*lo*/,const int* /*hi*/,
                                   const int* /*dom_lo*/, const int* /*dom_hi*/,
                                   const int* bc, int ng)
{
    bool thread_safe = true;
    if (!bf_ext_dir_threadsafe) {
        bool has_ext_dir = false;
        for (int i=0; i<2*AMREX_SPACEDIM*ng && !has_ext_dir; ++i) {
            has_ext_dir = ((bc[i]==BCType::ext_dir) || (bc[i]==BCType::ext_dir_cc));
        }
        if (has_ext_dir) { thread_safe = false; }
    }
    return thread_safe;
}

void
StateDescriptor::BndryFunc::operator () (Real* data,const int* lo,const int* hi,
                                         const int* dom_lo, const int* dom_hi,
                                         const Real* dx, const Real* grd_lo,
                                         const Real* time, const int* a_bc) const
{
    BL_ASSERT(m_func != nullptr || m_func3D != nullptr);

#ifdef AMREX_USE_OMP
    bool thread_safe = bf_thread_safety(lo, hi, dom_lo, dom_hi, a_bc, 1);
    if (thread_safe) {
#endif
        {
            if (m_func != nullptr) {
                m_func(data,AMREX_ARLIM(lo),AMREX_ARLIM(hi),dom_lo,dom_hi,dx,grd_lo,time,a_bc);
            } else {
                m_func3D(data,AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),AMREX_ARLIM_3D(dom_lo),AMREX_ARLIM_3D(dom_hi),
                         AMREX_ZFILL(dx),AMREX_ZFILL(grd_lo),time,a_bc);
            }
        }
#ifdef AMREX_USE_OMP
    } else {
#pragma omp critical (bndryfunc)
        {
            if (m_func != nullptr) {
                m_func(data,AMREX_ARLIM(lo),AMREX_ARLIM(hi),dom_lo,dom_hi,dx,grd_lo,time,a_bc);
            } else {
                m_func3D(data,AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),AMREX_ARLIM_3D(dom_lo),AMREX_ARLIM_3D(dom_hi),
                         AMREX_ZFILL(dx),AMREX_ZFILL(grd_lo),time,a_bc);
            }
        }
    }
#endif
}

void
StateDescriptor::BndryFunc::operator () (Real* data,const int* lo,const int* hi,
                                         const int* dom_lo, const int* dom_hi,
                                         const Real* dx, const Real* grd_lo,
                                         const Real* time, const int* a_bc, int ng) const
{
    BL_ASSERT(m_gfunc != nullptr || m_gfunc3D != nullptr);

    amrex::ignore_unused(ng);
#ifdef AMREX_USE_OMP
    bool thread_safe = bf_thread_safety(lo, hi, dom_lo, dom_hi, a_bc, ng);
    if (thread_safe) {
#endif
        {
            if (m_gfunc != nullptr) {
                m_gfunc(data,AMREX_ARLIM(lo),AMREX_ARLIM(hi),dom_lo,dom_hi,dx,grd_lo,time,a_bc);
            } else {
                m_gfunc3D(data,AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),AMREX_ARLIM_3D(dom_lo),AMREX_ARLIM_3D(dom_hi),
                          AMREX_ZFILL(dx),AMREX_ZFILL(grd_lo),time,a_bc);
            }
        }
#ifdef AMREX_USE_OMP
    } else {
#pragma omp critical (bndryfunc)
        {
            if (m_gfunc != nullptr) {
                m_gfunc(data,AMREX_ARLIM(lo),AMREX_ARLIM(hi),dom_lo,dom_hi,dx,grd_lo,time,a_bc);
            } else {
                m_gfunc3D(data,AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),AMREX_ARLIM_3D(dom_lo),AMREX_ARLIM_3D(dom_hi),
                          AMREX_ZFILL(dx),AMREX_ZFILL(grd_lo),time,a_bc);
            }
        }
    }
#endif
}

void
StateDescriptor::BndryFunc::operator () (Box const& bx, FArrayBox& data,
                                         const int dcomp, const int numcomp,
                                         Geometry const& geom, const Real time,
                                         const Vector<BCRec>& bcr, const int bcomp,
                                         const int scomp) const
{
    AMREX_ASSERT(m_funcfab != nullptr);
    m_funcfab(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
}

void
DescriptorList::clear ()
{
    desc.clear();
}

int
DescriptorList::size () const noexcept
{
    return static_cast<int>(desc.size());
}

void
DescriptorList::resetComponentBCs (int                               indx,
                                   int                               comp,
                                   const BCRec&                      bc,
                                   const StateDescriptor::BndryFunc& func)
{
    desc[indx]->resetComponentBCs(comp,bc,func);
}

void
DescriptorList::setComponent (int                               indx,
                              int                               comp,
                              const std::string&                nm,
                              const BCRec&                      bc,
                              const StateDescriptor::BndryFunc& func,
                              InterpBase*                       interp,
                              int                               max_map_start_comp,
                              int                               min_map_end_comp)
{
    desc[indx]->setComponent(comp,nm,bc,func,interp,max_map_start_comp,min_map_end_comp);
}

void
DescriptorList::setComponent (int                               indx,
                              int                               comp,
                              const Vector<std::string>&         nm,
                              const Vector<BCRec>&               bc,
                              const StateDescriptor::BndryFunc& func,
                              InterpBase*                       interp)
{
    BL_ASSERT(nm.size() == bc.size());
    BL_ASSERT(comp >= 0);
    for (int i = 0; i < nm.size(); i++)
    {
        const bool is_primary = (i == 0) ? true : false;

        desc[indx]->setComponent(comp+i,nm[i],bc[i],func,interp,is_primary,
                                 static_cast<int>(nm.size()));
    }
}

const StateDescriptor&
DescriptorList::operator[] (int k) const noexcept
{
    return *desc[k];
}

void
DescriptorList::addDescriptor (int                         indx,
                               IndexType                   typ,
                               StateDescriptor::TimeCenter ttyp,
                               int                         nextra,
                               int                         num_comp,
                               InterpBase*                 interp,
                               bool                        extrap,
                               bool                        a_store_in_checkpoint)
{
    if (indx >= desc.size()) {
        desc.resize(indx+1);
    }
    desc[indx] = std::make_unique<StateDescriptor>(typ,ttyp,indx,nextra,num_comp,interp,extrap,
                                                   a_store_in_checkpoint);
}

StateDescriptor::StateDescriptor (IndexType                   btyp,
                                  StateDescriptor::TimeCenter ttyp,
                                  int                         ident,
                                  int                         nextra,
                                  int                         num_comp,
                                  InterpBase*                 a_interp,
                                  bool                        a_extrap,
                                  bool                        a_store_in_checkpoint)
    :
    type(btyp),
    t_type(ttyp),
    id(ident),
    ncomp(num_comp),
    ngrow(nextra),
    mapper(a_interp),
    m_extrap(a_extrap),
    m_store_in_checkpoint(a_store_in_checkpoint)
{
    BL_ASSERT (num_comp > 0);

    names.resize(num_comp);
    bc.resize(num_comp);
    bc_func.resize(num_comp);
    mapper_comp.resize(num_comp);
    m_primary.resize(num_comp);
    m_groupsize.resize(num_comp);
    max_map_start_comp.resize(num_comp);
    min_map_end_comp.resize(num_comp);
}

void
StateDescriptor::resetComponentBCs (int              comp,
                                    const BCRec&     bcr,
                                    const BndryFunc& func)
{
    BL_ASSERT(comp >= 0 && comp < ncomp);

    bc_func[comp] = std::make_unique<BndryFunc>(func);
    bc[comp] = bcr;
}

IndexType
StateDescriptor::getType () const noexcept
{
    return type;
}

StateDescriptor::TimeCenter
StateDescriptor::timeType () const noexcept
{
    return t_type;
}

int
StateDescriptor::nComp () const noexcept
{
    return ncomp;
}

int
StateDescriptor::nExtra () const noexcept
{
    return ngrow;
}

InterpBase*
StateDescriptor::interp () const noexcept
{
    return mapper;
}

InterpBase*
StateDescriptor::interp (int i) const noexcept
{
    return mapper_comp[i] == 0 ? mapper : mapper_comp[i];
}

const std::string&
StateDescriptor::name (int i) const noexcept
{
    return names[i];
}

const BCRec&
StateDescriptor::getBC (int i) const noexcept
{
    return bc[i];
}

const Vector<BCRec>&
StateDescriptor::getBCs () const noexcept
{
    return bc;
}

bool
StateDescriptor::extrap () const noexcept
{
    return m_extrap;
}

bool
StateDescriptor::store_in_checkpoint () const noexcept
{
    return m_store_in_checkpoint;
}


const StateDescriptor::BndryFunc&
StateDescriptor::bndryFill (int i) const noexcept
{
    return *bc_func[i];
}

int
StateDescriptor::inRange (int sc, int nc) const noexcept
{
    return sc>=0 && sc+nc<=ncomp;
}

void
StateDescriptor::define (IndexType                   btyp,
                         StateDescriptor::TimeCenter ttyp,
                         int                         ident,
                         int                         nextra,
                         int                         num_comp,
                         InterpBase*                 a_interp,
                         bool                        a_extrap,
                         bool                        a_store_in_checkpoint)
{
    type     = btyp;
    t_type   = ttyp;
    id       = ident;
    ngrow    = nextra;
    ncomp    = num_comp;
    mapper   = a_interp;
    m_extrap = a_extrap;
    m_store_in_checkpoint = a_store_in_checkpoint;

    BL_ASSERT (num_comp > 0);

    names.resize(num_comp);
    bc.resize(num_comp);
    bc_func.resize(num_comp);
    mapper_comp.resize(num_comp);
    m_primary.resize(num_comp);
    m_groupsize.resize(num_comp);
    max_map_start_comp.resize(num_comp);
    min_map_end_comp.resize(num_comp);
}

void
StateDescriptor::setComponent (int                               comp,
                               const std::string&                nm,
                               const BCRec&                      bcr,
                               const StateDescriptor::BndryFunc& func,
                               InterpBase*                       a_interp,
                               int                               max_map_start_comp_,
                               int                               min_map_end_comp_)
{
    BL_ASSERT(comp >= 0 && comp < ncomp);
    bc_func[comp] = std::make_unique<BndryFunc>(func);

    names[comp]       = nm;
    bc[comp]          = bcr;
    mapper_comp[comp] = a_interp;
    m_primary[comp]    = false;
    m_groupsize[comp] = 0;

    if (max_map_start_comp_>=0 && min_map_end_comp_>=0)
    {
        BL_ASSERT(comp >= max_map_start_comp_ &&
                  comp <= min_map_end_comp_   &&
                  min_map_end_comp_ < ncomp);
        max_map_start_comp[comp] = max_map_start_comp_;
        min_map_end_comp[comp]   = min_map_end_comp_;
    }
    else
    {
        max_map_start_comp[comp] = comp;
        min_map_end_comp[comp]   = comp;
    }
}


void
StateDescriptor::setComponent (int                               comp,
                               const std::string&                nm,
                               const BCRec&                      bcr,
                               const StateDescriptor::BndryFunc& func,
                               InterpBase*                       a_interp,
                               bool                              a_primary,
                               int                               a_groupsize)
{
    setComponent(comp,nm,bcr,func,a_interp,-1,-1);

    m_primary[comp]    = a_primary;
    m_groupsize[comp] = a_groupsize;
}

void
StateDescriptor::dumpNames (std::ostream& os,
                            int           start_comp,
                            int           num_comp) const
{
    BL_ASSERT(start_comp >= 0 && start_comp+num_comp <= ncomp);

    for (int k = 0; k < num_comp; k++)
    {
        os << names[start_comp+k] << ' ';
    }
}

bool
StateDescriptor::identicalInterps (int a_scomp,
                                   int a_ncomp) const noexcept
{
    BL_ASSERT(a_scomp >= 0);
    BL_ASSERT(a_ncomp >= 1);

    InterpBase* map = interp(a_scomp);

    for (int i = a_scomp+1; i < a_scomp+a_ncomp; i++) {
        if (!(map == interp(i))) {
            return false;
        }
    }

    return true;
}

std::vector< std::pair<int,int> >
StateDescriptor::sameInterps (int a_scomp,
                              int a_ncomp) const
{
    BL_ASSERT(a_scomp >= 0);
    BL_ASSERT(a_ncomp >= 1);

    std::vector< std::pair<int,int> > range;

    InterpBase* map = interp(a_scomp);

    int SComp = a_scomp, NComp = 1;

    for (int i = a_scomp+1; i < a_scomp+a_ncomp; i++)
    {
        if (map == interp(i))
        {
            NComp++;
        }
        else
        {
            range.emplace_back(SComp,NComp);

            map   = interp(i);
            SComp = i;
            NComp = 1;
        }
    }

    range.emplace_back(SComp,NComp);

#ifdef AMREX_DEBUG
    int sum = 0;
    for (auto const& r : range) {
        sum += r.second;
    }
    BL_ASSERT(sum == a_ncomp);
#endif

    return range;
}

}
