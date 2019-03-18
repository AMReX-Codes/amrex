
#include <algorithm>
#include <string>
#include <iostream>

#include <AMReX_StateDescriptor.H>
#include <AMReX_Interpolater.H>
#include <AMReX_BCRec.H>

namespace amrex {

int StateDescriptor::bf_ext_dir_threadsafe = 0;

StateDescriptor::BndryFunc*
StateDescriptor::BndryFunc::clone () const
{
    return new BndryFunc(*this);
}

StateDescriptor::BndryFunc::~BndryFunc () {}

bool
StateDescriptor::bf_thread_safety (const int* lo,const int* hi,
				   const int* dom_lo, const int* dom_hi,
				   const int* bc, int ng)
{
    bool thread_safe = true;
    if (!bf_ext_dir_threadsafe) {
	bool has_ext_dir = false;
	for (int i=0; i<2*AMREX_SPACEDIM*ng && !has_ext_dir; ++i) {
	    has_ext_dir = bc[i]==BCType::ext_dir;
	}
	if (has_ext_dir) thread_safe = false;
    }
    return thread_safe;
}

void
StateDescriptor::BndryFunc::operator () (Real* data,const int* lo,const int* hi,
                                         const int* dom_lo, const int* dom_hi,
                                         const Real* dx, const Real* grd_lo,
                                         const Real* time, const int* bc) const
{
    BL_ASSERT(m_func != 0 || m_func3D != 0);

    bool thread_safe = bf_thread_safety(lo, hi, dom_lo, dom_hi, bc, 1);
    if (thread_safe) {
      if (m_func != 0)
	m_func(data,AMREX_ARLIM(lo),AMREX_ARLIM(hi),dom_lo,dom_hi,dx,grd_lo,time,bc);
      else
	m_func3D(data,AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),AMREX_ARLIM_3D(dom_lo),AMREX_ARLIM_3D(dom_hi),
                 AMREX_ZFILL(dx),AMREX_ZFILL(grd_lo),time,bc);
    } else {
#ifdef _OPENMP
#pragma omp critical (bndryfunc)
#endif
      if (m_func != 0)
	m_func(data,AMREX_ARLIM(lo),AMREX_ARLIM(hi),dom_lo,dom_hi,dx,grd_lo,time,bc);
      else
	m_func3D(data,AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),AMREX_ARLIM_3D(dom_lo),AMREX_ARLIM_3D(dom_hi),
                 AMREX_ZFILL(dx),AMREX_ZFILL(grd_lo),time,bc);
    }
}

void
StateDescriptor::BndryFunc::operator () (Real* data,const int* lo,const int* hi,
                                         const int* dom_lo, const int* dom_hi,
                                         const Real* dx, const Real* grd_lo,
                                         const Real* time, const int* bc, int ng) const
{
    BL_ASSERT(m_gfunc != 0 || m_gfunc3D != 0);

    bool thread_safe = bf_thread_safety(lo, hi, dom_lo, dom_hi, bc, ng);
    if (thread_safe) {
        if (m_gfunc != 0)
	  m_gfunc(data,AMREX_ARLIM(lo),AMREX_ARLIM(hi),dom_lo,dom_hi,dx,grd_lo,time,bc);
	else
	  m_gfunc3D(data,AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),AMREX_ARLIM_3D(dom_lo),AMREX_ARLIM_3D(dom_hi),
                    AMREX_ZFILL(dx),AMREX_ZFILL(grd_lo),time,bc);
    } else {
#ifdef _OPENMP
#pragma omp critical (bndryfunc)
#endif
        if (m_gfunc != 0)
	  m_gfunc(data,AMREX_ARLIM(lo),AMREX_ARLIM(hi),dom_lo,dom_hi,dx,grd_lo,time,bc);
	else
	  m_gfunc3D(data,AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),AMREX_ARLIM_3D(dom_lo),AMREX_ARLIM_3D(dom_hi),
                    AMREX_ZFILL(dx),AMREX_ZFILL(grd_lo),time,bc);
    }
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

DescriptorList::DescriptorList () noexcept
{}

void
DescriptorList::clear ()
{
    desc.clear();
}

int
DescriptorList::size () const noexcept
{
    return desc.size();
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
                              Interpolater*                     interp,
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
                              Interpolater*                     interp)
{
    for (int i = 0; i < nm.size(); i++)
    {
        const bool master = (i == 0) ? true : false;

        desc[indx]->setComponent(comp+i,nm[i],bc[i],func,interp,master,nm.size());
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
                               Interpolater*               interp,
                               bool                        extrap,
                               bool                        a_store_in_checkpoint)
{
    if (indx >= desc.size())
        desc.resize(indx+1);
    desc[indx].reset(new StateDescriptor(typ,ttyp,indx,nextra,num_comp,interp,extrap,a_store_in_checkpoint));
}  


StateDescriptor::StateDescriptor () noexcept
    :
    t_type(Point),
    id(-1),
    ncomp(0),
    ngrow(0),
    mapper(0),
    m_extrap(false),
    m_store_in_checkpoint(true)
{}

StateDescriptor::StateDescriptor (IndexType                   btyp,
                                  StateDescriptor::TimeCenter ttyp,
                                  int                         ident,
                                  int                         nextra, 
                                  int                         num_comp,
                                  Interpolater*               a_interp,
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
    m_master.resize(num_comp);
    m_groupsize.resize(num_comp);
    max_map_start_comp.resize(num_comp);
    min_map_end_comp.resize(num_comp);
}

StateDescriptor::~StateDescriptor ()
{
    mapper = 0;
}

void
StateDescriptor::resetComponentBCs (int              comp,
                                    const BCRec&     bcr,
                                    const BndryFunc& func)
{
    BL_ASSERT(comp >= 0 && comp < ncomp);

    bc_func[comp].reset(func.clone());
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

Interpolater*
StateDescriptor::interp () const noexcept
{
    return mapper;
}

Interpolater*
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
                         Interpolater*               a_interp,
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
    m_master.resize(num_comp);
    m_groupsize.resize(num_comp);
    max_map_start_comp.resize(num_comp);
    min_map_end_comp.resize(num_comp);
}

void
StateDescriptor::setComponent (int                               comp,
                               const std::string&                nm,
                               const BCRec&                      bcr,
                               const StateDescriptor::BndryFunc& func,
                               Interpolater*                     a_interp, 
                               int                               max_map_start_comp_,
                               int                               min_map_end_comp_)
{
    bc_func[comp].reset(func.clone());

    names[comp]       = nm;
    bc[comp]          = bcr;
    mapper_comp[comp] = a_interp;
    m_master[comp]    = false;
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
                               Interpolater*                     a_interp,
                               bool                              a_master,
                               int                               a_groupsize)
{
    setComponent(comp,nm,bcr,func,a_interp,-1,-1);

    m_master[comp]    = a_master;
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

void
StateDescriptor::setUpMaps (int&                use_default_map,
                            const Interpolater* default_map,
                            int                 start_comp,
                            int                 num_comp,
                            Interpolater**&     maps, 
                            int&                nmaps,
                            int*&               map_start_comp,
                            int*&               map_num_comp, 
                            int*&               max_start_comp,
                            int*&               min_end_comp) const
{
    BL_ASSERT(start_comp>=0 && start_comp+num_comp-1 < ncomp && num_comp>0);

    maps           = 0;
    map_start_comp = 0;
    map_num_comp   = 0;
    max_start_comp = 0;
    min_end_comp   = 0;
    //
    // First, count number of interpolaters needed and allocate.
    //
    Interpolater* map = mapper_comp[start_comp];
    if (!map) map = (Interpolater*) default_map;
    nmaps = 1; 
    int icomp = start_comp+1;

    use_default_map = 1;
    while (icomp < start_comp+num_comp)
    {
        Interpolater* mapper_icomp = mapper_comp[icomp];
        if (!mapper_icomp)
        {
            mapper_icomp = (Interpolater *) default_map;
        }
        else
        {
            use_default_map = 0;
        }
        if (map != mapper_icomp)
        {
            map = mapper_icomp;
            nmaps++;
        }
        icomp++;
    }

    if (use_default_map) return;

    maps           = new Interpolater*[nmaps];
    map_start_comp = new int[nmaps];
    map_num_comp   = new int[nmaps];
    min_end_comp   = new int[nmaps];
    max_start_comp = new int[nmaps];
    //
    // Now fill the slots.
    //
    int imap             = 0;

    if (mapper_comp[start_comp])
    {
        maps[imap] = mapper_comp[start_comp];
    }
    else
    {
        maps[imap] = (Interpolater *) default_map;
    }

    icomp                = start_comp+1;
    map_start_comp[imap] = start_comp;
    map_num_comp[imap]   = 1;
    min_end_comp[imap]   = min_map_end_comp[start_comp];
    max_start_comp[imap] = max_map_start_comp[start_comp];

    while (icomp < start_comp+num_comp)
    {
        Interpolater* mapper_icomp = mapper_comp[icomp];

        if (!mapper_icomp)
            mapper_icomp = (Interpolater *) default_map;

        if (maps[imap] != mapper_icomp)
        {
            imap++;

            BL_ASSERT (imap < nmaps);

            maps[imap]           = mapper_icomp;
            map_start_comp[imap] = icomp;
            map_num_comp[imap]   = 1;
            min_end_comp[imap]   = min_map_end_comp[icomp];
            max_start_comp[imap] = max_map_start_comp[icomp];

        }
        else
        {
            map_num_comp[imap]++;
            min_end_comp[imap]   = std::max(min_end_comp[imap],min_map_end_comp[icomp]);
            max_start_comp[imap] = std::min(max_start_comp[imap],max_map_start_comp[icomp]);
        }
        icomp++;
    }
}

void
StateDescriptor::cleanUpMaps (Interpolater**& maps, 
                              int*&           map_start_comp,
                              int*&           map_num_comp,
                              int*&           max_start_comp,
                              int*&           min_end_comp) const
{
    delete [] maps;
    delete [] map_start_comp;
    delete [] map_num_comp;
    delete [] max_start_comp;
    delete [] min_end_comp;
}

bool
StateDescriptor::identicalInterps (int a_scomp,
                                   int a_ncomp) const noexcept
{
    BL_ASSERT(a_scomp >= 0);
    BL_ASSERT(a_ncomp >= 1);

    Interpolater* map = interp(a_scomp);

    for (int i = a_scomp+1; i < a_scomp+a_ncomp; i++)
        if (!(map == interp(i)))
            return false;

    return true;
}

std::vector< std::pair<int,int> >
StateDescriptor::sameInterps (int a_scomp,
                              int a_ncomp) const
{
    BL_ASSERT(a_scomp >= 0);
    BL_ASSERT(a_ncomp >= 1);

    std::vector< std::pair<int,int> > range;

    Interpolater* map = interp(a_scomp);

    int SComp = a_scomp, NComp = 1;

    for (int i = a_scomp+1; i < a_scomp+a_ncomp; i++)
    {
        if (map == interp(i))
        {
            NComp++;
        }
        else
        {
            range.push_back(std::pair<int,int>(SComp,NComp));

            map   = interp(i);
            SComp = i;
            NComp = 1;
        }
    }

    range.push_back(std::pair<int,int>(SComp,NComp));

#ifdef AMREX_DEBUG
    int sum = 0;
    for (int i = 0; i < static_cast<int>(range.size()); i++)
        sum += range[i].second;
    BL_ASSERT(sum == a_ncomp);
#endif

    return range;
}

}

