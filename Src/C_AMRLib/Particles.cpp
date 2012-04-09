#include <Particles.H>
#include <ParmParse.H>
#include <limits>

void
ParticleBase::AssignDensityCoeffs (const ParticleBase& p,
                                   const Real*         plo,
                                   const Real*         dx,
                                   Real*               fracs,
                                   IntVect*            cells)
{
    const Real len[BL_SPACEDIM] = { D_DECL((p.m_pos[0]-plo[0])/dx[0] + 0.5,
                                           (p.m_pos[1]-plo[1])/dx[1] + 0.5,
                                           (p.m_pos[2]-plo[2])/dx[2] + 0.5) };

    const IntVect csect(D_DECL(floor(len[0]), floor(len[1]), floor(len[2])));

    Real frac[BL_SPACEDIM] = { D_DECL(len[0]-csect[0], len[1]-csect[1], len[2]-csect[2]) };

    for (int d = 0; d < BL_SPACEDIM; d++)
    {
        if (frac[d] > 1) frac[d] = 1;
        if (frac[d] < 0) frac[d] = 0;
    }

    IntVect cell = csect;
    //
    // Note that cells[0] contains csect.
    //
#if (BL_SPACEDIM == 1)
    // High
    cells[0] = cell;
    fracs[0] = frac[0];

    // Low
    cell[0]  = cell[0] - 1;
    cells[1] = cell;
    fracs[1] = (1-frac[0]);

#elif (BL_SPACEDIM == 2)
    // HH
    cells[0] = cell;
    fracs[0] = frac[0] * frac[1] ;
    
    // LH
    cell[0]  = cell[0] - 1;
    cells[1] = cell;
    fracs[1] = (1-frac[0]) * frac[1];
    
    // LL
    cell[1]  = cell[1] - 1;
    cells[2] = cell;
    fracs[2] = (1-frac[0]) * (1-frac[1]);
    
    // HL
    cell[0]  = cell[0] + 1;
    cells[3] = cell;
    fracs[3] = frac[0] * (1-frac[1]);

#elif (BL_SPACEDIM == 3)
    // HHH
    cells[0] = cell;
    fracs[0] = frac[0] * frac[1] * frac[2];

    // LHH
    cell[0]  = cell[0] - 1;
    cells[1] = cell;
    fracs[1] = (1-frac[0]) * frac[1] * frac[2];

    // LLH
    cell[1]  = cell[1] - 1;
    cells[2] = cell;
    fracs[2] = (1-frac[0]) * (1-frac[1]) * frac[2];
    
    // HLH
    cell[0]  = cell[0] + 1;
    cells[3] = cell;
    fracs[3] = frac[0] * (1-frac[1]) * frac[2];

    cell = csect;

    // HHL
    cell[2]  = cell[2] - 1;
    cells[4] = cell;
    fracs[4] = frac[0] * frac[1] * (1-frac[2]);
    
    // LHL
    cell[0]  = cell[0] - 1;
    cells[5] = cell;
    fracs[5] = (1-frac[0]) * frac[1] * (1-frac[2]);

    // LLL
    cell[1]  = cell[1] - 1;
    cells[6] = cell;
    fracs[6] = (1-frac[0]) * (1-frac[1]) * (1-frac[2]);
    
    // HLL
    cell[0]  = cell[0] + 1;
    cells[7] = cell;
    fracs[7] = frac[0] * (1-frac[1]) * (1-frac[2]);
#endif
}

bool
ParticleBase::CrseToFine (const ParticleBase& p,
                          const Real*         plo,
                          const Real*         dx,
                          const BoxArray&     cfba,
                          const IntVect*      cells,
                          bool*               which,
                          int*                fgrid)
{
    //
    // We're in AssignDensity(). We want to know whether or not updating
    // with a particle, will we cross a  crse->fine boundary of the level
    // with coarsened fine BoxArray "cfba".  "cells" are calculated from
    // AssignDensityCoeffs().
    //
    const int M = D_TERM(2,+2,+4);

    for (int i = 0; i < M; i++)
    {
        fgrid[i] = -1;
        which[i] =  false;
    }

    bool result = false;

    std::vector< std::pair<int,Box> > isects;

    for (int i = 0; i < M; i++)
    {
        isects = cfba.intersections(Box(cells[i],cells[i]));

        if (!isects.empty())
        {
            BL_ASSERT(isects.size() == 1);

            result   = true;
            which[i] = true;
            fgrid[i] = isects[0].first;  // The grid ID at fine level that we hit.
        }
    }

    return result;
}

bool
ParticleBase::FineToCrse (const ParticleBase& p,
                          const Real*         plo,
                          const Real*         dx,
                          const Box&          vbx,
                          const BoxArray&     ba,
                          const IntVect*      cells,
                          bool*               which)
{
    //
    // We're in AssignDensity(). We want to know whether or not updating
    // with a particle at "cell", whose valid box is "vbx", will we cross a
    // fine->crse boundary of the level with BoxArray "ba".  "cells" are
    // calculated from AssignDensityCoeffs().
    //
    const int M = D_TERM(2,+2,+4);

    for (int i = 0; i < M; i++)
        which[i] = false;

    Box ibx = BoxLib::grow(vbx,-1);

    BL_ASSERT(ibx.ok());

    if (ibx.contains(p.m_cell))
        //
        // We're strictly contained in our valid box.
        // We can't cross a fine->crse boundary.
        //
        return false;
    //
    // Otherwise ...
    //
    bool result = false;

    for (int i = 0; i < M; i++)
    {
        if (!ba.contains(cells[i]))
        {
            result   = true;
            which[i] = true;
        }
    }

    return result;
}

int
ParticleBase::MaxReaders ()
{
    const int Max_Readers_def = 64;

    static int Max_Readers;

    static bool first = true;

    if (first)
    {
        first = false;

        ParmParse pp("particles");

        Max_Readers = Max_Readers_def;

        pp.query("nreaders", Max_Readers);

        Max_Readers = std::min(ParallelDescriptor::NProcs(),Max_Readers);

        if (Max_Readers <= 0)
            BoxLib::Abort("particles.nreaders must be positive");
    }

    return Max_Readers;
}

const std::string&
ParticleBase::DataPrefix ()
{
    //
    // The actual particle data is stored in files of the form: DATA_nnnn.
    //
    static const std::string data("DATA_");

    return data;
}

const std::string&
ParticleBase::Version ()
{
    //
    // If we change the Checkpoint/Restart format we should increment this.
    //
    static const std::string version("Version_One_Dot_Zero");

    return version;
}

static int the_next_id = 1;

int
ParticleBase::NextID ()
{
    int next;

#ifdef _OPENMP
#pragma omp critical(nextid_lock)
#endif
    {
        if (the_next_id == std::numeric_limits<int>::max())
            BoxLib::Abort("ParticleBase::NextID() -- too many particles");

        next = the_next_id++;
    }

    return next;
}

void
ParticleBase::NextID (int nextid)
{
    the_next_id = nextid;
}

IntVect
ParticleBase::Index (const ParticleBase& p,
                     int                 lev,
                     const Amr*          amr)
{
    BL_ASSERT(amr != 0);
    BL_ASSERT(lev >= 0 && lev <= amr->finestLevel());

    IntVect iv;

    const Geometry& geom = amr->Geom(lev);

    D_TERM(iv[0]=floor((p.m_pos[0]-geom.ProbLo(0))/geom.CellSize(0));,
           iv[1]=floor((p.m_pos[1]-geom.ProbLo(1))/geom.CellSize(1));,
           iv[2]=floor((p.m_pos[2]-geom.ProbLo(2))/geom.CellSize(2)););

    iv += geom.Domain().smallEnd();

    return iv;
}

bool
ParticleBase::Where (ParticleBase& p,
                     const Amr*    amr,
                     bool          update)
{
    BL_ASSERT(amr != 0);

    if (update)
    {
        //
        // We have a valid particle whose position has changed slightly.
        // Try to update m_cell and m_grid smartly.
        //
        BL_ASSERT(p.m_id > 0);
        BL_ASSERT(p.m_grid >= 0 && p.m_grid < amr->boxArray(p.m_lev).size());

        IntVect iv = Index(p,p.m_lev,amr);

        if (p.m_cell == iv)
            //
            // The particle hasn't left its cell.
            //
            return true;

        if (p.m_lev == amr->finestLevel())
        {
            p.m_cell = iv;

            if (amr->boxArray(p.m_lev)[p.m_grid].contains(p.m_cell))
                //
                // It has left its cell but is still in the same grid.
                //
                return true;
        }
    }

    std::vector< std::pair<int,Box> > isects;

    for (int lev = amr->finestLevel(); lev >= 0; lev--)
    {
        IntVect iv = Index(p,lev,amr);

        isects = amr->boxArray(lev).intersections(Box(iv,iv));

        if (!isects.empty())
        {
            BL_ASSERT(isects.size() == 1);

            p.m_lev  = lev;
            p.m_grid = isects[0].first;
            p.m_cell = iv;

            return true;
        }
    }

    return false;
}

void
ParticleBase::PeriodicShift (ParticleBase& p,
                             const Amr*    amr)
{
    //
    // This routine should only be called when Where() returns false.
    //
    BL_ASSERT(amr != 0);
    //
    // We'll use level 0 stuff since ProbLo/ProbHi are the same for every level.
    //
    const Geometry& geom = amr->Geom(0);
    const Box&      dmn  = geom.Domain();
    IntVect         iv   = Index(p,0,amr);
    const Real      eps  = 1.e-13;

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        if (!geom.isPeriodic(i)) continue;

        if (iv[i] > dmn.bigEnd(i))
        {
            if (p.m_pos[i] == geom.ProbHi(i))
                //
                // Don't let particles lie exactly on the domain face.
                // Force the particle to be outside the domain so the
                // periodic shift will bring it back inside.
                //
                p.m_pos[i] += eps;

            p.m_pos[i] -= geom.ProbLength(i);

            BL_ASSERT(p.m_pos[i] >= geom.ProbLo(i));
        }
        else if (iv[i] < dmn.smallEnd(i))
        {
            if (p.m_pos[i] == geom.ProbLo(i))
                //
                // Don't let particles lie exactly on the domain face.
                // Force the particle to be outside the domain so the
                // periodic shift will bring it back inside.
                //
                p.m_pos[i] -= eps;

            p.m_pos[i] += geom.ProbLength(i);

            BL_ASSERT(p.m_pos[i] <= geom.ProbHi(i));
        }
    }
    //
    // The particle may still be outside the domain in the case
    // where we aren't periodic on the face out which it travelled.
    //
}

void
ParticleBase::Reset (ParticleBase& p,
                     const Amr*    amr,
                     bool          update)
{
    BL_ASSERT(amr != 0);

    if (!ParticleBase::Where(p,amr,update))
    {
        //
        // Here's where we need to deal with boundary conditions.
        //
        // Attempt to shift the particle back into the domain if it
        // crossed a periodic boundary.  Otherwise (for now) we
        // invalidate the particle.
        //
        ParticleBase::PeriodicShift(p,amr);

        if (!ParticleBase::Where(p,amr))
        {
#ifdef _OPENMP
#pragma omp critical(reset_lock)
#endif
            {
                std::cout << "Invalidating out-of-domain particle: " << p << '\n';
            }

            BL_ASSERT(p.m_id > 0);

            p.m_id = -p.m_id;
        }
    }
}

Real
ParticleBase::InterpDoit (const FArrayBox& fab,
                          const IntVect&   hi,
                          const Real*      frac,
                          int              comp)
{
    IntVect cell = hi;

    Real val = 0;

    const Real ifrac[BL_SPACEDIM] = { D_DECL((1-frac[0]), (1-frac[1]), (1-frac[2])) };

#if (BL_SPACEDIM == 1)
    // High
    val += fab(cell,comp) * frac[0];

    // Low
    cell[0] = cell[0] - 1;
    val += fab(cell,comp) * ifrac[0];

#elif (BL_SPACEDIM == 2)
    // HH
    val += fab(cell,comp) *    frac[0]  *    frac[1] ;

    // LH
    cell[0] = cell[0] - 1;
    val += fab(cell,comp) * ifrac[0] *    frac[1] ;

    // LL
    cell[1]   = cell[1] - 1;
    val += fab(cell,comp) * ifrac[0] * ifrac[1];

    // HL
    cell[0] = cell[0] + 1;
    val += fab(cell,comp) *    frac[0]  * ifrac[1];

#elif (BL_SPACEDIM == 3)
    // HHH
    val += fab(cell,comp) *    frac[0]  *    frac[1]  *    frac[2] ;
   
    // LHH
    cell[0] = cell[0] - 1;
    val += fab(cell,comp) * ifrac[0] *    frac[1]  *    frac[2] ;
   
    // LLH
    cell[1] = cell[1] - 1;
    val += fab(cell,comp) * ifrac[0] * ifrac[1] *    frac[2] ;
   
    // HLH
    cell[0] = cell[0] + 1;
    val += fab(cell,comp) *    frac[0]  * ifrac[1] *    frac[2] ;

    cell    = hi;
    cell[2] = cell[2] - 1;

    // HHL
    val += fab(cell,comp) *    frac[0]  *    frac[1]  *    ifrac[2] ;
   
    // LHL
    cell[0] = cell[0] - 1;
    val += fab(cell,comp) * ifrac[0] *    frac[1]  *    ifrac[2] ;
   
    // LLL
    cell[1] = cell[1] - 1;
    val += fab(cell,comp) * ifrac[0] * ifrac[1] *    ifrac[2] ;
   
    // HLL
    cell[0] = cell[0] + 1;
    val += fab(cell,comp) *    frac[0]  * ifrac[1] *    ifrac[2] ;

#endif

    return val;
}

void
ParticleBase::Interp (const ParticleBase& prt,
                      const Amr*          amr,
                      const FArrayBox&    fab,
                      const int*          idx,
                      Real*               val,
                      int                 cnt)
{
    BL_ASSERT(amr != 0);
    BL_ASSERT(idx != 0);
    BL_ASSERT(val != 0);

    const Geometry& gm = amr->Geom(prt.m_lev);
    const Real*     dx = gm.CellSize();

    const Real len[BL_SPACEDIM] = { D_DECL((prt.m_pos[0]-gm.ProbLo(0))/dx[0] + 0.5,
                                           (prt.m_pos[1]-gm.ProbLo(1))/dx[1] + 0.5,
                                           (prt.m_pos[2]-gm.ProbLo(2))/dx[2] + 0.5) };

    const IntVect hi(D_DECL(floor(len[0]), floor(len[1]), floor(len[2])));

    Real frac[BL_SPACEDIM] = { D_DECL(len[0]-hi[0], len[1]-hi[1], len[2]-hi[2]) };

    for (int d = 0; d < BL_SPACEDIM; d++)
    {
        if (frac[d] > 1) frac[d] = 1;
        if (frac[d] < 0) frac[d] = 0;
    }

    for (int i = 0; i < cnt; i++)
    {
        const int comp = idx[i];

        BL_ASSERT(comp >= 0 && comp < fab.nComp());

        val[i] = ParticleBase::InterpDoit(fab,hi,frac,comp);
    }
}

void
ParticleBase::GetGravity (const FArrayBox&    gfab,
                          const Amr*          amr,
                          const ParticleBase& p,
                          Real*               grav)
{
    BL_ASSERT(amr  != 0);
    BL_ASSERT(grav != 0);

    int idx[BL_SPACEDIM] = { D_DECL(0,1,2) };

    ParticleBase::Interp(p,amr,gfab,idx,grav,BL_SPACEDIM);
}

std::ostream&
operator<< (std::ostream& os, const ParticleBase& p)
{
    os << p.m_id   << ' '
       << p.m_cpu  << ' '
       << p.m_lev  << ' '
       << p.m_grid << ' '
       << p.m_cell << ' ';

    for (int i = 0; i < BL_SPACEDIM; i++)
        os << p.m_pos[i] << ' ';

    if (!os.good())
        BoxLib::Error("operator<<(ostream&,ParticleBase&) failed");

    return os;
}
