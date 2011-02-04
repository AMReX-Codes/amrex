#include "Particles.H"

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

#ifdef BL_USE_OMP
#pragma omp critical(nextid_lock)
#endif
    {
        next = the_next_id++;
    }

    return next;
}

void
ParticleBase::NextID (int nextid)
{
    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Restoring ParticleBase::NextID: " << nextid << std::endl;
    }
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
    const Real      eps  = 1.e13;

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
ParticleBase::GetGravity (const FArrayBox&    gfab,
                          const Amr*          amr,
                          int                 lev,
                          const ParticleBase& p,
                          Real*               grav)
{
    BL_ASSERT(amr  != 0);
    BL_ASSERT(grav != 0);

    const Geometry& geom = amr->Geom(lev);
    const Real*     dx   = geom.CellSize();

    const IntVect csect(D_DECL(floor((p.m_pos[0]-geom.ProbLo(0))/dx[0] + 0.5),
                               floor((p.m_pos[1]-geom.ProbLo(1))/dx[1] + 0.5),
                               floor((p.m_pos[2]-geom.ProbLo(2))/dx[2] + 0.5)));

    const Real frac[BL_SPACEDIM] = { D_DECL(-csect[0] + p.m_pos[0]/dx[0] + 0.5,
                                            -csect[1] + p.m_pos[1]/dx[1] + 0.5,
                                            -csect[2] + p.m_pos[2]/dx[2] + 0.5) };

    for (int dm = 0; dm < BL_SPACEDIM; dm++)
    {
        IntVect  cell = csect;
        grav[dm] = 0;

#if (BL_SPACEDIM == 1)
        // High
        grav[dm] += gfab(cell, dm) * frac[0];

        // Low
        cell[0]   = cell[0] - 1;
        grav[dm] += gfab(cell, dm) * (1-frac[0]);
#endif
            
#if (BL_SPACEDIM == 2)
        // HH
        grav[dm] += gfab(cell, dm) *    frac[0]  *    frac[1] ;

        // LH
        cell[0]   = cell[0] - 1;
        grav[dm] += gfab(cell, dm) * (1-frac[0]) *    frac[1] ;

        // LL
        cell[1]   = cell[1] - 1;
        grav[dm] += gfab(cell, dm) * (1-frac[0]) * (1-frac[1]);

        // HL
        cell[0]   = cell[0] + 1;
        grav[dm] += gfab(cell, dm) *    frac[0]  * (1-frac[1]);
#endif
 

#if (BL_SPACEDIM == 3)

        // HHH
        grav[dm] += gfab(cell, dm) *    frac[0]  *    frac[1]  *    frac[2] ;
   
        // LHH
        cell[0]   = cell[0] - 1;
        grav[dm] += gfab(cell, dm) * (1-frac[0]) *    frac[1]  *    frac[2] ;
   
        // LLH
        cell[1]   = cell[1] - 1;
        grav[dm] += gfab(cell, dm) * (1-frac[0]) * (1-frac[1]) *    frac[2] ;
   
        // HLH
        cell[0]   = cell[0] + 1;
        grav[dm] += gfab(cell, dm) *    frac[0]  * (1-frac[1]) *    frac[2] ;

        cell     = csect;
        cell[2]  = cell[2] - 1;

        // HHL
        grav[dm] += gfab(cell, dm) *    frac[0]  *    frac[1]  *    (1-frac[2]) ;
   
        // LHL
        cell[0]   = cell[0] - 1;
        grav[dm] += gfab(cell, dm) * (1-frac[0]) *    frac[1]  *    (1-frac[2]) ;
   
        // LLL
        cell[1]   = cell[1] - 1;
        grav[dm] += gfab(cell, dm) * (1-frac[0]) * (1-frac[1]) *    (1-frac[2]) ;
   
        // HLL
        cell[0]   = cell[0] + 1;
        grav[dm] += gfab(cell, dm) *    frac[0]  * (1-frac[1]) *    (1-frac[2]) ;
#endif
    }
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
