#include "Particles.H"

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
    const Geometry& geom    = amr->Geom(0);
    const Box&      dmn     = geom.Domain();
    IntVect         iv      = Index(p,0,amr);
    bool            shifted = false;

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        if (iv[i] > dmn.bigEnd(i))
        {
            shifted = true;

            p.m_pos[i] -= geom.ProbLength(i);

            BL_ASSERT(p.m_pos[i] >= geom.ProbLo(i));
        }
        else if (iv[i] < dmn.smallEnd(i))
        {
            shifted = true;

            p.m_pos[i] += geom.ProbLength(i);

            BL_ASSERT(p.m_pos[i] <= geom.ProbHi(i));
        }
    }
    BL_ASSERT(shifted);
    //
    // Modulo any arithmetic funny business the particle should
    // now be back in the problem domain.
    //
    if (!ParticleBase::Where(p,amr))
    {
        std::cout.precision(16);

        std::cout << "ParticleBase::PeriodicShift(): out-of-domain particle"
                  << p
                  << "\ngeom: "
                  << geom << std::endl;
        BoxLib::Abort();
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

    os << '\n';

    return os;
}
