#include <Particles.H>
#include <ParmParse.H>
#include <limits>

void
ParticleBase::CIC_Cells_Fracs_Basic (const ParticleBase& p,
                                     const Real*         plo,
                                     const Real*         dx,
                                     Real*               fracs,
                                     IntVect*            cells)
{
    BL_PROFILE("ParticleBase::CIC_Cells_Fracs_B()");
    //
    // "fracs" should be dimensioned: Real    fracs[D_TERM(2,+2,+4)]
    //
    // "cells" should be dimensioned: IntVect cells[D_TERM(2,+2,+4)]
    //
    const Real len[BL_SPACEDIM] = { D_DECL((p.m_pos[0]-plo[0])/dx[0] + Real(0.5),
                                           (p.m_pos[1]-plo[1])/dx[1] + Real(0.5),
                                           (p.m_pos[2]-plo[2])/dx[2] + Real(0.5)) };

    const IntVect cell(D_DECL(floor(len[0]), floor(len[1]), floor(len[2])));

    const Real frac[BL_SPACEDIM] = { D_DECL(len[0]-cell[0], len[1]-cell[1], len[2]-cell[2]) };

    ParticleBase::CIC_Fracs(frac, fracs);
    ParticleBase::CIC_Cells(cell, cells);
}

int
ParticleBase::CIC_Cells_Fracs (const ParticleBase& p,
                               const Real*         plo,
                               const Real*         dx_geom,
                               const Real*         dx_part,
                               Array<Real>&        fracs,
                               Array<IntVect>&     cells)
{
    BL_PROFILE("ParticleBase::CIC_Cells_Fracs()");
    if (dx_geom == dx_part)
    {
        const int M = D_TERM(2,+2,+4);
        fracs.resize(M);
        cells.resize(M);
        ParticleBase::CIC_Cells_Fracs_Basic(p,plo,dx_geom,fracs.dataPtr(),cells.dataPtr());
        return M;
    }
    //
    // The first element in fracs and cells is the lowest corner, the last is the highest.
    //
    const Real hilen[BL_SPACEDIM] = { D_DECL((p.m_pos[0]-plo[0]+dx_part[0]/2)/dx_geom[0],
                                             (p.m_pos[1]-plo[1]+dx_part[1]/2)/dx_geom[1],
                                             (p.m_pos[2]-plo[2]+dx_part[2]/2)/dx_geom[2]) };

    const Real lolen[BL_SPACEDIM] = { D_DECL((p.m_pos[0]-plo[0]-dx_part[0]/2)/dx_geom[0],
                                             (p.m_pos[1]-plo[1]-dx_part[1]/2)/dx_geom[1],
                                             (p.m_pos[2]-plo[2]-dx_part[2]/2)/dx_geom[2]) };

    const IntVect hicell(D_DECL(floor(hilen[0]), floor(hilen[1]), floor(hilen[2])));
    
    const IntVect locell(D_DECL(floor(lolen[0]), floor(lolen[1]), floor(lolen[2])));
    
    const Real cell_density = D_TERM(dx_geom[0]/dx_part[0],*dx_geom[1]/dx_part[1],*dx_geom[2]/dx_part[2]);
    
    const int M = D_TERM((hicell[0]-locell[0]+1),*(hicell[1]-locell[1]+1),*(hicell[2]-locell[2]+1));

    fracs.resize(M);
    cells.resize(M);
    //
    // This portion might be slightly inefficient. Feel free to redo it if need be.
    //
    int i = 0;
#if (BL_SPACEDIM == 1)
    for (int xi = locell[0]; xi <= hicell[0]; xi++)
    {
        cells[i][0] = xi;
        fracs[i] = (std::min(hilen[0]-xi,Real(1))-std::max(lolen[0]-xi,Real(0)))*cell_density;
        i++;
    }
#elif (BL_SPACEDIM == 2)
    for (int yi = locell[1]; yi <= hicell[1]; yi++)
    {
        const Real yf = std::min(hilen[1]-yi,Real(1))-std::max(lolen[1]-yi,Real(0));
        for (int xi = locell[0]; xi <= hicell[0]; xi ++)
        {
            cells[i][0] = xi;
            cells[i][1] = yi;
            fracs[i] = yf * (std::min(hilen[0]-xi,Real(1))-std::max(lolen[0]-xi,Real(0)))*cell_density;
            i++;
        }
    }
#elif (BL_SPACEDIM == 3)
    for (int zi = locell[2]; zi <= hicell[2]; zi++)
    {
        const Real zf = std::min(hilen[2]-zi,Real(1))-std::max(lolen[2]-zi,Real(0));
        for (int yi = locell[1]; yi <= hicell[1]; yi++)
        {
            const Real yf = std::min(hilen[1]-yi,Real(1))-std::max(lolen[1]-yi,Real(0));
            for (int xi = locell[0]; xi <= hicell[0]; xi++)
            {
                cells[i][0] = xi;
                cells[i][1] = yi;
                cells[i][2] = zi;
                fracs[i] = zf * yf * (std::min(hilen[0]-xi,Real(1))-std::max(lolen[0]-xi,Real(0))) * cell_density;
                i++;
            }
        }
    }
#endif

    return M;
}

bool
ParticleBase::CrseToFine (const BoxArray&       cfba,
                          const Array<IntVect>& cells,
                          Array<IntVect>&       cfshifts,
                          const Geometry&       gm,
                          Array<int>&           which,
                          Array<IntVect>&       pshifts)
{
    BL_PROFILE("ParticleBase::CrseToFine()");
    //
    // We're in AssignDensity(). We want to know whether or not updating
    // with a particle, will we cross a  crse->fine boundary of the level
    // with coarsened fine BoxArray "cfba".  "cells" are as calculated from
    // CIC_Cells_Fracs().
    //
    const int M = cells.size();

    which.resize(M);
    cfshifts.resize(M);

    for (int i = 0; i < M; i++)
        which[i] =  0;

    bool result = false;

    for (int i = 0; i < M; i++)
    {
        if (cfba.contains(cells[i]))
        {
            result      = true;
            which[i]    = 1;
            cfshifts[i] = IntVect::TheZeroVector();
        }
        else if (!gm.Domain().contains(cells[i]))
        {
            BL_ASSERT(gm.isAnyPeriodic());
            //
            // Can the cell be shifted into cfba?
            //
            const Box bx(cells[i],cells[i]);

            gm.periodicShift(bx, gm.Domain(), pshifts);

            if (!pshifts.empty())
            {
                BL_ASSERT(pshifts.size() == 1);

                const Box& dbx = bx - pshifts[0];

                BL_ASSERT(dbx.ok());

                if (cfba.contains(dbx))
                {
                    //
                    // Note that pshifts[0] is from the coarse perspective.
                    // We'll later need to multiply it by ref ratio to use
                    // at the fine level.
                    //
                    result      = true;
                    which[i]    = 1;
                    cfshifts[i] = pshifts[0];
                }
            }
        }
    }

    return result;
}

bool
ParticleBase::FineToCrse (const ParticleBase&                p,
                          int                                flev,
                          const ParGDBBase*                  gdb,
                          const Array<IntVect>&              fcells,
                          const BoxArray&                    fvalid,
                          const BoxArray&                    compfvalid_grown,
                          Array<IntVect>&                    ccells,
                          Array<Real>&                       cfracs,
                          Array<int>&                        which,
                          Array<int>&                        cgrid,
                          Array<IntVect>&                    pshifts,
                          std::vector< std::pair<int,Box> >& isects)
{
    BL_PROFILE("ParticleBase::FineToCrse()");
    BL_ASSERT(gdb != 0);
    BL_ASSERT(flev > 0);
    //
    // We're in AssignDensity(). We want to know whether or not updating
    // with a particle we'll cross a fine->crse boundary.  Note that crossing
    // a periodic boundary, where the periodic shift lies in our valid region,
    // is not considered a Fine->Crse crossing.
    //
    const int M = fcells.size();

    which.resize(M);
    cgrid.resize(M);
    ccells.resize(M);
    cfracs.resize(M);

    for (int i = 0; i < M; i++)
    {
        cgrid[i] = -1;
        which[i] =  0;
    }

    const Box& ibx = BoxLib::grow(gdb->ParticleBoxArray(flev)[p.m_grid],-1);

    BL_ASSERT(ibx.ok());

    if (ibx.contains(p.m_cell))
        //
        // We're strictly contained in our valid box.
        // We can't cross a fine->crse boundary.
        //
        return false;

    if (!compfvalid_grown.contains(p.m_cell))
        //
        // We're strictly contained in our "valid" region. Note that the valid
        // region contains any periodically shifted ghost cells that intersect
        // valid region.
        //
        return false;
    //
    // Otherwise ...
    //
    const Geometry& cgm = gdb->Geom(flev-1);
    const IntVect&  rr  = gdb->refRatio(flev-1);
    const BoxArray& cba = gdb->ParticleBoxArray(flev-1);

    ParticleBase::CIC_Cells_Fracs(p, cgm.ProbLo(), cgm.CellSize(), cfracs, ccells);

    bool result = false;

    for (int i = 0; i < M; i++)
    {
        IntVect ccell_refined = ccells[i]*rr;
        //
        // We've got to protect against the case when we're at the low
        // end of the domain because coarsening & refining don't work right
        // when indices go negative.
        //
        for (int dm = 0; dm < BL_SPACEDIM; dm++)
            ccell_refined[dm] = std::max(ccell_refined[dm], -1);

        if (!fvalid.contains(ccell_refined))
        {
            result   = true;
            which[i] = 1;

            Box cbx(ccells[i],ccells[i]);
    
            if (!cgm.Domain().contains(ccells[i]))
            {
                //
                // We must be at a periodic boundary.
                // Find valid box into which we can be periodically shifted.
                //
                BL_ASSERT(cgm.isAnyPeriodic());

                cgm.periodicShift(cbx, cgm.Domain(), pshifts);

                BL_ASSERT(pshifts.size() == 1);

                cbx -= pshifts[0];

                ccells[i] -= pshifts[0];
                BL_ASSERT(cbx.ok());
                BL_ASSERT(cgm.Domain().contains(cbx));
            }
            //
            // Which grid at the crse level do we need to update?
            //
            cba.intersections(cbx,isects,true,0);

            BL_ASSERT(!isects.empty());

            cgrid[i] = isects[0].first;  // The grid ID at crse level that we hit.
        }
    }

    return result;
}

void
ParticleBase::FineCellsToUpdateFromCrse (const ParticleBase&                p,
                                         int                                lev,
                                         const ParGDBBase*                  gdb,
                                         const IntVect&                     ccell,
                                         const IntVect&                     cshift,
                                         Array<int>&                        fgrid,
                                         Array<Real>&                       ffrac,
                                         Array<IntVect>&                    fcells,
                                         std::vector< std::pair<int,Box> >& isects)
{
    BL_PROFILE("ParticleBase::FineCellsToUpdateFromCrse()");
    BL_ASSERT(lev >= 0);
    BL_ASSERT(lev < gdb->finestLevel());

    const Box&      fbx = BoxLib::refine(Box(ccell,ccell),gdb->refRatio(lev));
    const BoxArray& fba = gdb->ParticleBoxArray(lev+1);
    const Real*     plo = gdb->Geom(lev).ProbLo();
    const Real*     dx  = gdb->Geom(lev).CellSize();
    const Real*     fdx = gdb->Geom(lev+1).CellSize();

    if (cshift == IntVect::TheZeroVector())
    {
        BL_ASSERT(fba.contains(fbx));
    }
    //
    // Instead of clear()ing these we'll do a resize(0).
    // This'll preserve their capacity so that we'll only need
    // to do any memory allocation when their capacity needs to increase.
    //
    fgrid.resize(0);
    ffrac.resize(0);
    fcells.resize(0);
    //
    // Which fine cells does particle "p" (that wants to update "ccell") do we
    // touch at the finer level?
    //
    for (IntVect iv = fbx.smallEnd(); iv <= fbx.bigEnd(); fbx.next(iv))
    {
        bool touches = true;

        for (int k = 0; k < BL_SPACEDIM; k++)
        {
            const Real celllo = iv[k]  * fdx[k] + plo[k];
            const Real cellhi = celllo + fdx[k];

            if ((p.m_pos[k] < celllo) && (celllo > (p.m_pos[k] + dx[k]/2)))
                touches = false;

            if ((p.m_pos[k] > cellhi) && (cellhi < (p.m_pos[k] - dx[k]/2)))
                touches = false;
        }

        if (touches)
        {
            fcells.push_back(iv);
        }
    }

    Real sum_fine = 0;
    //
    // We need to figure out the fine fractions and the fine grid needed updating.
    //
    for (int j = 0; j < fcells.size(); j++)
    {
        IntVect& iv = fcells[j];

        Real the_frac = 1;

        for (int k = 0; k < BL_SPACEDIM; k++)
        {
            const Real celllo = (iv[k] * fdx[k] + plo[k]);

            if (p.m_pos[k] <= celllo)
            {
                const Real isecthi = p.m_pos[k] + dx[k]/2;

                the_frac *= std::min((isecthi - celllo),fdx[k]);
            }
            else
            {
                const Real cellhi  = (iv[k]+1) * fdx[k] + plo[k];
                const Real isectlo = p.m_pos[k] - dx[k]/2;

                the_frac *= std::min((cellhi - isectlo),fdx[k]);
            }
        }

        ffrac.push_back(the_frac);

        sum_fine += the_frac;

        if (cshift != IntVect::TheZeroVector())
        {
            //
            // Update to the correct fine cell needing updating.
            // Note that "cshift" is from the coarse perspective.
            //
            const IntVect& fshift = cshift * gdb->refRatio(lev);
            //
            // Update fcells[j] to indicate a shifted fine cell needing updating.
            //
            iv -= fshift;
        }

        fba.intersections(Box(iv,iv),isects,true,0);

        BL_ASSERT(!isects.empty());

        fgrid.push_back(isects[0].first);
    }

    BL_ASSERT(ffrac.size() == fcells.size());
    BL_ASSERT(fgrid.size() == fcells.size());
    //
    // Now adjust the fine fractions so they sum to one.
    //
    for (int j = 0; j < ffrac.size(); j++)
        ffrac[j] /= sum_fine;
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

long
ParticleBase::MaxParticlesPerRead ()
{
    //
    // This is the maximum particles that "each" reader will attempt to read
    // before doing a Redistribute(). 
    //
    const long Max_Particles_Per_Read_def = 100000;

    static long Max_Particles_Per_Read;

    static bool first = true;

    if (first)
    {
        first = false;

        ParmParse pp("particles");

        Max_Particles_Per_Read = Max_Particles_Per_Read_def;

        pp.query("nparts_per_read", Max_Particles_Per_Read);

        if (Max_Particles_Per_Read <= 0)
            BoxLib::Abort("particles.nparts_per_read must be positive");
    }

    return Max_Particles_Per_Read;
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
    // Previous version strings:
    //
    //    "Version_One_Dot_Zero"
    //
    static const std::string version("Version_One_Dot_One");

    return version;
}

static int the_next_id = 1;

int
ParticleBase::NextID ()
{
    int next;

#ifdef _OPENMP
#pragma omp atomic capture
#endif
    next = the_next_id++;

    if (next == std::numeric_limits<int>::max())
	BoxLib::Abort("ParticleBase::NextID() -- too many particles");

    return next;
}

int
ParticleBase::UnprotectedNextID ()
{
    int next = the_next_id++;
    if (next == std::numeric_limits<int>::max())
	BoxLib::Abort("ParticleBase::NextID() -- too many particles");
    return next;
}

void
ParticleBase::NextID (int nextid)
{
    the_next_id = nextid;
}

IntVect
ParticleBase::Index (const ParticleBase& p,
                     const Geometry&     geom)
{
    IntVect iv;

    D_TERM(iv[0]=floor((p.m_pos[0]-geom.ProbLo(0))/geom.CellSize(0));,
           iv[1]=floor((p.m_pos[1]-geom.ProbLo(1))/geom.CellSize(1));,
           iv[2]=floor((p.m_pos[2]-geom.ProbLo(2))/geom.CellSize(2)););

    iv += geom.Domain().smallEnd();

    return iv;
}

bool
ParticleBase::Where (ParticleBase& p,
                     const ParGDBBase* gdb,
                     int           lev_min,
                     int           finest_level)
{
    BL_ASSERT(gdb != 0);

    if (finest_level == -1)
        finest_level = gdb->finestLevel();

    BL_ASSERT(finest_level <= gdb->finestLevel());

    std::vector< std::pair<int,Box> > isects;

    for (int lev = finest_level; lev >= lev_min; lev--)
    {
        const IntVect& iv = ParticleBase::Index(p,gdb->Geom(lev));

	if (lev == p.m_lev) { 
            // We may take a shortcut because the fact that we are here means 
            // this particle does not belong to any finer grids.
	    const BoxArray& ba = gdb->ParticleBoxArray(p.m_lev);
	    if (0 <= p.m_grid && p.m_grid < ba.size() && ba[p.m_grid].contains(iv)) {
		p.m_cell = iv;
		return true;
	    }
	}

        gdb->ParticleBoxArray(lev).intersections(Box(iv,iv),isects,true,0);

        if (!isects.empty())
        {
            p.m_lev  = lev;
            p.m_grid = isects[0].first;
            p.m_cell = iv;

            return true;
        }
    }
    return false;
}

bool
ParticleBase::PeriodicWhere (ParticleBase& p,
                             const ParGDBBase*    gdb,
                             int           lev_min,
                             int           finest_level)
{
    BL_ASSERT(gdb != 0);

    if (!gdb->Geom(0).isAnyPeriodic()) return false;

    if (finest_level == -1)
        finest_level = gdb->finestLevel();

    BL_ASSERT(finest_level <= gdb->finestLevel());
    //
    // Create a copy "dummy" particle to check for periodic outs.
    //
    ParticleBase p_prime = p;

    if (ParticleBase::PeriodicShift(p_prime, gdb))
    {
        std::vector< std::pair<int,Box> > isects;

        for (int lev = finest_level; lev >= lev_min; lev--)
        {
            const IntVect& iv = ParticleBase::Index(p_prime,gdb->Geom(lev));

            gdb->ParticleBoxArray(lev).intersections(Box(iv,iv),isects,true,0);

            if (!isects.empty())
            {
                D_TERM(p.m_pos[0] = p_prime.m_pos[0];,
                       p.m_pos[1] = p_prime.m_pos[1];,
                       p.m_pos[2] = p_prime.m_pos[2];);

                p.m_lev  = lev;
                p.m_grid = isects[0].first;
                p.m_cell = iv;

                return true;
            }
        }
    }

    return false;
}

bool
ParticleBase::RestrictedWhere (ParticleBase& p,
                               const ParGDBBase*    gdb,
                               int           ngrow)
{
    BL_ASSERT(gdb != 0);

    const IntVect& iv = ParticleBase::Index(p,gdb->Geom(p.m_lev));

    if (BoxLib::grow(gdb->ParticleBoxArray(p.m_lev)[p.m_grid], ngrow).contains(iv))
    {
        p.m_cell = iv;

        return true;
    }

    return false;
}

bool 
ParticleBase::SingleLevelWhere (ParticleBase& p, 
                                const ParGDBBase*    gdb,
                                int           level)
{
    BL_ASSERT(gdb != 0);

    const IntVect& iv = ParticleBase::Index(p,gdb->Geom(level));

    std::vector< std::pair<int,Box> > isects;

    gdb->ParticleBoxArray(level).intersections(Box(iv,iv),isects,true,0);

    if (!isects.empty())
    {
        p.m_lev  = level;
        p.m_grid = isects[0].first;
        p.m_cell = iv;

        return true;
    }

    return false;
}

bool
ParticleBase::PeriodicShift (ParticleBase& p,
                             const ParGDBBase*    gdb)
{
    BL_PROFILE("ParticleBase::PeriodicShift()");
    //
    // This routine should only be called when Where() returns false.
    //
    BL_ASSERT(gdb != 0);
    //
    // We'll use level 0 stuff since ProbLo/ProbHi are the same for every level.
    //
    const Geometry& geom    = gdb->Geom(0);
    const Box&      dmn     = geom.Domain();
    const IntVect&  iv      = ParticleBase::Index(p,gdb->Geom(0));
    bool            shifted = false;  

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
                p.m_pos[i] += .125*geom.CellSize(i);

            p.m_pos[i] -= geom.ProbLength(i);

            if (p.m_pos[i] <= geom.ProbLo(i))
                //
                // This can happen due to precision issues.
                //
                p.m_pos[i] += .125*geom.CellSize(i);

            BL_ASSERT(p.m_pos[i] >= geom.ProbLo(i));

            shifted = true;
        }
        else if (iv[i] < dmn.smallEnd(i))
        {
            if (p.m_pos[i] == geom.ProbLo(i))
                //
                // Don't let particles lie exactly on the domain face.
                // Force the particle to be outside the domain so the
                // periodic shift will bring it back inside.
                //
                p.m_pos[i] -= .125*geom.CellSize(i);

            p.m_pos[i] += geom.ProbLength(i);

            if (p.m_pos[i] >= geom.ProbHi(i))
                //
                // This can happen due to precision issues.
                //
                p.m_pos[i] -= .125*geom.CellSize(i);

            BL_ASSERT(p.m_pos[i] <= geom.ProbHi(i));

            shifted = true;
        }
    }
    //
    // The particle may still be outside the domain in the case
    // where we aren't periodic on the face out which it travelled.
    //
    return shifted;
}

void
ParticleBase::Reset (ParticleBase& p,
                     const ParGDBBase*    gdb,
                     bool          update,
		     bool          verbose)
{
    BL_PROFILE("ParticleBase::Reset()");
    BL_ASSERT(gdb != 0);

    bool ok = ParticleBase::Where(p,gdb);

    if (!ok && gdb->Geom(0).isAnyPeriodic())
    {
        // Attempt to shift the particle back into the domain if it
        // crossed a periodic boundary.
	ParticleBase::PeriodicShift(p,gdb);
	ok = ParticleBase::Where(p,gdb);
    }
    
    if (!ok) {
        // invalidate the particle.
	if (verbose) {
#ifdef _OPENMP
#pragma omp critical(reset_lock)
#endif
	    {
		std::cout << "Invalidating out-of-domain particle: " << p << '\n';
	    }
	}

	BL_ASSERT(p.m_id > 0);

	p.m_id = -p.m_id;
    }
}

Real
ParticleBase::InterpDoit (const FArrayBox& fab,
                          const IntVect&   cell,
                          const Real*      frac,
                          int              comp)
{
    BL_PROFILE("ParticleBase::InterpDoit(fcfc)");
    const int M = D_TERM(2,+2,+4);

    Real    fracs[M];
    IntVect cells[M];

    ParticleBase::CIC_Fracs(frac, fracs);
    ParticleBase::CIC_Cells(cell, cells);

    Real val = ParticleBase::InterpDoit(fab,fracs,cells,comp);

    return val;
}

Real
ParticleBase::InterpDoit (const FArrayBox& fab,
                          const Real*      fracs,
                          const IntVect*   cells,
                          int              comp)
{
    BL_PROFILE("ParticleBase::InterpDoit(ffcc)");
    const int M = D_TERM(2,+2,+4);

    Real val = 0;

    for (int i = 0; i < M; i++)
    {
        val += fab(cells[i],comp) * fracs[i];
    }

    return val;
}

void
ParticleBase::Interp (const ParticleBase& prt,
                      const Geometry&     geom,
                      const FArrayBox&    fab,
                      const int*          idx,
                      Real*               val,
                      int                 cnt)
{
    BL_PROFILE("ParticleBase::Interp()");
    BL_ASSERT(idx != 0);
    BL_ASSERT(val != 0);

    const int       M   = D_TERM(2,+2,+4);
    const Real*     plo = geom.ProbLo();
    const Real*     dx  = geom.CellSize();

    Real    fracs[M];
    IntVect cells[M];
    //
    // Get "fracs" and "cells".
    //
    ParticleBase::CIC_Cells_Fracs_Basic(prt, plo, dx, fracs, cells);

    for (int i = 0; i < cnt; i++)
    {
        BL_ASSERT(idx[i] >= 0 && idx[i] < fab.nComp());

        val[i] = ParticleBase::InterpDoit(fab,fracs,cells,idx[i]);
    }
}

void
ParticleBase::GetGravity (const FArrayBox&    gfab,
                          const Geometry&     geom,
                          const ParticleBase& p,
                          Real*               grav)
{
    BL_PROFILE("ParticleBase::GetGravity()");
    BL_ASSERT(grav != 0);

    int idx[BL_SPACEDIM] = { D_DECL(0,1,2) };

    ParticleBase::Interp(p,geom,gfab,idx,grav,BL_SPACEDIM);
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

