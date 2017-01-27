#ifndef _PARTICLES_H_
#define _PARTICLES_H_ 

#include <map>
#include <deque>
#include <vector>
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <array>
#include <memory>
#include <limits>

#include <AMReX_ParmParse.H>

#include <AMReX_ParGDB.H>
#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Array.H>
#include <AMReX_Utility.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_Particles_F.H>
#include <AMReX_RealBox.H>

#ifdef BL_LAZY
#include <AMReX_Lazy.H>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

namespace
{
    std::string aggregation_type   = "";
    int         aggregation_buffer = 1;
    const int   GhostParticleID    = 2000000000;
    const int   VirtualParticleID  = 1000000000;
}

  //
  // A struct used for communicating particle data accross processes
  // during multi-level operations.
  //
struct ParticleCommData
{
#ifdef BL_SINGLE_PRECISION_PARTICLES
    typedef float RealType;
#else
    typedef double RealType;
#endif
  int     m_lev;
  int     m_grid;
  IntVect m_cell;
  RealType m_data[1 + BL_SPACEDIM];
};

  //
  // A struct used for storing a particle's position in the AMR hierarchy.
  //
struct ParticleLocData
{
  int     m_lev;
  int     m_grid;
  IntVect m_cell;
};

template<int NR, int NI=0>
struct Particle
{
    //
    // The floating point type used for the particles.
    //
#ifdef BL_SINGLE_PRECISION_PARTICLES
    typedef float RealType;
#else
    typedef double RealType;
#endif

    //
    // The real data. We always have SPACEDIM position coordinates,
    // and optionally we have NR additional real attributes.
    //
    union rm_t
    {
      RealType pos[BL_SPACEDIM];
      RealType arr[BL_SPACEDIM+NR];
    };

    //
    // The integer data. We always have id and cpu, and optionally we
    // have NI additional integer attributes.
    //
    union im_t
    {
      struct {int id; int cpu;};
      int arr[2+NI];
    };

    rm_t rdata;
    im_t idata;

    static IntVect Index (const Particle<NR, NI>& p, const Geometry& geom);

    //
    // Checks a particles location on levels lev_min and higher.
    // Returns false if the particle does not exist on that level.
    //
    static bool Where (const Particle<NR, NI>& prt, const ParGDBBase* gdb, 
		       ParticleLocData& pld, int lev_min = 0, int finest_level = -1);

    //
    // Checks whether the particle has crossed a periodic boundary in such a way
    // that it is on levels lev_min and higher.
    //
    static bool PeriodicWhere (Particle<NR, NI>& prt, const ParGDBBase* gdb, 
			       ParticleLocData& pld, int lev_min = 0, int finest_level = -1);

    //
    // Checks whether a particle is within its grid (including grow cells).
    //
    static bool RestrictedWhere (const Particle<NR, NI>& p, const ParGDBBase* gdb, 
				 ParticleLocData& pld, int ngrow);

    //
    // Checks a particle's location on a specific level.
    // (Yes this is distict from the functionality provided above)
    //
    static bool SingleLevelWhere (const Particle<NR, NI>& p, const ParGDBBase* gdb, 
				  ParticleLocData& pld, int level);

    //
    // Updates a particle's location (Where), tries to periodic shift any particles
    // that have left the domain. May need work (see inline comments)
    //

    static void Reset (Particle<NR, NI>& prt, const ParGDBBase* gdb, bool update, bool verbose=true); 

    //
    // Returns true if the particle was shifted.
    //
    static bool PeriodicShift (Particle<NR, NI>& prt, const ParGDBBase* gdb);

    static Real InterpDoit (const FArrayBox& fab, const Real* fracs, const IntVect* cells, int comp);

    static Real InterpDoit (const FArrayBox& fab, const IntVect& hi, const Real* frac, int comp);

    static void Interp (const Particle<NR, NI>& prt, const Geometry& geom, const FArrayBox& fab, const int* idx, Real* val, int cnt);

    static const std::string& Version ();

    static const std::string& DataPrefix ();

    static void GetGravity (const FArrayBox& gfab, const Geometry& geom, const Particle<NR, NI>& p, Real* grav);

    static int MaxReaders ();

    static long MaxParticlesPerRead ();

    //
    // Returns the next particle ID for this processor.
    // Particle IDs start at 1 and are never reused.
    // The pair, consisting of the ID and the CPU on which the particle is "born",
    // is a globally unique identifier for a particle.  The maximum of this value
    // across all processors must be checkpointed and then restored on restart
    // so that we don't reuse particle IDs.
    //
    static int NextID ();

    // This version can only be used inside omp critical.
    static int UnprotectedNextID ();

    //
    // Reset on restart.
    //
    static void NextID (int nextid);

    //
    // Used by AssignDensity.
    //
    static bool CrseToFine (const BoxArray&       cfba, 
                            const Array<IntVect>& cells, 
                            Array<IntVect>&       cfshifts, 
                            const Geometry&       gm, 
                            Array<int>&           which, 
                            Array<IntVect>&       pshifts);

    static bool FineToCrse (const Particle<NR, NI>&            p, 
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
                            std::vector< std::pair<int,Box> >& isects);

    static void FineCellsToUpdateFromCrse (const Particle<NR, NI>&            p, 
                                           int lev, const ParGDBBase*         gdb, 
                                           const IntVect&                     ccell,
                                           const IntVect&                     cshift, 
                                           Array<int>&                        fgrid, 
                                           Array<Real>&                       ffrac, 
                                           Array<IntVect>&                    fcells, 
                                           std::vector< std::pair<int,Box> >& isects);

    static void CIC_Fracs (const Real* frac, Real* fracs);

    static void CIC_Cells (const IntVect& hicell, IntVect* cells);

    //
    // Old, *-based CIC for use in Interp.
    //
    static void CIC_Cells_Fracs_Basic (const Particle<NR, NI>& p, const Real* plo, const Real* dx, Real* fracs,  IntVect* cells);

    //
    // Wraps the arbitrary dx function.
    //
    static int CIC_Cells_Fracs (const Particle<NR, NI>& p, 
                                const Real*         plo, 
                                const Real*         dx, 
                                Array<Real>&        fracs,  
                                Array<IntVect>&     cells);
    //
    // Does CIC computations for arbitrary particle/grid dx's.
    //
    static int CIC_Cells_Fracs (const Particle<NR, NI>& p, 
                                const Real*         plo, 
                                const Real*         dx_geom, 
                                const Real*         dx_part, 
                                Array<Real>&        fracs,  
                                Array<IntVect>&     cells);
};

template <int NR, int NI>
void
Particle<NR, NI>::CIC_Cells_Fracs_Basic (const Particle<NR, NI>& p,
					 const Real*         plo,
					 const Real*         dx,
					 Real*               fracs,
					 IntVect*            cells)
{
    BL_PROFILE("Particle<NR, NI>::CIC_Cells_Fracs_B()");
    //
    // "fracs" should be dimensioned: Real    fracs[D_TERM(2,+2,+4)]
    //
    // "cells" should be dimensioned: IntVect cells[D_TERM(2,+2,+4)]
    //
    const Real len[BL_SPACEDIM] = { D_DECL((p.rdata.pos[0]-plo[0])/dx[0] + Real(0.5),
                                           (p.rdata.pos[1]-plo[1])/dx[1] + Real(0.5),
                                           (p.rdata.pos[2]-plo[2])/dx[2] + Real(0.5)) };

    const IntVect cell(D_DECL(floor(len[0]), floor(len[1]), floor(len[2])));

    const Real frac[BL_SPACEDIM] = { D_DECL(len[0]-cell[0], len[1]-cell[1], len[2]-cell[2]) };

    Particle<NR, NI>::CIC_Fracs(frac, fracs);
    Particle<NR, NI>::CIC_Cells(cell, cells);
}

template <int NR, int NI>
int
Particle<NR, NI>::CIC_Cells_Fracs (const Particle<NR, NI>& p,
				   const Real*         plo,
				   const Real*         dx_geom,
				   const Real*         dx_part,
				   Array<Real>&        fracs,
				   Array<IntVect>&     cells)
{
    BL_PROFILE("Particle<NR, NI>::CIC_Cells_Fracs()");
    if (dx_geom == dx_part)
    {
        const int M = D_TERM(2,+2,+4);
        fracs.resize(M);
        cells.resize(M);
        Particle<NR, NI>::CIC_Cells_Fracs_Basic(p,plo,dx_geom,fracs.dataPtr(),cells.dataPtr());
        return M;
    }
    //
    // The first element in fracs and cells is the lowest corner, the last is the highest.
    //
    const Real hilen[BL_SPACEDIM] = { D_DECL((p.rdata.pos[0]-plo[0]+dx_part[0]/2)/dx_geom[0],
                                             (p.rdata.pos[1]-plo[1]+dx_part[1]/2)/dx_geom[1],
                                             (p.rdata.pos[2]-plo[2]+dx_part[2]/2)/dx_geom[2]) };

    const Real lolen[BL_SPACEDIM] = { D_DECL((p.rdata.pos[0]-plo[0]-dx_part[0]/2)/dx_geom[0],
                                             (p.rdata.pos[1]-plo[1]-dx_part[1]/2)/dx_geom[1],
                                             (p.rdata.pos[2]-plo[2]-dx_part[2]/2)/dx_geom[2]) };

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

template <int NR, int NI>
bool
Particle<NR, NI>::CrseToFine (const BoxArray&       cfba,
			      const Array<IntVect>& cells,
			      Array<IntVect>&       cfshifts,
			      const Geometry&       gm,
			      Array<int>&           which,
			      Array<IntVect>&       pshifts)
{
    BL_PROFILE("Particle<NR, NI>::CrseToFine()");
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

template <int NR, int NI>
bool
Particle<NR, NI>::FineToCrse (const Particle<NR, NI>&            p,
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
    BL_PROFILE("Particle<NR, NI>::FineToCrse()");
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

    ParticleLocData pld;
    Particle<NR, NI>::Where(p, gdb, pld);

    const Box& ibx = amrex::grow(gdb->ParticleBoxArray(flev)[pld.m_grid],-1);

    BL_ASSERT(ibx.ok());

    if (ibx.contains(pld.m_cell))
        //
        // We're strictly contained in our valid box.
        // We can't cross a fine->crse boundary.
        //
        return false;

    if (!compfvalid_grown.contains(pld.m_cell))
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

    Particle<NR, NI>::CIC_Cells_Fracs(p, cgm.ProbLo(), cgm.CellSize(), cfracs, ccells);

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

template <int NR, int NI>
void
Particle<NR, NI>::FineCellsToUpdateFromCrse (const Particle<NR, NI>&                p,
					     int                                lev,
					     const ParGDBBase*                  gdb,
					     const IntVect&                     ccell,
					     const IntVect&                     cshift,
					     Array<int>&                        fgrid,
					     Array<Real>&                       ffrac,
					     Array<IntVect>&                    fcells,
					     std::vector< std::pair<int,Box> >& isects)
{
    BL_PROFILE("Particle<NR, NI>::FineCellsToUpdateFromCrse()");
    BL_ASSERT(lev >= 0);
    BL_ASSERT(lev < gdb->finestLevel());

    const Box&      fbx = amrex::refine(Box(ccell,ccell),gdb->refRatio(lev));
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

            if ((p.rdata.pos[k] < celllo) && (celllo > (p.rdata.pos[k] + dx[k]/2)))
                touches = false;

            if ((p.rdata.pos[k] > cellhi) && (cellhi < (p.rdata.pos[k] - dx[k]/2)))
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

            if (p.rdata.pos[k] <= celllo)
            {
                const Real isecthi = p.rdata.pos[k] + dx[k]/2;

                the_frac *= std::min((isecthi - celllo),fdx[k]);
            }
            else
            {
                const Real cellhi  = (iv[k]+1) * fdx[k] + plo[k];
                const Real isectlo = p.rdata.pos[k] - dx[k]/2;

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

template <int NR, int NI>
int
Particle<NR, NI>::MaxReaders ()
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
            amrex::Abort("particles.nreaders must be positive");
    }

    return Max_Readers;
}

template <int NR, int NI>
long
Particle<NR, NI>::MaxParticlesPerRead ()
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
            amrex::Abort("particles.nparts_per_read must be positive");
    }

    return Max_Particles_Per_Read;
}

template <int NR, int NI>
const std::string&
Particle<NR, NI>::DataPrefix ()
{
    //
    // The actual particle data is stored in files of the form: DATA_nnnn.
    //
    static const std::string data("DATA_");

    return data;
}

template <int NR, int NI>
const std::string&
Particle<NR, NI>::Version ()
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

template <int NR, int NI>
int
Particle<NR, NI>::NextID ()
{
    int next;

#ifdef _OPENMP
#pragma omp atomic capture
#endif
    next = the_next_id++;

    if (next == std::numeric_limits<int>::max())
	amrex::Abort("Particle<NR, NI>::NextID() -- too many particles");

    return next;
}

template <int NR, int NI>
int
Particle<NR, NI>::UnprotectedNextID ()
{
    int next = the_next_id++;
    if (next == std::numeric_limits<int>::max())
	amrex::Abort("Particle<NR, NI>::NextID() -- too many particles");
    return next;
}

template <int NR, int NI>
void
Particle<NR, NI>::NextID (int nextid)
{
    the_next_id = nextid;
}

template <int NR, int NI>
IntVect
Particle<NR, NI>::Index (const Particle<NR, NI>& p,
			 const Geometry&     geom)
{
    IntVect iv;

    D_TERM(iv[0]=floor((p.rdata.pos[0]-geom.ProbLo(0))*geom.InvCellSize(0));,
           iv[1]=floor((p.rdata.pos[1]-geom.ProbLo(1))*geom.InvCellSize(1));,
           iv[2]=floor((p.rdata.pos[2]-geom.ProbLo(2))*geom.InvCellSize(2)););

    iv += geom.Domain().smallEnd();

    return iv;
}

template <int NR, int NI>
bool
Particle<NR, NI>::Where (const Particle<NR, NI>&     p,
			 const ParGDBBase*           gdb,
			 ParticleLocData&            pld,
			 int                         lev_min,
			 int                         finest_level)
{
    BL_ASSERT(gdb != 0);

    if (finest_level == -1)
        finest_level = gdb->finestLevel();

    BL_ASSERT(finest_level <= gdb->finestLevel());

    std::vector< std::pair<int, Box> > isects;

    for (int lev = finest_level; lev >= lev_min; lev--)
    {
        const IntVect& iv = Particle<NR, NI>::Index(p, gdb->Geom(lev));

        gdb->ParticleBoxArray(lev).intersections(Box(iv, iv), isects, true, 0);

        if (!isects.empty())
        {
            pld.m_lev  = lev;
            pld.m_grid = isects[0].first;
            pld.m_cell = iv;

            return true;
        }
    }
    return false;
}

template <int NR, int NI>
bool
Particle<NR, NI>::PeriodicWhere (Particle<NR, NI>&           p,
				 const ParGDBBase*           gdb,
				 ParticleLocData&       pld,
				 int                         lev_min,
				 int                         finest_level)
{
    BL_ASSERT(gdb != 0);

    if (!gdb->Geom(0).isAnyPeriodic()) return false;

    if (finest_level == -1)
        finest_level = gdb->finestLevel();

    BL_ASSERT(finest_level <= gdb->finestLevel());

    //
    // Create a copy "dummy" particle to check for periodic outs.
    //
    Particle<NR, NI> p_prime = p;

    if (Particle<NR, NI>::PeriodicShift(p_prime, gdb))
    {
        std::vector< std::pair<int,Box> > isects;

        for (int lev = finest_level; lev >= lev_min; lev--)
        {
            const IntVect& iv = Particle<NR, NI>::Index(p_prime,gdb->Geom(lev));

            gdb->ParticleBoxArray(lev).intersections(Box(iv,iv),isects,true,0);

            if (!isects.empty())
            {
                D_TERM(p.rdata.pos[0] = p_prime.rdata.pos[0];,
                       p.rdata.pos[1] = p_prime.rdata.pos[1];,
                       p.rdata.pos[2] = p_prime.rdata.pos[2];);

                pld.m_lev  = lev;
                pld.m_grid = isects[0].first;
                pld.m_cell = iv;

                return true;
            }
        }
    }

    return false;
}

template <int NR, int NI>
bool
Particle<NR, NI>::RestrictedWhere (const Particle<NR, NI>&     p,
				   const ParGDBBase*           gdb,
				   ParticleLocData&       pld,
				   int                         ngrow)
{
    BL_ASSERT(gdb != 0);

    const IntVect& iv = Particle<NR, NI>::Index(p,gdb->Geom(pld.m_lev));

    if (amrex::grow(gdb->ParticleBoxArray(pld.m_lev)[pld.m_grid], ngrow).contains(iv))
    {
        pld.m_cell = iv;

        return true;
    }

    return false;
}

template <int NR, int NI>
bool 
Particle<NR, NI>::SingleLevelWhere (const Particle<NR, NI>&     p, 
				    const ParGDBBase*     gdb,
				    ParticleLocData& pld,
				    int                   level)
{
    BL_ASSERT(gdb != 0);

    const IntVect& iv = Particle<NR, NI>::Index(p,gdb->Geom(level));

    std::vector< std::pair<int,Box> > isects;

    gdb->ParticleBoxArray(level).intersections(Box(iv,iv),isects,true,0);

    if (!isects.empty())
    {
        pld.m_lev  = level;
        pld.m_grid = isects[0].first;
        pld.m_cell = iv;

        return true;
    }

    return false;
}

template <int NR, int NI>
bool
Particle<NR, NI>::PeriodicShift (Particle<NR, NI>& p,
				 const ParGDBBase*    gdb)
{
    BL_PROFILE("Particle<NR, NI>::PeriodicShift()");
    //
    // This routine should only be called when Where() returns false.
    //
    BL_ASSERT(gdb != 0);
    //
    // We'll use level 0 stuff since ProbLo/ProbHi are the same for every level.
    //
    const Geometry& geom    = gdb->Geom(0);
    const Box&      dmn     = geom.Domain();
    const IntVect&  iv      = Particle<NR, NI>::Index(p,gdb->Geom(0));
    bool            shifted = false;  

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        if (!geom.isPeriodic(i)) continue;

        if (iv[i] > dmn.bigEnd(i))
        {
            if (p.rdata.pos[i] == geom.ProbHi(i))
                //
                // Don't let particles lie exactly on the domain face.
                // Force the particle to be outside the domain so the
                // periodic shift will bring it back inside.
                //
                p.rdata.pos[i] += .125*geom.CellSize(i);

            p.rdata.pos[i] -= geom.ProbLength(i);

            if (p.rdata.pos[i] <= geom.ProbLo(i))
                //
                // This can happen due to precision issues.
                //
                p.rdata.pos[i] += .125*geom.CellSize(i);

            BL_ASSERT(p.rdata.pos[i] >= geom.ProbLo(i));

            shifted = true;
        }
        else if (iv[i] < dmn.smallEnd(i))
        {
            if (p.rdata.pos[i] == geom.ProbLo(i))
                //
                // Don't let particles lie exactly on the domain face.
                // Force the particle to be outside the domain so the
                // periodic shift will bring it back inside.
                //
                p.rdata.pos[i] -= .125*geom.CellSize(i);

            p.rdata.pos[i] += geom.ProbLength(i);

            if (p.rdata.pos[i] >= geom.ProbHi(i))
                //
                // This can happen due to precision issues.
                //
                p.rdata.pos[i] -= .125*geom.CellSize(i);

            BL_ASSERT(p.rdata.pos[i] <= geom.ProbHi(i));

            shifted = true;
        }
    }
    //
    // The particle may still be outside the domain in the case
    // where we aren't periodic on the face out which it travelled.
    //
    return shifted;
}

template <int NR, int NI>
void
Particle<NR, NI>::Reset (Particle<NR, NI>& p,
			 const ParGDBBase*    gdb,
			 bool          update,
			 bool          verbose)
{
    BL_PROFILE("Particle<NR, NI>::Reset()");
    BL_ASSERT(gdb != 0);

    ParticleLocData pld;
    bool ok = Particle<NR, NI>::Where(p, gdb, pld);

    if (!ok && gdb->Geom(0).isAnyPeriodic())
    {
        // Attempt to shift the particle back into the domain if it
        // crossed a periodic boundary.
      Particle<NR, NI>::PeriodicShift(p, gdb, pld);
      ok = Particle<NR, NI>::Where(p, gdb, pld);
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

	BL_ASSERT(p.idata.id > 0);

	p.idata.id = -p.idata.id;
    }
}

template <int NR, int NI>
Real
Particle<NR, NI>::InterpDoit (const FArrayBox& fab,
			      const IntVect&   cell,
			      const Real*      frac,
			      int              comp)
{
    BL_PROFILE("Particle<NR, NI>::InterpDoit(fcfc)");
    const int M = D_TERM(2,+2,+4);

    Real    fracs[M];
    IntVect cells[M];

    Particle<NR, NI>::CIC_Fracs(frac, fracs);
    Particle<NR, NI>::CIC_Cells(cell, cells);

    Real val = Particle<NR, NI>::InterpDoit(fab,fracs,cells,comp);

    return val;
}

template <int NR, int NI>
Real
Particle<NR, NI>::InterpDoit (const FArrayBox& fab,
			      const Real*      fracs,
			      const IntVect*   cells,
			      int              comp)
{
    BL_PROFILE("Particle<NR, NI>::InterpDoit(ffcc)");
    const int M = D_TERM(2,+2,+4);

    Real val = 0;

    for (int i = 0; i < M; i++)
    {
        val += fab(cells[i],comp) * fracs[i];
    }

    return val;
}

template <int NR, int NI>
void
Particle<NR, NI>::Interp (const Particle<NR, NI>& prt,
			  const Geometry&     geom,
			  const FArrayBox&    fab,
			  const int*          idx,
			  Real*               val,
			  int                 cnt)
{
    BL_PROFILE("Particle<NR, NI>::Interp()");
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
    Particle<NR, NI>::CIC_Cells_Fracs_Basic(prt, plo, dx, fracs, cells);

    for (int i = 0; i < cnt; i++)
    {
        BL_ASSERT(idx[i] >= 0 && idx[i] < fab.nComp());

        val[i] = Particle<NR, NI>::InterpDoit(fab,fracs,cells,idx[i]);
    }
}

template <int NR, int NI>
void
Particle<NR, NI>::GetGravity (const FArrayBox&    gfab,
			      const Geometry&     geom,
			      const Particle<NR, NI>& p,
			      Real*               grav)
{
    BL_PROFILE("Particle<NR, NI>::GetGravity()");
    BL_ASSERT(grav != 0);

    int idx[BL_SPACEDIM] = { D_DECL(0,1,2) };

    Particle<NR, NI>::Interp(p,geom,gfab,idx,grav,BL_SPACEDIM);
}

template <int NR, int NI>
std::ostream&
operator<< (std::ostream& os, const Particle<NR, NI>& p)
{
    os << p.idata.id   << ' '
       << p.idata.cpu  << ' ';

    for (int i = 0; i < BL_SPACEDIM + NR; i++)
        os << p.rdata.arr[i] << ' ';

    for (int i = 2; i < 2 + NI; i++)
        os << p.idata.arr[i] << ' ';

    if (!os.good())
        amrex::Error("operator<<(ostream&,Particle<NR, NI>&) failed");

    return os;
}


template <int NR, int NI = 0, class C = std::deque<Particle<NR,NI> > >
class ParticleContainer
{
public:
    //
    // The type of Particles we hold.
    //
    typedef Particle<NR,NI> ParticleType;
    //
    // We want to store the particles on a level by level and grid by grid basis.  This will
    // make accessing them and doing operations on them more memory efficient since most of our
    // operations on particles are done on a level by level basis or grid by grid basis.
    //
    using PBox = C;
    //
    // A level of particles is stored in a map indexed by the grid number.
    //
    using PMap = std::map<int,PBox>;

    ParticleContainer ()
	: m_verbose(1), m_gdb(nullptr), allow_particles_near_boundary(false) {}

    ParticleContainer (ParGDBBase* gdb)
        :
        m_verbose(1), m_gdb(gdb), allow_particles_near_boundary(false) { }

    ParticleContainer (const Geometry            & geom, 
		       const DistributionMapping & dmap,
		       const BoxArray            & ba)
	:
        m_verbose(1), 
	allow_particles_near_boundary(false),
	m_gdb_object(geom,dmap,ba)
    {
	m_gdb = & m_gdb_object;
    }

    ParticleContainer (const Array<Geometry>            & geom, 
		       const Array<DistributionMapping> & dmap,
		       const Array<BoxArray>            & ba,
		       const Array<int>                 & rr)
	:
        m_verbose(1), 
	allow_particles_near_boundary(false),
	m_gdb_object(geom,dmap,ba,rr)
    {
	m_gdb = & m_gdb_object;
    }

    ~ParticleContainer () {}

    void Define (ParGDBBase* gdb)
    {
	m_gdb = gdb;
    }

    void Define (const Geometry            & geom, 
		 const DistributionMapping & dmap,
		 const BoxArray            & ba)
    {
	m_gdb_object = ParGDB(geom, dmap, ba);
	m_gdb = &m_gdb_object;
    }

    void Define (const Array<Geometry>            & geom, 
		 const Array<DistributionMapping> & dmap,
		 const Array<BoxArray>            & ba,
		 const Array<int>                 & rr)
    {
	m_gdb_object = ParGDB(geom, dmap, ba, rr);
	m_gdb = &m_gdb_object;
    }

    void SetParticleBoxArray (int lev, const BoxArray& new_ba)
	{ m_gdb->SetParticleBoxArray(lev, new_ba); }

    void SetParticleDistributionMap (int lev, const DistributionMapping& new_dmap)
	{ m_gdb->SetParticleDistributionMap(lev, new_dmap); }

    const BoxArray& ParticleBoxArray (int lev) const 
	{ return m_gdb->ParticleBoxArray(lev); }

    const DistributionMapping& ParticleDistributionMap (int lev) const 
	{ return m_gdb->ParticleDistributionMap(lev); }

    const ParGDBBase* GetParGDB () const { return m_gdb; }

    void InitFromAsciiFile (const std::string& file, int extradata, const IntVect* Nrep = 0);
    void InitFromBinaryFile (const std::string& file, int extradata);
    void InitFromBinaryMetaFile (const std::string& file, int extradata);
    void InitRandom (long icount, unsigned long iseed, Real particleMass, bool serialize = false, RealBox bx = RealBox());
    void InitOnePerCell (Real x_off, Real y_off, Real z_off,
                         Real particleMass, MultiFab& particle_mf);
    void InitNRandomPerCell (int n_per_cell, Real particleMass, MultiFab& particle_mf);

    void addOneParticle (int id_in, int cpu_in, 
                         std::vector<Real>& xloc, std::vector<Real>& attributes); 
    void addNParticles (int n, D_DECL(Array<Real>& x, Array<Real>& y, Array<Real>& z), 
                               D_DECL(Array<Real>& u, Array<Real>& v, Array<Real>& w), 
                        int na, Array<Real>& attributes);

    void GetParticleIDs        (Array<int> & part_ids);
    void GetParticleCPU        (Array<int> & part_cpu);
    void GetParticleLocations  (Array<Real>& part_locs);
    void GetParticleData       (Array<Real>& part_data, int start_comp, int num_comp);

    void SetParticleLocations  (Array<Real>& part_data);

    void MoveRandom ();
    void MoveRandom (int level);

    void Increment (MultiFab& mf, int level);
    long IncrementWithTotal (MultiFab& mf, int level);

    // rho_index: rho index in rdata
    Real sumParticleMass (int rho_index, int level) const;

    // Set the flag that allows particles to live near the domain boundary and throw away 
    //     the part of their contribution in AssignDensity that is outside the domain.
    void SetAllowParticlesNearBoundary(bool value);
 
    void Redistribute (bool where_already_called = false,
		       bool full_where           = false, // w.z.: why is default not true?
		       int  lev_min              = 0, 
		       int  nGrow                = 0);

    void RedistributeMPI (PMap& not_ours);
    //
    // OK checks that all particles are in the right places (for some value of right)
    //
    // These flags are used to do proper checking for subcycling particles
    // the default values are fine for non-subcycling methods
    //
    bool OK (bool full_where = false, int lev_min = 0 , int ngrow = 0, int finest_level = -1) const;

    void ByteSpread () const;
    //
    // Returns # of particles at specified the level.
    //
    // If "only_valid" is true it only counts valid particles.
    //
    long NumberOfParticlesAtLevel (int level, bool only_valid = true, bool only_local = false) const;
    Array<long> NumberOfParticlesInGrid  (int level, bool only_valid = true, bool only_local = false) const;
    //
    // Returns # of particles at all levels
    //
    // If "only_valid" is true it only counts valid particles.
    //
    long TotalNumberOfParticles (bool only_valid=true, bool only_local=false) const;

    //
    // The Following methods are for managing Virtual and Ghost Particles.
    //
    // Removes all particles at a given level
    //
    void RemoveParticlesAtLevel (int level);
    // 
    void RemoveParticlesNotAtFinestLevel ();
    //
    // Creates virtual particles for a given level that represent
    // in some capacity all particles at finer levels
    //
    void CreateVirtualParticles (int level, PBox& virts) const;
    // 
    // Create ghost particles for a given level that are copies of particles
    // near coarse->fine boundaries in level-1
    //
    void CreateGhostParticles (int level, int ngrow, PBox& ghosts) const;
    //
    // Add particles from a pbox to the grid at this level
    //
    void AddParticlesAtLevel (int level, PBox& virts, bool where_already_called = false);

    void Checkpoint (const std::string& dir, const std::string& name, bool is_checkpoint = true) const;

    void Restart (const std::string& dir, const std::string& file, bool is_checkpoint = true);

    void WritePlotFile (const std::string& dir, const std::string& name) const;

    void WriteAsciiFile          (const std::string& file);

    int Verbose () { return m_verbose; }

    void SetVerbose (int verbose) { m_verbose = verbose; }

    const PMap& GetParticles(int lev) const { return m_particles[lev]; }
    PMap& GetParticles(int lev) { return m_particles[lev]; }

    // 
    // Functions depending the layout of the data.  Use with caution.
    //
    // Multi-level version
    void AssignDensity (int rho_index, bool sub_cycle, Array<std::unique_ptr<MultiFab> >& mf, 
			int lev_min = 0, int ncomp = 1, int finest_level = -1) const;
    // Single-level version
    void AssignDensitySingleLevel (int rho_index, MultiFab& mf, int level,
				   int ncomp=1, int particle_lvl_offset = 0) const;
    void AssignCellDensitySingleLevel (int rho_index, MultiFab& mf, int level,
				       int ncomp=1, int particle_lvl_offset = 0) const;
    void NodalDepositionSingleLevel   (int rho_index, MultiFab& mf, int level,
				       int ncomp=1, int particle_lvl_offset = 0) const;
    //
    void moveKick (MultiFab& acceleration, int level, Real timestep, 
		   Real a_new = 1.0, Real a_half = 1.0,
		   int start_comp_for_accel = -1);

protected:

    bool OnSameGrids (int level, const MultiFab& mf) const { return m_gdb->OnSameGrids(level, mf); }

    //
    // Helper function for Checkpoint() and WritePlotFile().
    //
    void WriteParticles (int            level,
                         std::ofstream& ofs,
                         int            fnum,
                         Array<int>&    which,
                         Array<int>&    count,
                         Array<long>&   where,
                         bool           is_checkpoint) const;
    //
    // Helper functions for Restart().
    //
    void Restart_Doit (const std::string& fullname,
                       std::ifstream&     HdrFile,
                       const std::string& how,
                       bool           is_checkpoint);

    void ReadParticles_DoublePrecision (int            cnt,
                                        int            grd,
                                        int            lev,
                                        bool           is_checkpoint,
                                        std::ifstream& ifs);

    void ReadParticles_SinglePrecision (int            cnt,
                                        int            grd,
                                        int            lev,
                                        bool           is_checkpoint,
                                        std::ifstream& ifs);
    //
    // The data.
    //
    int         m_verbose;


    ParGDBBase* m_gdb;
    bool allow_particles_near_boundary;
    ParGDB      m_gdb_object;
    Array<PMap> m_particles;

private:
    void AssignDensityDoit (int rho_index, 
			    Array<std::unique_ptr<MultiFab> >& mf, 
			    std::map<int, std::deque<ParticleCommData> >& data,
			    int ncomp, int lev_min = 0) const;
};

template <int NR, int NI, class C>
long
ParticleContainer<NR,NI,C>::TotalNumberOfParticles (bool only_valid, bool only_local) const
{
    long nparticles = 0;
    for (int lev = 0; lev <= m_gdb->finestLevel(); lev++) {
        nparticles += NumberOfParticlesAtLevel(lev,only_valid,false);
    }
    if (!only_local) {
	ParallelDescriptor::ReduceLongSum(nparticles);
    }
    return nparticles;
}

template <int NR, int NI, class C>
Array<long>
ParticleContainer<NR,NI,C>::NumberOfParticlesInGrid (int  lev, 
						   bool only_valid,
						   bool only_local) const
{
    int ngrids = m_gdb->ParticleBoxArray(lev).size();
    Array<long> nparticles(ngrids,0);

    if (lev >= 0 && lev < int(m_particles.size()))
    {
        const PMap& pmap = m_particles[lev];

        for (const auto& kv : pmap)
        {
            int         pgrd = kv.first;
            const PBox& pbox = kv.second;

            if (only_valid)
            {
		for (const auto& p : pbox)
		{
		    if (p.idata.id > 0) ++nparticles[pgrd];
		}
            }
            else
            {
                nparticles[pgrd] = pbox.size();
            }
        }

        if (!only_local)
            ParallelDescriptor::ReduceLongSum(&nparticles[0],ngrids);
    }

    return nparticles;
}

template <int NR, int NI, class C>
long
ParticleContainer<NR,NI,C>::NumberOfParticlesAtLevel (int  lev,
						      bool only_valid,
						      bool only_local) const
{
    long nparticles = 0;

    if (lev >= 0 && lev < int(m_particles.size()))
    {
        const PMap& pmap = m_particles[lev];

	for (const auto& kv : pmap)
        {
            const PBox& pbox = kv.second;

            if (only_valid)
            {
		for (const auto& p : pbox)
		{
		    if (p.idata.id > 0) ++nparticles;
		}
            }
            else
            {
                nparticles += pbox.size();
            }
        }
    }

    if (!only_local)
        ParallelDescriptor::ReduceLongSum(nparticles);

    return nparticles;
}

//
// This includes both valid and invalid particles since invalid particles still take up space.
//

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::ByteSpread () const
{
    long cnt = 0;

    for (const auto& pmap : m_particles) {
	for (const auto& kv : pmap) {
	    cnt += kv.second.size();
	}
    }

    long mn = cnt, mx = mn;

    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    const std::size_t sz = sizeof(ParticleType);

#ifdef BL_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
    ParallelDescriptor::ReduceLongMin(mn, IOProc);
    ParallelDescriptor::ReduceLongMax(mx, IOProc);
    ParallelDescriptor::ReduceLongSum(cnt,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "ParticleContainer<NR,NI,C> byte spread across MPI nodes: ["
                  << mn*sz
		  << " (" << mn << ")"
                  << " ... "
                  << mx*sz
		  << " (" << mx << ")"
                  << "] total particles: (" << cnt << ")\n";
    }
#ifdef BL_LAZY
    });
#endif
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::addOneParticle (int                  id_in,
					    int                  cpu_in, 
					    std::vector<double>& xloc, 
					    std::vector<double>& attributes)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::addOneParticle()");
    if (m_particles.size() == 0)
    {
       m_particles.resize(m_gdb->finestLevel()+1);
    }

    ParticleType p;

    p.idata.id  = id_in;
    p.idata.cpu = cpu_in;

    if (p.idata.id <= 0)
        amrex::Abort("Particle ID's must be > 0 in addOneParticle");
 
    if (ParallelDescriptor::MyProc() != p.idata.cpu)
        amrex::Abort("cpu_in must equal MyProc() in addOneParticle");

    for (int i = 0; i < BL_SPACEDIM; i++)
       p.rdata.pos[i] = xloc[i];

    for (int i = BL_SPACEDIM; i < BL_SPACEDIM + NR; i++)
       p.rdata.arr[i] = attributes[i];
    //
    // It's possible that particles are defined to live right on the boundary,
    // so if that's the case we move them ever so slightly inside the domain 
    // instead.
    //
    const Geometry& geom = m_gdb->Geom(0);

    const Real  Delta[BL_SPACEDIM] = { D_DECL(Real(.125)*geom.CellSize(0),
                                              Real(.125)*geom.CellSize(1),
                                              Real(.125)*geom.CellSize(2)) };
    //
    // If particle is right on a domain boundary then move it just inside the boundary.
    //
    for (int d = 0; d < BL_SPACEDIM; d++)
    {
        if (p.rdata.pos[d] <= geom.ProbLo(d)) p.rdata.pos[d] += Delta[d];
        if (p.rdata.pos[d] >= geom.ProbHi(d)) p.rdata.pos[d] -= Delta[d];
    }

    ParticleLocData pld;
    if (!ParticleType::Where(p, m_gdb, pld))
    {
        ParticleType::PeriodicShift(p, m_gdb);

        if (!ParticleType::Where(p, m_gdb, pld))
        {
            for (int d = 0; d < BL_SPACEDIM; d++)
            {
                std::cout << "BAD PARTICLE POS(" << d << ") " << p.rdata.pos[d] << std::endl;
            }
            amrex::Abort("ParticleContainer<NR,NI,C>::addOneParticle(): invalid particle");
        }
    }

    m_particles[pld.m_lev][pld.m_grid].push_back(p);

    // Note that we will need to call Redistribute once we are done adding particles this way.
    // The Where call above assigns a particle to a grid (pld.m_grid) based on the particle location.
    // However, the particle may not currently live on the processor that owns that grid.
    // The Redistribute routine should ensure that the particle ends up on the right processor.
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::addNParticles(int n_part, Array<double>& x, Array<double>& y, 
#if (BL_SPACEDIM == 3)
                                                      Array<double>& z, 
#endif
                                                      Array<double>& u, Array<double>& v, 
#if (BL_SPACEDIM == 3)
                                                      Array<double>& w, 
#endif
                                          int n_attr, Array<double>& attributes)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::addNParticles()");

    BL_ASSERT(OK());
    BL_ASSERT(NR >= BL_SPACEDIM+n_attr);

    if (m_particles.size() == 0)
       m_particles.resize(m_gdb->finestLevel()+1);

    ParticleType p;

    for (long j = 0; j < n_part; j++)
    {
       ParticleType p;

       p.rdata.pos[j  ] = x[j]; // x
       p.rdata.pos[j+1] = y[j]; // y
#if (BL_SPACEDIM == 3)
       p.rdata.pos[j+2] = z[j]; // z
#endif

       p.rdata.arr[BL_SPACEDIM + j]     = u[j]; // u
       p.rdata.arr[BL_SPACEDIM + j + 1] = v[j]; // v
#if (BL_SPACEDIM == 3)
       p.rdata.arr[BL_SPACEDIM + j + 2] = w[j]; // w
#endif

       ParticleLocData pld;
       if (!ParticleType::Where(p, m_gdb, pld))
              amrex::Abort("ParticleContainer<NR,NI,C>::addNParticles(): invalid particle location");

       for (int i = 0; i < n_attr; i++)
          p.rdata.arr[BL_SPACEDIM+i] = attributes[n_attr*j + i];

       p.idata.id  = ParticleType::NextID();
       p.idata.cpu = ParallelDescriptor::MyProc();

       m_particles[pld.m_lev][pld.m_grid].push_back(p);
    }

    Redistribute();
}

//
// The ParticleContainer<NR,NI,C>::Init... routines have been moved into a separate file.
//
#include "AMReX_ParticleInit.H"

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::MoveRandom ()
{
    //
    // Move particles randomly at all levels
    //
    for (int lev = 0; lev < int(m_particles.size()); lev++)
    {
       MoveRandom(lev);
    }
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::MoveRandom (int lev)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::MoveRandom(lev)");
    BL_ASSERT(OK());
    BL_ASSERT(m_gdb != 0);
    // 
    // Move particles up to FRAC*CellSize distance in each coordinate direction.
    //
    const Real FRAC = 0.25;

    static bool first = true;

    static Array<amrex::mt19937> rn;

    if (first)
    {
        first = false;
        //
        // Build and initialize a random number generator per thread.
        //
        int tnum = 1;

#ifdef _OPENMP
        tnum = omp_get_max_threads();
#endif
        rn.resize(tnum);

        for (int i = 0; i < tnum; i++)
        {
            //
            // We want to give each thread across all MPI processes a unique non-zero seed.
            //
            const unsigned long seedbase = 1+tnum*ParallelDescriptor::MyProc();

            rn[i] = amrex::mt19937(seedbase+i);
        }
    }

    PMap&       pmap              = m_particles[lev];
    const Real* dx                = m_gdb->Geom(lev).CellSize();
    const Real  dist[BL_SPACEDIM] = { D_DECL(FRAC*dx[0], FRAC*dx[1], FRAC*dx[2]) };

    for (auto& kv : pmap)
    {
        PBox&     pbox = kv.second;
        const int n    = pbox.size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.idata.id <= 0) continue;

#ifdef _OPENMP
            int tid = omp_get_thread_num();
#else
            int tid = 0;
#endif
            for (int i = 0; i < BL_SPACEDIM; i++)
            {
                p.rdata.pos[i] += dist[i]*(2*rn[tid].d_value()-1);
            }

            ParticleType::Reset(p,m_gdb,true);
        }
    }

    Redistribute(true);
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::Increment (MultiFab& mf,
				     int       lev) 
{
    IncrementWithTotal(mf,lev);
}

template <int NR, int NI, class C>
long
ParticleContainer<NR,NI,C>::IncrementWithTotal (MultiFab& mf,
						int       lev)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::IncrementWithTotal(lev)");
    BL_ASSERT(OK());

    if (m_particles.empty()) return 0;

    BL_ASSERT(lev >= 0 && lev < int(m_particles.size()));

    const PMap& pmap = m_particles[lev];
  
    long num_particles_in_domain = 0;

    MultiFab* mf_pointer;

    if (OnSameGrids(lev, mf))
    {
        // If we are already working with the internal mf defined on the
        // particle_box_array, then we just work with this.
        mf_pointer = &mf;
    }
    else
    {
        // If mf is not defined on the particle_box_array, then we need
        // to make a temporary mf_pointer here and copy it into mf at the end.
        mf_pointer = new MultiFab(m_gdb->ParticleBoxArray(lev),
				  m_gdb->ParticleDistributionMap(lev),
				  mf.nComp(),mf.nGrow());
				  
    }

    for (const auto& kv : pmap)
    {
        const int   grid = kv.first;
        const PBox& pbox = kv.second;
        FArrayBox&  fab  = (*mf_pointer)[grid];

	for (const auto& p : pbox)
        {
	  ParticleLocData pld;
	  Particle<NR, NI>::Where(p, m_gdb, pld);

            if (p.idata.id > 0)
            {
                BL_ASSERT(pld.m_grid == grid);

                fab(pld.m_cell) += 1;
                num_particles_in_domain += 1;
            }
        }
    }

    // If mf is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf.   I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf) 
    {
	mf.copy(*mf_pointer,0,0,mf.nComp());  
	delete mf_pointer;
    }

    ParallelDescriptor::ReduceLongSum(num_particles_in_domain);

    return num_particles_in_domain;
}

template <int NR, int NI, class C>
Real
ParticleContainer<NR,NI,C>::sumParticleMass (int rho_index, int lev) const
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::sumParticleMass(lev)");
    BL_ASSERT(NR >= 1);
    BL_ASSERT(lev >= 0 && lev < int(m_particles.size()));

    Real msum = 0;

    const PMap& pmap = m_particles[lev];

    for (const auto& kv : pmap)
    {
        const PBox& pbox = kv.second;

	for (const auto& p : pbox)
        {
            if (p.idata.id > 0)
            {
                msum += p.rdata.arr[BL_SPACEDIM+rho_index];
            }
        }
    }

    ParallelDescriptor::ReduceRealSum(msum);

    return msum;
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::GetParticleIDs (Array<int>& part_ids)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::GetParticleIDs()");
    //
    // This gives us the starting point into the part_ids array
    // If only one processor (or no MPI), then that's all we need.
    //
    int cnt = 0;

#if BL_USE_MPI
    Array<long> cnts(ParallelDescriptor::NProcs());

    // This returns the number of particles on this processor
    long lcnt = TotalNumberOfParticles(true,true);

    // This accumulates the "lcnt" values into "cnts"
    MPI_Gather(&lcnt,1,              
               ParallelDescriptor::Mpi_typemap<long>::type(),
               cnts.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<long>::type(),
               ParallelDescriptor::IOProcessorNumber(),
               ParallelDescriptor::Communicator());

    ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());

    for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); iproc++) {
        cnt += cnts[iproc];
    }

    std::cout << "PROC CNT " << ParallelDescriptor::MyProc() << " " << cnt << std::endl;

    // Each particle takes up 1 int so no need to multiply cnt by anything
#endif

    // This is the total number of particles on *all* processors
    long npart = TotalNumberOfParticles(true,false);

    // Locations
    part_ids.resize(npart,0);

    for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
    {
        const PMap& pmap = m_particles[lev];

	for (const auto& kv : pmap)
        {
            const PBox& pbx = kv.second;
            const int   n   = pbx.size();
    
	    for (const auto& p : pbx)
            {
                if (p.idata.id > 0)
                {
                    // Load the ID
                    part_ids[cnt] = p.idata.id;

                    // Update counter
                    cnt++;
                }
            }
        }
    }

    ParallelDescriptor::ReduceIntSum(part_ids.dataPtr(),part_ids.size()); 
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::GetParticleCPU (Array<int>& part_cpu)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::GetParticleCPU()");
    //
    // This gives us the starting point into the part_cpu array
    // If only one processor (or no MPI), then that's all we need.
    //
    int cnt = 0;

#if BL_USE_MPI
    Array<long> cnts(ParallelDescriptor::NProcs());

    // This returns the number of particles on this processor
    long lcnt = TotalNumberOfParticles(true,true);

    // This accumulates the "lcnt" values into "cnts"
    MPI_Gather(&lcnt,1,              
               ParallelDescriptor::Mpi_typemap<long>::type(),
               cnts.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<long>::type(),
               ParallelDescriptor::IOProcessorNumber(),
               ParallelDescriptor::Communicator());

    ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());

    for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); iproc++)
        cnt += cnts[iproc];

    // Each particle takes up 1 int so no need to multiply cnt by anything
#endif

    // This is the total number of particles on *all* processors
    long npart = TotalNumberOfParticles(true,false);

    // Locations
    part_cpu.resize(npart,0);

    for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
    {
        const PMap& pmap = m_particles[lev];

	for (const auto& kv : pmap)
        {
            const PBox& pbx = kv.second;
            const int   n   = pbx.size();
    
	    for (const auto& p : pbx)
            {
                if (p.idata.id > 0)
                {
                    // Load the ID
                    part_cpu[cnt] = p.idata.cpu;

                    // Update counter
                    cnt++;
                }
            }
        }
    }

    ParallelDescriptor::ReduceIntSum(part_cpu.dataPtr(),part_cpu.size()); 
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::GetParticleLocations (Array<Real>& part_data)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::GetParticleLocations()");
    //
    // This gives us the starting point into the part_data array
    // If only one processor (or no MPI), then that's all we need.
    //
    int cnt = 0;

#if BL_USE_MPI
    Array<long> cnts(ParallelDescriptor::NProcs());

    // This returns the number of particles on this processor
    long lcnt = TotalNumberOfParticles(true,true);

    // This accumulates the "lcnt" values into "cnts"
    MPI_Gather(&lcnt,1,              
               ParallelDescriptor::Mpi_typemap<long>::type(),
               cnts.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<long>::type(),
               ParallelDescriptor::IOProcessorNumber(),
               ParallelDescriptor::Communicator());

    ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());

    for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); iproc++)
        cnt += cnts[iproc];

    // Each particle takes up BL_SPACEDIM Reals
    cnt *= (BL_SPACEDIM);
#endif

    // This is the total number of particles on *all* processors
    long npart = TotalNumberOfParticles(true,false);

    // Locations
    part_data.resize(BL_SPACEDIM*npart,0);

    for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
    {
        const PMap& pmap = m_particles[lev];

	for (const auto& kv : pmap)
        {
            const PBox& pbx = kv.second;
            const int   n   = pbx.size();
    
	    for (const auto& p : pbx)
	    {
                if (p.idata.id > 0)
                {
                    // Load positions
                    for (int d=0; d < BL_SPACEDIM; d++)
                        part_data[cnt+d] = p.rdata.pos[d];

                    // Update counter
                    cnt += BL_SPACEDIM;
                }
            }
        }
    }

    ParallelDescriptor::ReduceRealSum(part_data.dataPtr(),part_data.size()); 
}


template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::GetParticleData (Array<Real>& part_data, int start_comp, int num_comp)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::GetParticleData()");
    //
    // This gives us the starting point into the part_data array
    // If only one processor (or no MPI), then that's all we need.
    //
    int cnt = 0;

    //
    // Make sure we don't try to get more than we have
    //
    if (start_comp + num_comp > BL_SPACEDIM + NR)
        amrex::Error("Tried to grab too many components in GetParticleData!!");

#if BL_USE_MPI
    Array<long> cnts(ParallelDescriptor::NProcs());

    // This returns the number of particles on this processor
    long lcnt = TotalNumberOfParticles(true,true);

    // This accumulates the "lcnt" values into "cnts"
    MPI_Gather(&lcnt,1,              
               ParallelDescriptor::Mpi_typemap<long>::type(),
               cnts.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<long>::type(),
               ParallelDescriptor::IOProcessorNumber(),
               ParallelDescriptor::Communicator());

    ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());

    for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); iproc++)
        cnt += cnts[iproc];

    // Each particle takes up num_comp Reals
    cnt*= num_comp;
#endif

    // This is the total number of particles on *all* processors
    long npart = TotalNumberOfParticles(true,false);

    part_data.resize(num_comp*npart,0);

    for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
    {
        const PMap& pmap = m_particles[lev];

	for (const auto& kv : pmap)
        {
            const PBox& pbx = kv.second;
            const int   n   = pbx.size();
    
	    for (const auto& p : pbx)
	    {
                if (p.idata.id > 0)
                {
                    // Load particle data, whatever it is.
                    for (int d = 0; d < num_comp; d++)
                      part_data[cnt+d] = p.rdata.arr[BL_SPACEDIM + start_comp + d];

                    // Update counter
                    cnt += num_comp;
                }
            }
        }
    }

    ParallelDescriptor::ReduceRealSum(part_data.dataPtr(),part_data.size()); 
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::SetAllowParticlesNearBoundary (bool value)
{
    allow_particles_near_boundary = value; 
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::SetParticleLocations (Array<Real>& part_data)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::SetParticleLocations()");
   // This gives us the starting point into the part_data array
   // If only one processor (or no MPI), then that's all we need
   int cnt = 0;

#if BL_USE_MPI
   Array<long> cnts(ParallelDescriptor::NProcs());

   // This returns the number of particles on this processor
   long lcnt = TotalNumberOfParticles(true,true);

   // This accumulates the "lcnt" values into "cnts"
   MPI_Gather(&lcnt,1,              
              ParallelDescriptor::Mpi_typemap<long>::type(),
              cnts.dataPtr(),
              1,
              ParallelDescriptor::Mpi_typemap<long>::type(),
              ParallelDescriptor::IOProcessorNumber(),
              ParallelDescriptor::Communicator());

   ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());

   for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); iproc++)
       cnt += cnts[iproc];

   // Each particle takes up BL_SPACEDIM Reals
   cnt*= BL_SPACEDIM;
#endif

   // This is the total number of particles on *all* processors
   long npart = TotalNumberOfParticles(true,false);

   // Mass + locations
   if (part_data.size() != npart*BL_SPACEDIM)
       amrex::Abort("Sending in wrong size part_data to SetParticleLocations");

   for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
   {
       PMap& pmap = m_particles[lev];

       for (auto& kv : pmap)
       {
           PBox&     pbx = kv.second;
           const int n   = pbx.size();
    
	   for (auto& p : pbx)
           {
              if (p.idata.id > 0)
              {
                  // Load positions
                  for (int d=0; d < BL_SPACEDIM; d++)
                     p.rdata.pos[d] = part_data[cnt+d];

                  // Update counter
                  cnt += BL_SPACEDIM;
              }
           }
       }
    }
}


template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::AddParticlesAtLevel (int   level,
					       PBox& virts,
					       bool  where_already_called)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::AddParticlesAtLevel()");
    if (int(m_particles.size()) < level+1)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<NR,NI,C>::AddParticlesAtLevel resizing m_particles from "
                      << m_particles.size()
                      << " to "
                      << level+1 << '\n';
        }
        m_particles.resize(level + 1);
    }

    const int MyProc = ParallelDescriptor::MyProc();
    //
    // The valid particles that we don't own.
    //
    PMap not_ours;

    while (!virts.empty())
    {
        ParticleType& p = virts.back();
	ParticleLocData pld;

        if (p.idata.id > 0)
        {
            if (!where_already_called)
            {
                //
                // Put the particle in this level.
                //
                p.m_lev = level;

                if (!ParticleType::SingleLevelWhere(p, m_gdb, pld, level))
                    //
                    // Virtuals shouldn't be in Ghost cells.
                    //
                    amrex::Abort("ParticleContainer<NR,NI,C>::AddParticlesAtLevel(): Can't add outside of domain\n");
            }
            else
            {
                BL_ASSERT(pld.m_lev == level);
            }

            const int who = m_gdb->ParticleDistributionMap(pld.m_lev)[pld.m_grid];

            if (who == MyProc)
            {
                m_particles[pld.m_lev][pld.m_grid].push_back(p);
            }
            else
            {
                not_ours[who].push_back(p);
            }
        }

        virts.pop_back();
    }

    if (ParallelDescriptor::NProcs() == 1)
    {
        BL_ASSERT(not_ours.empty());
    }
    else
    {
        RedistributeMPI(not_ours);
    }
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::RemoveParticlesAtLevel (int level)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::RemoveParticlesAtLevel()");
    if (level >= int(this->m_particles.size()))
        return;

    if (!this->m_particles[level].empty())
    {
        PMap().swap(this->m_particles[level]);
    }
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::RemoveParticlesNotAtFinestLevel ()
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::RemoveParticlesNotAtFinestLevel()");
    BL_ASSERT(this->m_gdb->finestLevel()+1 == int(this->m_particles.size()));

    int cnt = 0;

    for (int lev = 0; lev < this->m_gdb->finestLevel(); lev++)
    {
        PMap& pmap = this->m_particles[lev];

        if (!pmap.empty())
        {
	    for (const auto& kv : pmap)
            {
                cnt += kv.second.size();
            }

            PMap().swap(pmap);
        }
    }
    //
    // Print how many particles removed on each processor if any were removed.
    //
    if (this->m_verbose > 1)
    {
        int maxcnt = cnt;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceIntMax(maxcnt);

        if (maxcnt > 0)
        {
            for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
            {
                if (ParallelDescriptor::MyProc() == i)
                {
                    if (cnt > 0)
                    {
                        std::cout << "Processor "
                                  << i
                                  << " removed "
                                  << cnt
                                  << " particles not in finest level" << std::endl;
                    }
                }
                ParallelDescriptor::Barrier();
            }
        }
#ifdef BL_LAZY
	});
#endif
    }
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::CreateVirtualParticles (int   level,
						    PBox& virts) const
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::CreateVirtualParticles()");
    BL_ASSERT(level > 0);
    BL_ASSERT(virts.empty());

    if (level >= int(m_particles.size()))
        //
        // This level could exist and simply have no particles.
        //
        return;
    //
    // Read these from the parm file if we haven't done so yet.
    //
    if (aggregation_type == "")
    {
        ParmParse pp("particles");
        aggregation_type = "None";
        pp.query("aggregation_type",aggregation_type);
        aggregation_buffer = 2;
        pp.query("aggregation_buffer",aggregation_buffer);
    }
    //
    // Create a buffer so that particles near the cf border are not aggregated.
    //
    BoxArray buffer = amrex::complementIn(m_gdb->Geom(level).Domain(), m_gdb->ParticleBoxArray(level));

    buffer.grow(aggregation_buffer);

    const PMap& pmap = m_particles[level];

    for (const auto& kv : pmap)
    {
        const PBox& pbox = kv.second;
        //
        // Map for use in Cell aggregation.
        //
        std::map<IntVect,ParticleType> agg_map;

	for (auto it = pbox.cbegin(); it != pbox.cend(); ++it)
        {
	  ParticleLocData pld;
	  Particle<NR, NI>::Where(*it, m_gdb, pld);
            if (buffer.contains(pld.m_cell))
            {
                //
                // It's in the no-aggregation buffer.
                //
                virts.push_back(*it);
                //
                // Set its id to indicate that it's a virt.
                //
                virts.back().idata.id = VirtualParticleID;
            }
            else
            {
                if (aggregation_type == "None")
                {
                    //
                    // No aggregation.  Simply clone the particle.
                    //
                    virts.push_back(*it);
                    //
                    // Set its id to indicate that it's a virt.
                    //
                    virts.back().idata.id = VirtualParticleID;
                }
                else if (aggregation_type == "Cell")
                {
                    //
                    // Note that Cell aggregation assumes that p.rdata.arr[BL_SPACEDIM] is mass and
                    // that all other components should be combined in a mass-weighted
                    // average.
                    //
                    auto agg_map_it = agg_map.find(pld.m_cell);

                    if (agg_map_it == agg_map.end())
                    {
                        //
                        // Add the particle.
                        //
                        ParticleType p = *it;
                        //
                        // Set its id to indicate that it's a virt.
                        //
                        p.idata.id = VirtualParticleID;
                        agg_map[pld.m_cell] = p;
                    }
                    else
                    {
                        BL_ASSERT(agg_map_it != agg_map.end());
                        const ParticleType&  pnew       = *it;
                        ParticleType&        pold       = agg_map_it->second;
                        const Real           old_mass   = pold.rdata.arr[BL_SPACEDIM];
                        const Real           new_mass   = pnew.rdata.arr[BL_SPACEDIM];
                        const Real           total_mass = old_mass + new_mass;
                        //
                        // Set the position to the center of mass.
                        //
                        for (int i = 0; i < BL_SPACEDIM; i++)
                        {
                            pold.rdata.pos[i] = (old_mass*pold.rdata.pos[i] + new_mass*pnew.rdata.pos[i])/total_mass;
                        }
                        BL_ASSERT(ParticleType::Index(pold,m_gdb->Geom(level)) == pld.m_cell);
                        //
                        // Set the metadata (presumably velocity) to the mass-weighted average.
                        //
                        for (int i = BL_SPACEDIM + 1; i < BL_SPACEDIM + NR; i++)
                        {
                            pold.rdata.arr[i] = (old_mass*pold.rdata.arr[i] + new_mass*pnew.rdata.arr[i])/total_mass;
                        }
                        pold.rdata.arr[BL_SPACEDIM] = total_mass;
                    }
                }
                else if (aggregation_type == "Flow")
                {
                    amrex::Abort("Flow aggregation not implemented");
                }
                else 
                {
                    amrex::Abort("Unknown Particle Aggregation mode");
                }
            }
        }
        if (aggregation_type == "Cell")
        {
            //
            // Add the aggregated particles to the virtuals.
            //
	    for (const auto& kv : agg_map)
            {
                virts.push_back(kv.second);
            }
        }
    }
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::CreateGhostParticles (int   level,
						int   ngrow,
						PBox& ghosts) const
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::CreateGhostParticles()");
    BL_ASSERT(ghosts.empty());
    BL_ASSERT(level < m_gdb->finestLevel());

    if (level >= int(m_particles.size()))
        //
        // This level could exist and simply have no particles.
        //
        return;

    const BoxArray& fine = m_gdb->ParticleBoxArray(level + 1);
    
    std::vector< std::pair<int,Box> > isects;

    const PMap& pmap = m_particles[level];

    for (const auto& kv : pmap)
    {
        const PBox& pbox = kv.second;

	for (auto it = pbox.cbegin(); it != pbox.cend(); ++it)
        {
            //
            // Find particle location on the finer level.
            //
            const IntVect& iv = ParticleType::Index(*it,m_gdb->Geom(level+1));
            //
            // Is it in the grown finer level?
            //
            fine.intersections(Box(iv,iv),isects,false,ngrow);
            //
            // Here we add the particle to each potential grid.
            //
	    for (const auto& isec : isects)
            {
                //
                // Create a copy.
                //
                ParticleType p = *it;
		ParticleLocData pld;
		
		ParticleType::SingleLevelWhere(p, m_gdb, pld, level + 1);

                //
                // Set its id to indicate that it's a ghost.
                //
                p.idata.id = GhostParticleID;

                //
                // Store it in the PBox.
                //
                ghosts.push_back(p);
            }
        }
    }
}
    
//
// This redistributes valid particles and discards invalid ones.
//

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::Redistribute (bool where_already_called,
					  bool full_where,
  					  int  lev_min,
					  int  nGrow)
{
    BL_PROFILE("ParticleContainer::Redistribute()");
    const int MyProc   = ParallelDescriptor::MyProc();
    Real      strttime = ParallelDescriptor::second();

    // Note - fix this later. It doesn't make sense to have this option
    // now that the particles don't store the lev, grid, and cell.
    where_already_called = false;

    //
    // On startup there are cases where Redistribute() could be called
    // with a given finestLevel() where that AmrLevel has yet to be defined.
    //
    int theEffectiveFinestLevel = m_gdb->finestLevel();

    while (!m_gdb->LevelDefined(theEffectiveFinestLevel))
        theEffectiveFinestLevel--;

    if (int(m_particles.size()) < theEffectiveFinestLevel+1)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<NR,NI,C>::Redistribute() resizing m_particles from "
                      << m_particles.size()
                      << " to "
                      << theEffectiveFinestLevel+1 << '\n';
        }
        m_particles.resize(theEffectiveFinestLevel+1);
    }
    //
    // The valid particles that we don't own.
    //
    PMap not_ours;

    for (int lev = lev_min, nlevs = m_particles.size(); lev < nlevs; lev++)
    {
        PMap& pmap = m_particles[lev];

	for (auto pmap_it = pmap.begin(); pmap_it != pmap.end(); /* no ++ */ )
        {
            const int grid = pmap_it->first;
            PBox&     pbox = pmap_it->second;

	    auto first = pbox.begin();
	    auto last  = pbox.end();

	    if (first != last) 
	    {
		for (auto it = first; it != last; ++it)
		{
		  ParticleLocData pld;
		  ParticleType& p = *it;

		  if (p.idata.id > 0)
		    {
		      if (!where_already_called)
			{
			  if (!ParticleType::Where(p, m_gdb, pld, lev_min, theEffectiveFinestLevel))
			    {                                
			      if (full_where) // Lengthier checks for subcycling.
				{
				  if (!ParticleType::PeriodicWhere(p, m_gdb, pld, lev_min, theEffectiveFinestLevel))
				    {
				      if (lev_min != 0) // RestrictedWhere should be unnecessary at top level.
					{
					  if (!ParticleType::RestrictedWhere(p, m_gdb, pld, nGrow))
					    amrex::Abort("ParticleContainer<NR,NI,C>::Redistribute(): invalid particle at non-coarse step");
					}
				      else
					{
					  //
					  // The particle has left the domain; invalidate it.
					  // This typically only happens on a coarse timestep.
					  //
					  p.idata.id = -p.idata.id;
					}
				    }
				}
			      else
				{
				  std::cout << "Bad Particle: " << p << '\n';
				  amrex::Abort("ParticleContainer<NR,NI,C>::Redistribute(): invalid particle in basic check");
				}
			    }
			}

			if (p.idata.id > 0)
			{
			    //
			    // The owner of the particle is the CPU owning the finest grid
			    // in state data that contains the particle.
			    //
			    const int who = m_gdb->ParticleDistributionMap(pld.m_lev)[pld.m_grid];
			    
			    if (who == MyProc)
			    {
				if (pld.m_lev != lev || pld.m_grid != grid)
				{
				    //
				    // We own it but must shift it to another place.
				    //
				    m_particles[pld.m_lev][pld.m_grid].push_back(p);
				    //
				    // Invalidate the particle so we can reclaim its space.
				    //
				    p.idata.id = -p.idata.id;
				}
			    }
			    else
			    {
				not_ours[who].push_back(p);
				//
				// Invalidate the particle so we can reclaim its space.
				//
				p.idata.id = -p.idata.id;
			    }
			}
		    }

		    if (p.idata.id > 0) // this is a valid particle
		    {
			if (it != first) *first = p;
			++first;
		    }
		}
		pbox.erase(first, last);
	    }
	    //
	    // Remove any map entries for which the particle container is now empty.
	    //
            if (pmap_it->second.empty())
            {
                pmap.erase(pmap_it++);
            }
            else
            {
                ++pmap_it;
            }
        }
    }

    if (int(m_particles.size()) > theEffectiveFinestLevel+1)
    {
        //
        // Looks like we lost an AmrLevel on a regrid.
        //
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<NR,NI,C>::Redistribute() resizing m_particles from "
                      << m_particles.size()
                      << " to "
                      << theEffectiveFinestLevel+1 << '\n';
        }
        BL_ASSERT(int(m_particles.size()) >= 2);
        BL_ASSERT(m_particles[m_particles.size()-1].empty());

        m_particles.resize(theEffectiveFinestLevel+1);
    }

    if (ParallelDescriptor::NProcs() == 1)
    {
        BL_ASSERT(not_ours.empty());
    }
    else
    {
        RedistributeMPI(not_ours);
    }

    BL_ASSERT(OK(full_where, lev_min, nGrow, theEffectiveFinestLevel));

    if (m_verbose > 0)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ByteSpread();

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
        if (ParallelDescriptor::IOProcessor())
            std::cout << "ParticleContainer<NR,NI,C>::Redistribute() time: " << stoptime << "\n\n";
#ifdef BL_LAZY
	});
#endif
    }
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::RedistributeMPI (PMap& not_ours)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::RedistributeMPI()");
#if BL_USE_MPI
    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();
    //
    // We may now have particles that are rightfully owned by another CPU.
    //
    Array<int> Snds(NProcs,0), Rcvs(NProcs,0);

    int NumSnds = 0, NumRcvs = 0;

    for (const auto& kv : not_ours)
    {
        NumSnds       += kv.second.size();
        Snds[kv.first] = kv.second.size();
    }

    ParallelDescriptor::ReduceIntMax(NumSnds);

    if (NumSnds == 0)
        //
        // There's no parallel work to do.
        //
        return;

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());

    BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 Rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 ParallelDescriptor::Communicator()) );
    BL_ASSERT(Rcvs[MyProc] == 0);

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    typedef std::map<int,int> IntIntMap;

    IntIntMap SndCnts, RcvCnts, rOffset;

    for (int i = 0; i < NProcs; i++)
        if (Snds[i] > 0)
            SndCnts[i] = Snds[i];

    for (int i = 0; i < NProcs; i++)
    {
        if (Rcvs[i] > 0)
        {
            RcvCnts[i] = Rcvs[i];
            rOffset[i] = NumRcvs;
            NumRcvs   += Rcvs[i];
        }
    }
    //
    // Don't need these anymore.
    //
    Array<int>().swap(Snds);
    Array<int>().swap(Rcvs);
    //
    // We'll store the particles we're to receive in a PMap indexed by proc # of receiver.
    //
    PMap nparticles;

    for (const auto& kv : RcvCnts)
    {
        nparticles[kv.first].resize(kv.second);
    }
    Array<int>         owner(RcvCnts.size());
    Array<int>         index(RcvCnts.size());
    Array<MPI_Status>  stats(RcvCnts.size());
    Array<MPI_Request> rreqs(RcvCnts.size());
    //
    // First send/recv the integer parts of the particles.
    //
    {
        const int SeqNum     = ParallelDescriptor::SeqNum();
        const int iChunkSize = 2 + NI;
        //
        // Allocate data for rcvs as one big chunk.
        //
        Array<int> recvdata(NumRcvs * iChunkSize);
        //
        // Post receives.
        //
        int idx = 0;
	for (auto it = RcvCnts.cbegin(); it != RcvCnts.cend(); ++it, ++idx)
        {
            const int Who = it->first;
            const int Cnt = it->second   * iChunkSize;
            const int Idx = rOffset[Who] * iChunkSize;

            BL_ASSERT(Cnt > 0);
            BL_ASSERT(Who >= 0 && Who < NProcs);
            BL_ASSERT(Cnt < std::numeric_limits<int>::max());

            owner[idx] = Who;
            rreqs[idx] = ParallelDescriptor::Arecv(&recvdata[Idx],Cnt,Who,SeqNum).req();
        }
        //
        // Send the integer data.
        //
        Array<int> senddata;

	for (const auto& kv : SndCnts)
        {
            const int Who = kv.first;
            const int Cnt = kv.second * iChunkSize;

            BL_ASSERT(Cnt > 0);
            BL_ASSERT(Who >= 0 && Who < NProcs);
            BL_ASSERT(Cnt < std::numeric_limits<int>::max());

            senddata.resize(Cnt);

            const PBox& pbox = not_ours[Who];

            int ioff = 0;
	    for (const auto& p : pbox)
            {
                BL_ASSERT(p.idata.id > 0);

                senddata[ioff+0] = p.idata.id;
                senddata[ioff+1] = p.idata.cpu;

 		for (int i = 0; i < NI; ++i) {
		    senddata[i+ioff+2] = p.idata.arr[i+2];
		}

                ioff += iChunkSize;
            }

            ParallelDescriptor::Send(senddata.dataPtr(),Cnt,Who,SeqNum);
        }
        //
        // Free up this memory ...
        //
        Array<int>().swap(senddata);
        //
        // Now receive and unpack the integer data.
        //
        for (int NWaits = rreqs.size(), completed; NWaits > 0; NWaits -= completed)
        {
            ParallelDescriptor::Waitsome(rreqs, completed, index, stats);

            for (int k = 0; k < completed; k++)
            {
                const int  Who  = owner[index[k]];
                const int  Idx  = rOffset[Who] * iChunkSize;
                const int* rcvp = &recvdata[Idx];
                PBox&      pbox = nparticles[Who];

                BL_ASSERT(int(pbox.size()) == RcvCnts[Who]);

		for (auto& p : pbox)
                {
		  BL_ASSERT(rcvp != 0);

		  p.idata.id   = rcvp[0];
		  p.idata.cpu  = rcvp[1];

		  for (int i = 0; i < NI; ++i) {
		    p.idata.arr[i+2] = rcvp[i+2];
		  }

		  rcvp += iChunkSize;
                }
            }
        }
    }
    //
    // Next send/recv the Real parts of the particles.
    //
    {
        const int SeqNum     = ParallelDescriptor::SeqNum();
        const int rChunkSize = BL_SPACEDIM+NR;
        //
        // Allocate data for rcvs as one big chunk.
        //
        Array<typename ParticleType::RealType> recvdata(NumRcvs * rChunkSize);
        //
        // Post receives.
        //
        int idx = 0;
	for (auto it = RcvCnts.cbegin(); it != RcvCnts.cend(); ++it, ++idx)
        {
            const int Who = it->first;
            const int Cnt = it->second   * rChunkSize;
            const int Idx = rOffset[Who] * rChunkSize;

            BL_ASSERT(Cnt > 0);
            BL_ASSERT(Who >= 0 && Who < NProcs);
            BL_ASSERT(Cnt < std::numeric_limits<int>::max());

            rreqs[idx] = ParallelDescriptor::Arecv(&recvdata[Idx],Cnt,Who,SeqNum).req();
        }
        //
        // Send the Real data.
        //
        Array<typename ParticleType::RealType> senddata;

	for (const auto& kv : SndCnts)
        {
            const int Who = kv.first;
            const int Cnt = kv.second * rChunkSize;

            BL_ASSERT(Cnt > 0);
            BL_ASSERT(Who >= 0 && Who < NProcs);
            BL_ASSERT(Cnt < std::numeric_limits<int>::max());
            
            senddata.resize(Cnt);

            PBox& pbox = not_ours[Who];

            int ioff = 0;
	    for (const auto& p : pbox)
            {
                BL_ASSERT(p.idata.id > 0);

                D_TERM(senddata[ioff+0] = p.rdata.pos[0];,
                       senddata[ioff+1] = p.rdata.pos[1];,
                       senddata[ioff+2] = p.rdata.pos[2];);

                ioff += BL_SPACEDIM;

                for (int j = 0; j < NR; j++)
                    senddata[ioff+j] = p.rdata.arr[BL_SPACEDIM + j];

                ioff += NR;
            }

	    PBox().swap(pbox);

            ParallelDescriptor::Send(senddata.dataPtr(),Cnt,Who,SeqNum);
        }
        //
        // Free up this memory ...
        //
        Array<typename ParticleType::RealType>().swap(senddata);
        //
        // Now receive and unpack the Real data.
        //
        for (int NWaits = rreqs.size(), completed; NWaits > 0; NWaits -= completed)
        {
            ParallelDescriptor::Waitsome(rreqs, completed, index, stats);

            for (int k = 0; k < completed; k++)
            {
                const int                              Who  = owner[index[k]];
                const int                              Idx  = rOffset[Who] * rChunkSize;
                const typename ParticleType::RealType* rcvp = &recvdata[Idx];
                PBox&                                  pbox = nparticles[Who];

                BL_ASSERT(int(pbox.size()) == RcvCnts[Who]);

		ParticleLocData pld;
		for (auto& p : pbox)
                {
                    BL_ASSERT(rcvp != 0);

                    D_TERM(p.rdata.pos[0] = rcvp[0];,
                           p.rdata.pos[1] = rcvp[1];,
                           p.rdata.pos[2] = rcvp[2];);

                    rcvp += BL_SPACEDIM;

                    for (int j = 0; j < NR; j++)
                        p.rdata.arr[BL_SPACEDIM + j] = rcvp[j];

                    rcvp += NR;

		    ParticleType::Where(p, m_gdb, pld);

                    m_particles[pld.m_lev][pld.m_grid].push_back(p);
                }

		PBox().swap(pbox);
            }
        }
    }
#endif /*BL_USE_MPI*/
}

template <int NR, int NI, class C>
bool
ParticleContainer<NR,NI,C>::OK (bool full_where,
				int  lev_min,
				int  ngrow,
				int  finest_level) const
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::OK()");
    if (finest_level == -1)
        finest_level = m_gdb->finestLevel();

    BL_ASSERT(finest_level <= m_gdb->finestLevel());
    //
    // Check that the integer data in each valid particle is what it should be.
    // This includes checking that particles are in the proper place in the particle
    // container based on what Where() says they should be.
    //
    // Particles are copied to avoid accidentally moving them with where.
    //
    for (int lev = lev_min, nlevs=m_particles.size(); lev < nlevs; lev++)
    {
        const PMap& pmap = m_particles[lev];

	for (const auto& kv : pmap)
        {
            const int   grid = kv.first;
            const PBox& pbox = kv.second;

	    for (auto it = pbox.cbegin(); it != pbox.cend(); ++it)
            {
                //
                // Yes I want to make a copy of the particle.
                //
                ParticleType p = *it;

                if (p.idata.id > 0)
                {
		  
		  ParticleLocData pld;
		  Particle<NR, NI>::Where(p, m_gdb, pld);

		  const int     llev  = pld.m_lev;
		  const int     lgrid = pld.m_grid;
		  const IntVect cell  = pld.m_cell;

		  if (!ParticleType::Where(p, m_gdb, pld, lev_min, finest_level))
                    {
		      if (full_where)
                        {
			  if (!ParticleType::PeriodicWhere(p, m_gdb, pld, lev_min, finest_level)) 
                            {
			      if (!ParticleType::RestrictedWhere(p, m_gdb, pld, ngrow))
				return false;
                            }
                        }
		      else
                        {
			  return false;
                        }
                    }
                    if ((lev  != pld.m_lev  || lev  != llev)  ||
                        (grid != pld.m_grid || grid != lgrid) || cell != pld.m_cell)
                    {
                        std::cout << "PARTICLE NUMBER " << p.idata.id << '\n';

                        std::cout << "POS IS ";
                        for (int i = 0; i < BL_SPACEDIM; i++)
                            std::cout << p.rdata.pos[i] << ' ';

                        if (lev != pld.m_lev || lev != llev)
                           std::cout << "BAD LEV  " << lev  << " " << pld.m_lev << '\n';

                        if (grid != pld.m_grid || grid != lgrid)
                           std::cout << "BAD GRID " << grid << " " << pld.m_grid << '\n';

                        if (cell != pld.m_cell)
                           std::cout << "BAD CELL " << cell << " " << pld.m_cell << '\n';

                        return false;
                    }
                }
            }
        }
    }

    return true;
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::Checkpoint (const std::string& dir,
					const std::string& name,
					bool               is_checkpoint) const
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::Checkpoint()");
    BL_ASSERT(OK());

    BL_ASSERT(sizeof(typename ParticleType::RealType) == 4 || sizeof(typename ParticleType::RealType) == 8);

    const int  MyProc   = ParallelDescriptor::MyProc();
    const int  NProcs   = ParallelDescriptor::NProcs();
    const int  IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real strttime = ParallelDescriptor::second();
    //
    // We store the particles in a subdirectory of "dir".
    //
    std::string pdir = dir;

    if (!pdir.empty() && pdir[pdir.size()-1] != '/')
        pdir += '/';

    pdir += name;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(pdir, 0755))
            amrex::CreateDirectoryFailed(pdir);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();
    //
    // The header contains the info we need to read back in the particles.
    //
    // Only the I/O processor writes to the header file.
    //
    std::ofstream HdrFile;

    long nparticles = 0;

    for (const auto& pmap : m_particles)
    {
	for (const auto& kv : pmap)
        {
            const PBox& pbox = kv.second;

	    for (const auto& p : pbox)
            {
                if (p.idata.id > 0)
                    //
                    // Only count (and checkpoint) valid particles.
                    //
                    nparticles++;
            }
        }
    }

    ParallelDescriptor::ReduceLongSum(nparticles,IOProc);

    int maxnextid = ParticleType::NextID();

    ParticleType::NextID(maxnextid);

    ParallelDescriptor::ReduceIntMax(maxnextid,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        std::string HdrFileName = pdir;

        if (!HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/')
            HdrFileName += '/';

        HdrFileName += "Header";

        HdrFile.open(HdrFileName.c_str(), std::ios::out|std::ios::trunc);

        if (!HdrFile.good())
            amrex::FileOpenFailed(HdrFileName);
        //
        // First thing written is our Checkpoint/Restart version string.
        // 
        // We append "_single" or "_double" to the version string indicating
        // whether we're using "float" or "double" floating point data in the
        // particles so that we can Restart from the checkpoint files.
        //
        if (sizeof(typename ParticleType::RealType) == 4)
        {
            HdrFile << ParticleType::Version() << "_single" << '\n';
        }
        else
        {
            HdrFile << ParticleType::Version() << "_double" << '\n';
        }
        //
        // BL_SPACEDIM and N for sanity checking.
        //
        HdrFile << BL_SPACEDIM << '\n';

        HdrFile << NR << '\n';
        //
        // The total number of particles.
        //
        HdrFile << nparticles << '\n';
        //
        // The value of nextid that we need to restore on restart.
        //
        HdrFile << maxnextid << '\n';
        //
        // Then the finest level of the AMR hierarchy.
        //
        HdrFile << m_gdb->finestLevel() << '\n';
        //
        // Then the number of grids at each level.
        //
        for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
        {
            HdrFile << m_gdb->ParticleBoxArray(lev).size() << '\n';
        }
    }
    //
    // We want to write the data out in parallel.
    //
    // We'll allow up to nOutFiles active writers at a time.
    //
    int nOutFiles(64);
    ParmParse pp("particles");
    pp.query("particles_nfiles",nOutFiles);
    if(nOutFiles == -1) {
      nOutFiles = NProcs;
    }
    nOutFiles = std::max(1, std::min(nOutFiles,NProcs));

    for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
    {
        const bool gotsome = (NumberOfParticlesAtLevel(lev) > 0);
        //
        // We store the particles at each level in their own subdirectory.
        //
        std::string LevelDir = pdir;

        if (gotsome)
        {
            if (!LevelDir.empty() && LevelDir[LevelDir.size()-1] != '/')
                LevelDir += '/';

            LevelDir = amrex::Concatenate(LevelDir + "Level_", lev, 1);

            if (ParallelDescriptor::IOProcessor())
                if (!amrex::UtilCreateDirectory(LevelDir, 0755))
                    amrex::CreateDirectoryFailed(LevelDir);
            //
            // Force other processors to wait till directory is built.
            //
            ParallelDescriptor::Barrier();
        }

	MFInfo info;
	info.SetAlloc(false);
	MultiFab state(m_gdb->ParticleBoxArray(lev),
		       m_gdb->ParticleDistributionMap(lev),
		       1,0,info);
        //
        // We eventually want to write out the file name and the offset
        // into that file into which each grid of particles is written.
        //
        Array<int>  which(state.size(),0);
        Array<int > count(state.size(),0);
        Array<long> where(state.size(),0);

        if (gotsome)
        {
            const int   FileNumber   = MyProc % nOutFiles;
            std::string FullFileName = LevelDir;

            FullFileName += '/';
            FullFileName += ParticleType::DataPrefix();
            FullFileName += amrex::Concatenate("", FileNumber, 4);

            std::ofstream ParticleFile;

            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

            ParticleFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

            const int nSets = ((NProcs + (nOutFiles - 1)) / nOutFiles);
            const int mySet = (MyProc / nOutFiles);

            for (int iSet = 0; iSet < nSets; ++iSet)
            {
                if (mySet == iSet)
                {
                    //
                    // Write all the data at this level to the file.
                    //
                    if (iSet == 0)
                        //
                        // First set.
                        //
                        ParticleFile.open(FullFileName.c_str(),
                                          std::ios::out|std::ios::trunc|std::ios::binary);
                    else
                    {
                        ParticleFile.open(FullFileName.c_str(),
                                          std::ios::out|std::ios::app|std::ios::binary);
                        //
                        // Set to the end of the file.
                        //
                        ParticleFile.seekp(0, std::ios::end);
                    }

                    if (!ParticleFile.good())
                        amrex::FileOpenFailed(FullFileName);
                    //
                    // Write out all the valid particles we own at the specified level.
                    // Do it grid block by grid block remembering the seek offset
                    // for the start of writing of each block of data.
                    //
                    WriteParticles(lev, ParticleFile, FileNumber, which, count, where, is_checkpoint);

                    ParticleFile.flush();

                    ParticleFile.close();

                    if (!ParticleFile.good())
                        amrex::Abort("ParticleContainer<NR,NI,C>::Checkpoint(): problem writing ParticleFile");

                    int iBuff = 0, wakeUpPID = (MyProc + nOutFiles), tag = (MyProc % nOutFiles);

                    if (wakeUpPID < NProcs)
                    {
                        ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
                    }
                }

                if (mySet == (iSet + 1))
                {
                    //
                    // Next set waits.
                    //
                    int iBuff, waitForPID = (MyProc - nOutFiles), tag = (MyProc % nOutFiles);

                    ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
                }
            }

            ParallelDescriptor::ReduceIntSum (which.dataPtr(), which.size(), IOProc);
            ParallelDescriptor::ReduceIntSum (count.dataPtr(), count.size(), IOProc);
            ParallelDescriptor::ReduceLongSum(where.dataPtr(), where.size(), IOProc);
        }

        if (ParallelDescriptor::IOProcessor())
        {
            for (int j = 0; j < state.size(); j++)
            {
                //
                // We now write the which file, the particle count, and the
                // file offset into which the data for each grid was written,
                // to the header file.
                //
                HdrFile << which[j] << ' ' << count[j] << ' ' << where[j] << '\n';
            }

            if (gotsome)
            {
                //
                // Unlink any zero-length data files.
                //
                Array<long> cnt(nOutFiles,0);

                for (int i = 0, N=count.size(); i < N; i++)
                    cnt[which[i]] += count[i];

                for (int i = 0, N=cnt.size(); i < N; i++)
                {
                    if (cnt[i] == 0)
                    {
                        std::string FullFileName = LevelDir;

                        FullFileName += '/';
                        FullFileName += ParticleType::DataPrefix();
                        FullFileName += amrex::Concatenate("", i, 4);

                        amrex::UnlinkFile(FullFileName.c_str());
                    }
                }
            }
        }
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            HdrFile.flush();

            HdrFile.close();

            if (!HdrFile.good())
                amrex::Abort("ParticleContainer<NR,NI,C>::Checkpoint(): problem writing HdrFile");

            std::cout << "ParticleContainer<NR,NI,C>::Checkpoint() time: " << stoptime << '\n';
        }
    }
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::WritePlotFile (const std::string& dir,
					 const std::string& name) const
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::WritePlotFile()");
    BL_ASSERT(OK());
    bool is_checkpoint = false;

    // For yt we need exactly the chk particle format so would need to set is_checkpoint = true
    // Anyway, it's not too bad to have particle ids on disk,
    // think of merger trees or backtracing of particles for nested ics
    // is_checkpoint = true; 
    Checkpoint(dir,name,is_checkpoint);
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::WriteParticles (int            lev,
					  std::ofstream& ofs,
					  int            fnum,
					  Array<int>&    which,
					  Array<int>&    count,
					  Array<long>&   where,
					  bool           is_checkpoint) const
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::WriteParticles()");
    const PMap&     pmap  = m_particles[lev];

    MFInfo info;
    info.SetAlloc(false);
    MultiFab state(m_gdb->ParticleBoxArray(lev),
		   m_gdb->ParticleDistributionMap(lev),
		   1,0,info);

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const int grid = mfi.index();
        //
        // Only write out valid particles.
        //
        int cnt = 0;

        auto pmap_it = pmap.find(grid);

        if (pmap_it != pmap.end())
        {
	    for (const auto& p : pmap_it->second)
            {
                if (p.idata.id > 0)
                    cnt++;
            }
        }

        which[grid] = fnum;
        count[grid] = cnt;
        where[grid] = VisMF::FileOffset(ofs);

        if (cnt == 0) continue;

        const PBox& pbox = pmap_it->second;

        if (is_checkpoint)
        {
            //
            // First write out the integer data in binary.
            // We do not need to write out the m_lev and m_grid
            // info since it's implicit in how the particles
            // are stored.  We can easily recreate them on restart.
            //
            const int iChunkSize = 2;

#ifdef BL_LOWMEMPWRITE
	    int maxItemsToWrite(8192);
	    int cntBufSize(maxItemsToWrite * iChunkSize);
	    int nItems(cnt), nItemsToWrite(0);
            Array<int> istuff(cntBufSize);
            auto it  = pbox.cbegin();
            auto End = pbox.cend();

	    while(nItems > 0) {
              int *iptr = istuff.dataPtr();
	      int itemCount(0);
              for( ; it != End && itemCount < maxItemsToWrite; ++it) {
                if(it->idata.id > 0) {

                    iptr[0] = it->idata.id;
                    iptr[1] = it->idata.cpu;

                    iptr += iChunkSize;
                }
		++itemCount;
              }
	      nItemsToWrite = nItems > maxItemsToWrite ? maxItemsToWrite : nItems;
              ofs.write((char *) istuff.dataPtr(), nItemsToWrite * iChunkSize * sizeof(int));
	      nItems -= nItemsToWrite;
	    }
#else
            Array<int> istuff(cnt*iChunkSize);

            int* iptr = istuff.dataPtr();

	    for (auto it = pbox.cbegin(); it != pbox.cend(); ++it)
            {
                if (it->idata.id > 0)
                {
                    iptr[0] = it->idata.id;
                    iptr[1] = it->idata.cpu;

                    iptr += iChunkSize;
                }
            }

            ofs.write((char*)istuff.dataPtr(),istuff.size()*sizeof(int));
#endif
        }


        //
        // Write the Real data in binary.
        //
        const int rChunkSize = BL_SPACEDIM+NR;

#ifdef BL_LOWMEMPWRITE
	int maxItemsToWrite(8192);
	int cntBufSize(maxItemsToWrite * rChunkSize);
	int nItems(cnt), nItemsToWrite(0);
        Array<ParticleType::RealType> rstuff(cntBufSize);
        auto it  = pbox.cbegin();
        auto End = pbox.cend();

	while(nItems > 0) {
          ParticleType::RealType *rptr = rstuff.dataPtr();
	  int itemCount(0);
          for( ; it != End && itemCount < maxItemsToWrite; ++it) {
            if(it->idata.id > 0) {
                for (int i = 0; i < BL_SPACEDIM + NR; i++)
                  rptr[i] = it->m.rdata.arr[i];
                rptr += rChunkSize;
            }
	    ++itemCount;
          }

	  nItemsToWrite = nItems > maxItemsToWrite ? maxItemsToWrite : nItems;
          ofs.write((char *) rstuff.dataPtr(), nItemsToWrite * rChunkSize * sizeof(typename ParticleType::RealType));
	  nItems -= nItemsToWrite;
	}
#else
        Array<typename ParticleType::RealType> rstuff(cnt*rChunkSize);

        typename ParticleType::RealType* rptr = rstuff.dataPtr();

	for (auto it = pbox.cbegin(); it != pbox.cend(); ++it)
        {
            if (it->idata.id > 0)
            {
	      for (int i = 0; i < BL_SPACEDIM + NR; i++)
		rptr[i] = it->m.rdata.arr[i];
	      rptr += rChunkSize;
            }
        }

        ofs.write((char*)rstuff.dataPtr(),rstuff.size()*sizeof(typename ParticleType::RealType));
#endif
    }
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::Restart (const std::string& dir,
				     const std::string& file,
				     bool is_checkpoint)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::Restart()");
    BL_ASSERT(!dir.empty());
    BL_ASSERT(!file.empty());

    const int  IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real strttime = ParallelDescriptor::second();

    std::string fullname = dir;

    if (!fullname.empty() && fullname[fullname.size()-1] != '/')
        fullname += '/';

    fullname += file;
    //
    // The header contains the info we need to read back in the particles.
    //
    // Only the IO processor reads the header file.
    //
    // It'll then broadcast() stuff of interest to all CPUs.
    //
    std::ifstream HdrFile;

    std::string HdrFileName = fullname;

    if (!HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/')
        HdrFileName += '/';

    HdrFileName += "Header";

    HdrFile.open(HdrFileName.c_str(), std::ios::in);

    if (!HdrFile.good())
        amrex::FileOpenFailed(HdrFileName);
    //
    // First value should be the version string.
    //
    Array<char> vbuf(128);

    std::string version;

    if (ParallelDescriptor::IOProcessor())
    {
        HdrFile >> version;

        BL_ASSERT(!version.empty());
        BL_ASSERT(vbuf.size() > version.size());

        for (size_t i = 0; i < version.size(); ++i)
            vbuf[i] = version[i];

        vbuf[version.size()] = '\0';
    }

    ParallelDescriptor::Bcast(vbuf.dataPtr(), vbuf.size(), IOProc);
    //
    // What do our version strings mean?
    //
    // "Version_One_Dot_Zero" -- hard-wired to write out in double precision.
    // 
    // "Version_One_Dot_One" -- can write out either as either single or double precision.
    //
    // Appended to the latter version string are either "_single" or "_double" to
    // indicate how the particles were written.
    //
    version = vbuf.dataPtr();

    if (version.find("Version_One_Dot_Zero") != std::string::npos)
    {
        Restart_Doit(fullname,HdrFile,"double",is_checkpoint);
    }
    else if (version.find("Version_One_Dot_One") != std::string::npos)
    {
        if (version.find("_single") != std::string::npos)
        {
            Restart_Doit(fullname,HdrFile,"single",is_checkpoint);
        }
        else if (version.find("_double") != std::string::npos)
        {
            Restart_Doit(fullname,HdrFile,"double",is_checkpoint);
        }
        else
        {
            std::string msg("ParticleContainer<NR,NI,C>::Restart(): bad version string: ");
            msg += version;
            amrex::Error(version.c_str());
        }
    }
    else
    {
        std::string msg("ParticleContainer<NR,NI,C>::Restart(): unknown version string: ");
        msg += version;
        amrex::Abort(msg.c_str());
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<NR,NI,C>::Restart() time: " << stoptime << '\n';
        }
    }
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::Restart_Doit (const std::string& fullname,
					std::ifstream&     HdrFile,
					const std::string& how,
					bool is_checkpoint)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::RestartDoit()");
    BL_ASSERT(!fullname.empty());

    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    //
    // Next value should be BL_SPACEDIM;
    //
    int dm;

    if (ParallelDescriptor::IOProcessor())
    {
        HdrFile >> dm;

        if (dm != BL_SPACEDIM)
            amrex::Abort("ParticleContainer<NR,NI,C>::Restart(): dm != BL_SPACEDIM");
    }
    ParallelDescriptor::Bcast(&dm, 1, IOProc);
    //
    // Next value should be our "N".
    //
    int n;

    if (ParallelDescriptor::IOProcessor())
    {
        HdrFile >> n;

        if (n != NR)
            amrex::Abort("ParticleContainer<NR,NI,C>::Restart(): n != N");
    }
    ParallelDescriptor::Bcast(&n, 1, IOProc);

    long nparticles;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // The total number of particles.
        //
        HdrFile >> nparticles;

        BL_ASSERT(nparticles >= 0);
    }
    ParallelDescriptor::Bcast(&nparticles, 1, IOProc);

    int maxnextid;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // The value of nextid that we need to restore.
        //
        HdrFile >> maxnextid;

        BL_ASSERT(maxnextid > 0);
    }
    ParallelDescriptor::Bcast(&maxnextid, 1, IOProc);
    //
    // Don't forget to restore it!!!
    //
    ParticleType::NextID(maxnextid);
    //
    // Then the finest level of the AMR hierarchy.
    //
    int finest_level;

    if (ParallelDescriptor::IOProcessor())
    {
        HdrFile >> finest_level;

        BL_ASSERT(finest_level >= 0);
    }
    ParallelDescriptor::Bcast(&finest_level, 1, IOProc);
    //
    // Then the number of grids at each level.
    //
    Array<int> ngrids(finest_level+1);

    BL_ASSERT(finest_level == m_gdb->finestLevel());

    if (ParallelDescriptor::IOProcessor())
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            HdrFile >> ngrids[lev];

            BL_ASSERT(ngrids[lev] > 0);
            BL_ASSERT(ngrids[lev] == int(m_gdb->ParticleBoxArray(lev).size()));
        }
    }
    ParallelDescriptor::Bcast(ngrids.dataPtr(), ngrids.size(), IOProc);
    //
    // The rest of HdrFile consists of triples of the form:
    //
    //   which count offset
    //
    // One for each grid at each level from 0 -> finest_level.
    //
    // We rebuild the filename from which and level.
    //
    for (int lev = 0; lev <= finest_level; lev++)
    {
        //
        // Read in the which, count & offset info for this level.
        //
        Array<int>  which(ngrids[lev]);
        Array<int>  count(ngrids[lev]);
        Array<long> where(ngrids[lev]);

        if (ParallelDescriptor::IOProcessor())
        {
            for (int i = 0; i < ngrids[lev]; i++)
            {
                HdrFile >> which[i] >> count[i] >> where[i];
            }
        }
        ParallelDescriptor::Bcast(which.dataPtr(), which.size(), IOProc);
        ParallelDescriptor::Bcast(count.dataPtr(), count.size(), IOProc);
        ParallelDescriptor::Bcast(where.dataPtr(), where.size(), IOProc);

        m_particles.resize(m_gdb->finestLevel()+1);

	MFInfo info;
	info.SetAlloc(false);
	MultiFab state(m_gdb->ParticleBoxArray(lev),
		       m_gdb->ParticleDistributionMap(lev),
		       1,0,info);

        for (MFIter mfi(state); mfi.isValid(); ++mfi)
        {
            const int grid = mfi.index();

            if (count[grid] <= 0) continue;
            //
            // The file names in the header file are relative.
            //
            std::string name = fullname;

            if (!name.empty() && name[name.size()-1] != '/')
                name += '/';

            name += "Level_";
            name += amrex::Concatenate("", lev, 1);
            name += '/';
            name += ParticleType::DataPrefix();
            name += amrex::Concatenate("", which[grid], 4);

            std::ifstream ParticleFile;

            ParticleFile.open(name.c_str(), std::ios::in);

            if (!ParticleFile.good())
                amrex::FileOpenFailed(name);

            ParticleFile.seekg(where[grid], std::ios::beg);

            if (how == "single")
            {
                ReadParticles_SinglePrecision(count[grid],grid,lev,is_checkpoint,ParticleFile);
            }
            else if (how == "double")
            {
                ReadParticles_DoublePrecision(count[grid],grid,lev,is_checkpoint,ParticleFile);
            }
            else
            {
                std::string msg("ParticleContainer<NR,NI,C>::Restart_Doit(): bad parameter: ");
                msg += how;
                amrex::Error(msg.c_str());
            }
                
            ParticleFile.close();

            if (!ParticleFile.good())
                amrex::Abort("ParticleContainer<NR,NI,C>::Restart_Doit(): problem reading particles");
        }
    }

    BL_ASSERT(OK());        
}

//
// This one stores real data as doubles.
//

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::ReadParticles_DoublePrecision (int            cnt,
							 int            grd,
							 int            lev,
							 bool           is_checkpoint,
							 std::ifstream& ifs)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::ReadParticles_DoublePrecision()");
    BL_ASSERT(cnt > 0);
    BL_ASSERT(lev < int(m_particles.size()));
    BL_ASSERT(lev >= 0 && lev <= m_gdb->finestLevel());
    BL_ASSERT(grd >= 0 && grd < m_gdb->ParticleBoxArray(lev).size());
    //
    // First read in the integer data in binary.  We do not store
    // the m_lev and m_grid data on disk.  We can easily recreate
    // that given the structure of the checkpoint file.
    //
    const int iChunkSize = 2;

    Array<int> istuff(cnt*iChunkSize);

    if (is_checkpoint)
        ifs.read((char*)istuff.dataPtr(),istuff.size()*sizeof(int));
    //
    // Then the double data in binary.
    //
    const int rChunkSize = BL_SPACEDIM+NR;

    Array<double> rstuff(cnt*rChunkSize);

    ifs.read((char*)rstuff.dataPtr(),rstuff.size()*sizeof(double));
    //
    // Now reassemble the particles.
    //
    int*            iptr = istuff.dataPtr();
    double*         rptr = rstuff.dataPtr();
    PBox&           pbox = m_particles[lev][grd];
    const Geometry& geom = m_gdb->Geom(0);

    const Real ProbLo[BL_SPACEDIM] = { D_DECL(geom.ProbLo(0), geom.ProbLo(1), geom.ProbLo(2)) };
    const Real ProbHi[BL_SPACEDIM] = { D_DECL(geom.ProbHi(0), geom.ProbHi(1), geom.ProbHi(2)) };
    const Real  Delta[BL_SPACEDIM] = { D_DECL(Real(.125)*geom.CellSize(0),
                                              Real(.125)*geom.CellSize(1),
                                              Real(.125)*geom.CellSize(2)) };
    ParticleType p;

    // If we are restarting from a plotfile instead of a checkpoint file, then we do not
    //    read in the particle id's, so we need to reset the id counter to zero and renumber them
    if (!is_checkpoint)
    {
        int maxnextid = 1;
        ParticleType::NextID(maxnextid);
    }

    for (int i = 0; i < cnt; i++)
    {
        if (is_checkpoint)
        {
            p.idata.id   = iptr[0];
            p.idata.cpu  = iptr[1];
        }
        else
        {
	  ParticleLocData pld;
	  if (!ParticleType::Where(p, m_gdb, pld))
            {
                ParticleType::PeriodicShift(p, m_gdb);

                if (!ParticleType::Where(p, m_gdb, pld))
                {
		  std::cout << "RESTART:BAD PARTICLE ID WOULD BE " << ParticleType::NextID() << '\n';

		  for (int d = 0; d < BL_SPACEDIM; d++)
                    {
		      std::cout << "RESTART:BAD PARTICLE POS(" << d << ") " << p.rdata.pos[d] << std::endl;
                    }

		  amrex::Abort("ParticleContainer<NR,NI,C>::ReadParticles_DoublePrecision(): invalid particle");
                }
            }
	  
	  p.idata.id   = ParticleType::NextID();
	  p.idata.cpu  = ParallelDescriptor::MyProc();
        }

        BL_ASSERT(p.idata.id > 0);

        iptr += iChunkSize;

        D_TERM(p.rdata.pos[0] = rptr[0];,
               p.rdata.pos[1] = rptr[1];,
               p.rdata.pos[2] = rptr[2];);
        //
        // If we're reading in doubles and storing'm in floats we have
        // to make sure the particles stay in the domain.
        //
        for (int d = 0; d < BL_SPACEDIM; d++)
        {
            if (p.rdata.pos[d] <= ProbLo[d]) p.rdata.pos[d] += Delta[d];
            if (p.rdata.pos[d] >= ProbHi[d]) p.rdata.pos[d] -= Delta[d];
        }

        for (int i = 0; i < NR; i++)
            p.rdata.arr[BL_SPACEDIM+i] = rptr[BL_SPACEDIM+i];

        rptr += rChunkSize;

        pbox.push_back(p);
    }
}

//
// This one stores real data as floats.
//

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::ReadParticles_SinglePrecision (int            cnt,
							   int            grd,
							   int            lev,
							   bool           is_checkpoint,
							   std::ifstream& ifs)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::ReadParticles_SinglePrecision()");
    BL_ASSERT(cnt > 0);
    BL_ASSERT(lev < int(m_particles.size()));
    BL_ASSERT(lev >= 0 && lev <= m_gdb->finestLevel());
    BL_ASSERT(grd >= 0 && grd < m_gdb->ParticleBoxArray(lev).size());
    //
    // First read in the integer data in binary.  We do not store
    // the m_lev and m_grid data on disk.  We can easily recreate
    // that given the structure of the checkpoint file.
    //
    const int iChunkSize = 2;

    Array<int> istuff(cnt*iChunkSize);

    if (is_checkpoint)
        ifs.read((char*)istuff.dataPtr(),istuff.size()*sizeof(int));
    //
    // Then the float data in binary.
    //
    const int rChunkSize = BL_SPACEDIM+NR;

    Array<float> rstuff(cnt*rChunkSize);

    ifs.read((char*)rstuff.dataPtr(),rstuff.size()*sizeof(float));
    //
    // Now reassemble the particles.
    //
    int*   iptr = istuff.dataPtr();
    float* rptr = rstuff.dataPtr();
    PBox&  pbox = m_particles[lev][grd];

    ParticleLocData pld;
    ParticleType p;

    // If we are restarting from a plotfile instead of a checkpoint file, then we do not
    //    read in the particle id's, so we need to reset the id counter to zero and renumber them
    if (!is_checkpoint)
    {
        int maxnextid = 1;
        ParticleType::NextID(maxnextid);
    }

    for (int i = 0; i < cnt; i++)
    {
        p.idata.id   = iptr[0];
        p.idata.cpu  = iptr[1];

        if (is_checkpoint)
        {
            p.idata.id    = iptr[0];
            p.idata.cpu   = iptr[1];
        }
        else
        {
	  if (!ParticleType::Where(p, m_gdb, pld))
            {
	      ParticleType::PeriodicShift(p, m_gdb);

	      if (!ParticleType::Where(p, m_gdb, pld))
                {
                    std::cout << "RESTART:BAD PARTICLE ID WOULD BE " << ParticleType::NextID() << '\n';

                    for (int d = 0; d < BL_SPACEDIM; d++)
                    {
                        std::cout << "RESTART:BAD PARTICLE POS(" << d << ") " << p.rdata.pos[d] << std::endl;
                    }

                    amrex::Abort("ParticleContainer<NR,NI,C>::ReadParticles_SinglePrecision(): invalid particle");
                }
            }

            p.idata.id   = ParticleType::NextID();
            p.idata.cpu  = ParallelDescriptor::MyProc();
        }

        BL_ASSERT(p.idata.id > 0);

        iptr += iChunkSize;

        D_TERM(p.rdata.pos[0] = rptr[0];,
               p.rdata.pos[1] = rptr[1];,
               p.rdata.pos[2] = rptr[2];);

        for (int i = 0; i < NR; i++)
            p.rdata.arr[BL_SPACEDIM + i] = rptr[BL_SPACEDIM + i];

        rptr += rChunkSize;

        pbox.push_back(p);
    }
}

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::WriteAsciiFile (const std::string& filename)
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::WriteAsciiFile()");
    BL_ASSERT(!filename.empty());

    const Real strttime = ParallelDescriptor::second();
    //
    // Count # of valid particles.
    //
    long nparticles = 0;

    for (const auto& pmap : m_particles)
    {
	for (const auto& kv : pmap)
        {
            const PBox& pbox = kv.second;

	    for (const auto& p : pbox)
            {
                if (p.idata.id > 0)
                    //
                    // Only count (and checkpoint) valid particles.
                    //
                    nparticles++;
            }
        }
    }
    //
    // And send count to I/O processor.
    //
    ParallelDescriptor::ReduceLongSum(nparticles,ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Have I/O processor open file and write out particle count.
        //
        std::ofstream File;

        File.open(filename.c_str(), std::ios::out|std::ios::trunc);

        if (!File.good())
            amrex::FileOpenFailed(filename);

        File << nparticles << '\n';
            
        File.flush();

        File.close();

        if (!File.good())
            amrex::Abort("ParticleContainer<NR,NI,C>::WriteAsciiFile(): problem writing file");
    }

    ParallelDescriptor::Barrier();

    const int MyProc = ParallelDescriptor::MyProc();

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (MyProc == i)
        {
            //
            // Each CPU opens the file for appending and adds its particles.
            //
            std::ofstream File;

            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

            File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

            File.open(filename.c_str(), std::ios::out|std::ios::app);

            File.precision(15);

            if (!File.good())
                amrex::FileOpenFailed(filename);
            
            for (const auto& pmap : m_particles)
            {
		for (const auto& kv : pmap)
                {
                    const PBox& pbox = kv.second;

		    for (auto it = pbox.cbegin(); it != pbox.cend(); ++it)
                    {
                        if (it->idata.id > 0)
                        {
                            D_TERM(File << it->rdata.pos[0] << ' ',
                                        << it->rdata.pos[1] << ' ',
                                        << it->rdata.pos[2] << ' ');

                            for (int i = BL_SPACEDIM; i < BL_SPACEDIM + NR; i++)
                            {
                                char ws = (i == BL_SPACEDIM + NR - 1) ? '\n' : ' ';

                                File << it->rdata.arr[i] << ws;
                            }
                        }
                    }
                }
            }

            File.flush();

            File.close();

            if (!File.good())
                amrex::Abort("ParticleContainer<NR,NI,C>::WriteAsciiFile(): problem writing file");

        }

        ParallelDescriptor::Barrier();
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<NR,NI,C>::WriteAsciiFile() time: " << stoptime << '\n';
        }
    }
}

template<int NR, int NI>
inline
void
Particle<NR, NI>::CIC_Fracs (const Real* frac, Real* fracs)
{
    //
    // "frac"  should be dimensioned: Real frac[BL_SPACEDIM]
    //
    // "fracs" should be dimensioned: Real fracs[D_TERM(2,+2,+4)]
    //
#if (BL_SPACEDIM == 1)
    // High
    fracs[0] = frac[0];

    // Low
    fracs[1] = (1-frac[0]);

#elif (BL_SPACEDIM == 2)
    // HH
    fracs[0] = frac[0] * frac[1] ;
    
    // LH
    fracs[1] = (1-frac[0]) * frac[1];
    
    // LL
    fracs[2] = (1-frac[0]) * (1-frac[1]);
    
    // HL
    fracs[3] = frac[0] * (1-frac[1]);

#elif (BL_SPACEDIM == 3)
    // HHH
    fracs[0] = frac[0] * frac[1] * frac[2];

    // LHH
    fracs[1] = (1-frac[0]) * frac[1] * frac[2];

    // LLH
    fracs[2] = (1-frac[0]) * (1-frac[1]) * frac[2];
    
    // HLH
    fracs[3] = frac[0] * (1-frac[1]) * frac[2];

    // HHL
    fracs[4] = frac[0] * frac[1] * (1-frac[2]);
    
    // LHL
    fracs[5] = (1-frac[0]) * frac[1] * (1-frac[2]);

    // LLL
    fracs[6] = (1-frac[0]) * (1-frac[1]) * (1-frac[2]);
    
    // HLL
    fracs[7] = frac[0] * (1-frac[1]) * (1-frac[2]);
#endif
}

template<int NR, int NI>
inline
void
Particle<NR, NI>::CIC_Cells (const IntVect& hicell, IntVect* cells)
{
    //
    // "cells" should be dimensioned: IntVect cells[D_TERM(2,+2,+4)]
    //
    IntVect cell = hicell;

#if (BL_SPACEDIM == 1)
    // High
    cells[0] = cell;

    // Low
    cell[0]  = cell[0] - 1;
    cells[1] = cell;

#elif (BL_SPACEDIM == 2)
    // HH
    cells[0] = cell;
    
    // LH
    cell[0]  = cell[0] - 1;
    cells[1] = cell;
    
    // LL
    cell[1]  = cell[1] - 1;
    cells[2] = cell;
    
    // HL
    cell[0]  = cell[0] + 1;
    cells[3] = cell;

#elif (BL_SPACEDIM == 3)
    // HHH
    cells[0] = cell;

    // LHH
    cell[0]  = cell[0] - 1;
    cells[1] = cell;

    // LLH
    cell[1]  = cell[1] - 1;
    cells[2] = cell;
    
    // HLH
    cell[0]  = cell[0] + 1;
    cells[3] = cell;

    cell = hicell;

    // HHL
    cell[2]  = cell[2] - 1;
    cells[4] = cell;
    
    // LHL
    cell[0]  = cell[0] - 1;
    cells[5] = cell;

    // LLL
    cell[1]  = cell[1] - 1;
    cells[6] = cell;
    
    // HLL
    cell[0]  = cell[0] + 1;
    cells[7] = cell;
#endif
}

template<int NR, int NI>
inline
int
Particle<NR, NI>::CIC_Cells_Fracs (const Particle<NR, NI>& p,
				   const Real*         plo,
				   const Real*         dx,
				   Array<Real>&        fracs,
				   Array<IntVect>&     cells)
{
    return Particle<NR, NI>::CIC_Cells_Fracs(p,plo,dx,dx,fracs,cells);
}

//
// This is the multi-level version.
//
// The Array should be empty on input.
//
// There'll be finest_level+1 of them.
//
template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::AssignDensity (int rho_index, bool sub_cycle,
					   Array<std::unique_ptr<MultiFab> >& mf_to_be_filled, 
					   int lev_min, int ncomp, int finest_level) const
{
    if (rho_index != 0) amrex::Abort("AssignDensity only works if rho_index = 0");

    BL_PROFILE("ParticleContainer<NR,NI,C>::AssignDensity()");
    BL_ASSERT(NR >= 1);
    BL_ASSERT(NR >= ncomp);
    BL_ASSERT(ncomp == 1 || ncomp == BL_SPACEDIM+1);

    if (finest_level == -1)
    {
        finest_level = m_gdb->finestLevel();
    }
    while (!m_gdb->LevelDefined(finest_level))
    {
        finest_level--;
    }
    //
    // The size of the returned multifab is limited by lev_min and 
    // finest_level. In the following code, lev is the real level, 
    // lev_index is the corresponding index for mf. 
    //

    // Create the space for mf_to_be_filled, regardless of whether we'll need a temporary mf
    mf_to_be_filled.resize(finest_level+1-lev_min);
    for (int lev = lev_min; lev <= finest_level; lev++)
    { 
        const int lev_index = lev - lev_min;
        mf_to_be_filled[lev_index].reset(new MultiFab(m_gdb->boxArray(lev),
						      m_gdb->DistributionMap(lev),
						      ncomp, 1));
	mf_to_be_filled[lev_index]->setVal(0.0);
    }

    // Test whether the grid structure of the boxArray is the same
    //       as the ParticleBoxArray at all levels 
    bool all_grids_the_same = true; 
    for (int lev = lev_min; lev <= finest_level; lev++) {
        if (!OnSameGrids(lev, *mf_to_be_filled[lev-lev_min])) {
	    all_grids_the_same = false;
	    break;
	}
    }

    Array<std::unique_ptr<MultiFab> > mf_part;
    if (!all_grids_the_same)
    { 
        // Create the space for the temporary, mf_part
        mf_part.resize(finest_level+1-lev_min);
        for (int lev = lev_min; lev <= finest_level; lev++)
        {
            const int lev_index = lev - lev_min;
            mf_part[lev_index].reset(new MultiFab(m_gdb->ParticleBoxArray(lev), 
						  m_gdb->ParticleDistributionMap(lev),
						  ncomp, 1));
	    mf_part[lev_index]->setVal(0.0);
        }
    }

    auto & mf = (all_grids_the_same) ? mf_to_be_filled : mf_part;

    if (finest_level == 0)
    {
        //
        // Just use the far simpler single-level version.
        //
        AssignDensitySingleLevel(rho_index, *mf[0],0,ncomp);
        //
        // I believe that we don't need any information in ghost cells so we don't copy those.
        //
        if ( ! all_grids_the_same) {
            mf_to_be_filled[0]->copy(*mf[0],0,0,ncomp);
	}
        return;
    }
    
    //
    // This is the "data" needed by other MPI procs.
    //
    std::map<int, std::deque<ParticleCommData> > data;

    const Real stime = ParallelDescriptor::second();
    //
    // Minimum M required.
    //
    const int M = D_TERM(2,+2,+4);

    Array<int>     cgrid(M);
    Array<int>    cwhich(M),  fwhich(M);
    Array<Real>    fracs(M),  cfracs(M);
    Array<IntVect> cells(M),  ccells(M), cfshifts(M);

    ParticleCommData pb;
    //
    // I'm going to allocate these badboys here & pass'm into routines that use'm.
    // This should greatly cut down on memory allocation/deallocation.
    //
    Array<IntVect>                    pshifts(27);
    std::vector< std::pair<int,Box> > isects;
    Array<int>                        fgrid(M);
    Array<Real>                       ffracs(M);
    Array<IntVect>                    fcells;
    //
    // "fvalid" contains all the valid region of the MultiFab at this level, together
    // with any ghost cells lying outside the domain, that can be periodically shifted into the
    // valid region.  "compfvalid" is the complement of the "fvalid", while "compfvalid_grown" is 
    // "compfvalid" grown by one.  Using these we can figure out whether or not a cell is in the
    // valid region of our MultiFab as well as whether or not we're at a Fine->Crse boundary.
    //
    for (int lev = lev_min; lev <= finest_level; lev++)
    {
        const Geometry& gm        = m_gdb->Geom(lev);
        const Geometry& gm_fine   = (lev < finest_level) ? m_gdb->Geom(lev+1) : gm;
        const Geometry& gm_coarse = (lev > 0) ? m_gdb->Geom(lev-1) : gm;
        const Box&      dm        = gm.Domain();
        const Real*     dx        = gm.CellSize();
        const Real*     plo       = gm.ProbLo();
        const Real*     dx_fine   = (lev < finest_level) ? m_gdb->Geom(lev+1).CellSize() : dx;
        const Real*     dx_coarse = (lev > 0) ? m_gdb->Geom(lev-1).CellSize() : dx;
        const int       lev_index = lev - lev_min;
        const BoxArray& grids     = mf[lev_index]->boxArray();
        const int       dgrow     = (lev == 0) ? 1 : m_gdb->MaxRefRatio(lev-1);

        BoxArray compfvalid, compfvalid_grown, fvalid = mf[lev_index]->boxArray();
        //
        // Do we have Fine->Crse overlap on a periodic boundary?
        // We want to add all ghost cells that can be shifted into valid region.
        //
        BoxList valid;

        for (int i = 0; i < grids.size(); i++)
        {
            if (gm.isAnyPeriodic())
            {
                const Box& dest = amrex::grow(grids[i],dgrow);

                if ( ! dm.contains(dest))
                {
                    for (int j = 0; j < grids.size(); j++)
                    {
                        BL_ASSERT(dm.contains(grids[j]));

                        gm.periodicShift(dest, grids[j], pshifts);

			for (const auto& kiv : pshifts)
                        {
                            const Box& sbx = grids[j] + kiv;
                            const Box& dbx = dest & sbx;

                            BL_ASSERT(dbx.ok());

                            valid.push_back(dbx);
                        }
                    }
                }
            }
        }
        if (valid.isNotEmpty())
        {
            //
            // We've got some Fine->Crse periodic overlap.
            // Don't forget to add the valid boxes too.
            //
            for (int i = 0; i < grids.size(); i++) {
                valid.push_back(grids[i]);
	    }
            fvalid = BoxArray(valid);
            fvalid.removeOverlap();
        }
        //
        // If we're at a lev < finestLevel, this is the coarsened fine BoxArray.
        // We use this for figuring out Crse->Fine issues.
        //
        BoxArray ccba;
        if (lev > 0)
        {
            ccba = m_gdb->boxArray(lev);
            ccba.coarsen(m_gdb->refRatio(lev-1));
        }
        BoxArray cfba;
        if (lev < finest_level)
        {
            cfba = m_gdb->boxArray(lev+1);
            cfba.coarsen(m_gdb->refRatio(lev));

            BL_ASSERT(mf[lev_index]->boxArray().contains(cfba));
        }
        //
        // This is cfba with any shifted ghost cells.
        //
        BoxArray cfvalid = cfba;

        if (lev < finest_level)
        {
            BoxList cvalid;

            const BoxArray& cgrids = mf[lev_index]->boxArray();

            for (int i = 0; i < cfba.size(); i++)
            {
                if (gm.isAnyPeriodic())
                {
                    const Box& dest = amrex::grow(cfba[i],mf[lev_index]->nGrow());

                    if ( ! dm.contains(dest))
                   { 
                        for (int j = 0; j < cgrids.size(); j++)
                        {
                            BL_ASSERT(dm.contains(cgrids[j]));

                            gm.periodicShift(dest, cgrids[j], pshifts);

			    for (const auto& kiv : pshifts)
                            {
                                const Box& sbx = cfba[i] - kiv;

                                cvalid.push_back(sbx);
                            }
                        }
                    }
                }
            }
            if (cvalid.isNotEmpty())
            {
                //
                // We've got some Fine->Crse periodic overlap.
                // Don't forget to add the valid boxes too.
                //
                for (int i = 0; i < cfba.size(); i++) {
                    cvalid.push_back(cfba[i]);
		}
                cfvalid = BoxArray(cvalid);
                cfvalid.removeOverlap();
            }
        }
        //
        // The "+1" is so we enclose the valid region together with any
        //  ghost cells that can be periodically shifted into valid.
        //
        compfvalid = amrex::complementIn(amrex::grow(dm,dgrow+1), fvalid);

        compfvalid_grown = compfvalid;
        compfvalid_grown.grow(1);
        compfvalid_grown.removeOverlap();
            
        if (gm.isAnyPeriodic() && ! gm.isAllPeriodic())
        {
            amrex::Error("AssignDensity: problem must be periodic in no or all directions");
        }
        //
        // If we're at a lev > 0, this is the coarsened BoxArray.
        // We use this for figuring out Fine->Crse issues.
        //
        BoxArray cba;
        if (lev > 0)
        {
            cba = m_gdb->boxArray(lev);
            cba.coarsen(m_gdb->refRatio(lev-1));
        }
        //
        // Do the grids at this level cover the full domain? If they do
        // there can be no Fine->Crse interactions at this level.
        //
        const bool GridsCoverDomain = fvalid.contains(m_gdb->Geom(lev).Domain());
        
	for (const auto& kv : m_particles[lev])
        {
            const PBox& pbx = kv.second;
            FArrayBox&  fab = (*mf[lev_index])[kv.first];

	    for (const auto& p : pbx)
            {
                if (p.idata.id <= 0) {
		  continue;
		}
                //
                // Get "fracs" and "cells" for the particle "p" at this level.
                //
                const int M = ParticleType::CIC_Cells_Fracs(p, plo, dx, fracs, cells);
                //
                // If this is not fully periodic then we have to be careful that no
                // particle's support leaves the domain. We test this by checking the low
                // and high corners respectively.
                //
                if ( ! gm.isAllPeriodic() && ! allow_particles_near_boundary) {
                    if ( ! gm.Domain().contains(cells[0]) || ! gm.Domain().contains(cells[M-1])) {
                        amrex::Error("AssignDensity: if not periodic, all particles must stay away from the domain boundary");
		    }
		}
                //
                // This section differs based on whether we subcycle.
                // Without subcycling we use the "stretchy" support for particles.
                // With subcycling a particles support is strictly defined 
                // by its resident level.
                //
                if (sub_cycle)
                {
                    bool isFiner    = false;
                    bool isBoundary = false;
                    //
                    // First sum the mass in the valid region
                    //
                    for (int i = 0; i < M; i++)
                    {
                        if (cfvalid.contains(cells[i]))
                        {
                            //
                            // Some part of the particle's mass lies in a 
                            // finer region; we'll deal with it shortly.
                            //
                            isFiner    = true;
                            isBoundary = true;
                            continue;
                        }
                        if ( ! fvalid.contains(cells[i]))
                        {
                            //
                            // We're out of the valid region.
                            //
                            isBoundary = true;
                            continue;
                        }
                        //
                        // Sum up mass in first component.
                        //
                        {
                            fab(cells[i],0) += p.rdata.arr[BL_SPACEDIM] * fracs[i];
                        }
                        //
                        // Sum up momenta in next components.
                        //

                        // If the domain is not periodic and we want to let particles
                        //    live near the boundary but "throw away" the contribution that 
                        //    does not fall into the domain ...
                        if ( ! gm.isAllPeriodic() && allow_particles_near_boundary &&
			     ! gm.Domain().contains(cells[i]))
			{
			  continue;
			}

                        for (int n = 1; n < ncomp; n++) {
                            fab(cells[i],n) += p.rdata.arr[BL_SPACEDIM + n] * p.rdata.arr[BL_SPACEDIM] * fracs[i];
			}
                    }
                    //
                    // Deal with mass that doesn't belong at this level.
                    // Here we assume proper nesting so that only one special case can
                    // be true for a given particle.
                    //
                    if (isBoundary)
                    {
                        if (isFiner)
                        {
                            BL_ASSERT(lev < finest_level);
                            //
                            // We're at a coarse->fine interface
                            //
                            // get fine cells/fracs
                            //
                            const int MF = ParticleType::CIC_Cells_Fracs(p, plo, dx_fine ,dx, ffracs, fcells);

                            for (int j = 0; j < MF; j++)
                            {
                                //
                                // Make sure this fine cell is valid. Check for periodicity.
                                //
                                const Box bx(fcells[j],fcells[j]);
                                gm_fine.periodicShift(bx, gm_fine.Domain(), pshifts);
                                if ( ! pshifts.empty())
                                {
                                    BL_ASSERT(int(pshifts.size()) == 1);
                                    fcells[j] = fcells[j] - pshifts[0];
                                }
                                mf[lev_index + 1]->boxArray().intersections(Box(fcells[j],fcells[j]),isects,true,0);
                                if (isects.size() == 0) {
                                    continue;
				}
                                const int grid = isects[0].first; 
                                const int who  = mf[lev_index+1]->DistributionMap()[grid];

                                if (who == ParallelDescriptor::MyProc())
                                {
                                    //
                                    // Sum up mass in first component.
                                    //
                                    {
                                        (*mf[lev_index+1])[grid](fcells[j],0) += p.rdata.arr[BL_SPACEDIM] * ffracs[j];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        (*mf[lev_index+1])[grid](fcells[j],n) += p.rdata.arr[BL_SPACEDIM+n] * p.rdata.arr[BL_SPACEDIM] * ffracs[j];
				    }
                                }
                                else
                                {

				  pb.m_lev  = lev+1;
				  pb.m_grid = grid;
				  pb.m_cell = fcells[j];

				  //
				  // Sum up mass in first component.
				  //
				  {
				    pb.m_data[0] = p.rdata.arr[BL_SPACEDIM] *  ffracs[j];
				  }
				  
				  //
				  // Sum up momenta in next components.
				  //
				  for (int n = 1; n < ncomp; n++) {
				    pb.m_data[n] = p.rdata.arr[BL_SPACEDIM+n] * p.rdata.arr[BL_SPACEDIM] * ffracs[j];
				  }
				  
				  data[who].push_back(pb);
                                }
                            }
                        }
                        else if (lev_index > 0)
                        {
                            //
                            // We must be at a fine->coarse interface.
                            //
                            const int MC = ParticleType::CIC_Cells_Fracs(p, plo, dx_coarse, dx, cfracs, ccells);
                            for (int j = 0; j < MC; j++)
                            {
                                //
                                // Make sure this coarse cell isn't in this level's valid region.
                                // This may not matter.
                                //
                                if (cba.contains(ccells[j]))
                                    continue;
                                //
                                // Check for periodicity.
                                //
                                const Box bx(ccells[j],ccells[j]);
                                gm_coarse.periodicShift(bx, gm_coarse.Domain(), pshifts);

                                if ( ! pshifts.empty())
                                {
                                    BL_ASSERT(int(pshifts.size()) == 1);
                                    ccells[j] = ccells[j] - pshifts[0]; 
                                }
                                //
                                // Find its resident grid.
                                //
                                mf[lev_index - 1]->boxArray().intersections(Box(ccells[j],ccells[j]),isects,true,0);
                                if (isects.size() == 0) {
                                    continue;
				}
                                const int grid = isects[0].first;
                                const int who  = mf[lev_index-1]->DistributionMap()[grid];
                                if (who == ParallelDescriptor::MyProc())
                                {
                                    //
                                    // Sum up mass in first component.
                                    //
                                    {
                                        (*mf[lev_index-1])[grid](ccells[j],0) += p.rdata.arr[BL_SPACEDIM] * cfracs[j];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        (*mf[lev_index-1])[grid](ccells[j],n) += p.rdata.arr[BL_SPACEDIM+n] * p.rdata.arr[BL_SPACEDIM] * cfracs[j];
				    }
                                }
                                else
                                {

				  pb.m_lev  = lev-1;
				  pb.m_grid = grid;
				  pb.m_cell = ccells[j];

                                  //
				  // Sum up mass in first component.
				  //
				  {
				    pb.m_data[0] = p.rdata.arr[BL_SPACEDIM] * cfracs[j];
				  }
                                  
				  //
				  // Sum up momenta in next components.
				  //
				  for (int n = 1; n < ncomp; n++) {
				    pb.m_data[n] = p.rdata.arr[BL_SPACEDIM+n] * p.rdata.arr[BL_SPACEDIM] * cfracs[j];
				  }
				  
				  data[who].push_back(pb);
                                }
                            }
                        }
                        else
                        {
                            // The mass is below levels we care about. Ignore it.
                        }
                    }
                }
                else 
                {
                    bool AnyCrseToFine = false;
                    if (lev < finest_level) {
                        AnyCrseToFine = ParticleType::CrseToFine(cfba,cells,cfshifts,gm,cwhich,pshifts);
		    }
                    //
                    // lev_index > 0 means that we don't do F->C for lower levels
                    // This may mean that the mass fraction is off.
                    //
                    bool AnyFineToCrse = false;
                    if (lev_index > 0 && !GridsCoverDomain)
                        AnyFineToCrse = ParticleType::FineToCrse(p,lev,m_gdb,cells,fvalid,compfvalid_grown,ccells,cfracs,fwhich,cgrid,pshifts,isects);

                    BL_ASSERT(!(AnyCrseToFine && AnyFineToCrse));

                    if ( ! AnyCrseToFine && ! AnyFineToCrse)
                    {
                        //
                        // By far the most common case.  Just do it!
                        //
                        for (int i = 0; i < M; i++)
                        {

                            // If the domain is not periodic and we want to let particles
                            //    live near the boundary but "throw away" the contribution that 
                            //    does not fall into the domain ...
                            if (! gm.isAllPeriodic() && allow_particles_near_boundary && ! gm.Domain().contains(cells[i]))
			    {
			      continue;
			    }
                            //
                            // Sum up mass in first component.
                            //
                            {
                                fab(cells[i],0) += p.rdata.arr[BL_SPACEDIM] * fracs[i];
                            }
                            //
                            // Sum up momenta in next components.
                            //
                            for (int n = 1; n < ncomp; n++) {
                                fab(cells[i],n) += p.rdata.arr[BL_SPACEDIM+n] * p.rdata.arr[BL_SPACEDIM] * fracs[i];
			    }
                        }
                    }
                    else if (AnyFineToCrse)
                    {
                        Real sum_crse = 0, sum_fine = 0;

                        for (int i = 0; i < M; i++)
                        {
                            if (fwhich[i])
                            {
                                //
                                // We're at a Fine->Crse boundary.
                                //
                                BL_ASSERT(cgrid[i] >= 0);
                                BL_ASSERT(cgrid[i] < mf[lev_index-1]->size());
                                //
                                // Here we need to update the crse region.  The coarse
                                // region is always going to be updated if we have a
                                // particle in a cell bordering a Fine->Crse boundary.
                                //
                                const int who = mf[lev_index-1]->DistributionMap()[cgrid[i]];

                                if (who == ParallelDescriptor::MyProc())
                                {
                                    if ( ! (*mf[lev_index-1])[cgrid[i]].box().contains(ccells[i])) {
				      continue;
				    }

                                    // If the domain is not periodic and we want to let particles
                                    //    live near the boundary but "throw away" the contribution that 
                                    //    does not fall into the domain ...
                                    if (! gm_coarse.isAllPeriodic() && allow_particles_near_boundary &&
				        ! gm_coarse.Domain().contains(ccells[i]))
				    {
				      continue;
				    }

                                    //
                                    // Sum up mass in first component.
                                    //
                                    {
                                        (*mf[lev_index-1])[cgrid[i]](ccells[i],0) += p.rdata.arr[BL_SPACEDIM] * cfracs[i];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        (*mf[lev_index-1])[cgrid[i]](ccells[i],n) += p.rdata.arr[BL_SPACEDIM+n] * p.rdata.arr[BL_SPACEDIM] * cfracs[i];
				    }
                                }
                                else
                                {
				  pb.m_lev  = lev-1;
				  pb.m_grid = cgrid[i];
				  pb.m_cell = ccells[i];

                                  //
				  // Sum up mass in first component.
				  //
				  {
				    pb.m_data[0] = p.rdata.arr[BL_SPACEDIM] * cfracs[i];
				  }

				  //
				  // Sum up momenta in next components.
				  //
				  for (int n = 1; n < ncomp; n++) {
				    pb.m_data[n] = p.rdata.arr[BL_SPACEDIM+n] * p.rdata.arr[BL_SPACEDIM] * cfracs[i];
				  }
				  data[who].push_back(pb);
                                }

                                sum_crse += cfracs[i];
                            }
                        }
                        //
                        // We've updated the Crse cells.  Now we have to update the fine
                        // cells in such a way that the total amount of mass we move
                        // around is precisely p.rdata.arr[BL_SPACEDIM]. In other words, the fractions
                        // we use at crse and fine have to sum to zero.  In the fine
                        // case, we have to account for the case where one or more of the
                        // cell indices is not in the valid region of the box containing 
                        // the particle.
                        //
                        sum_fine = 0;
                        for (int i = 0; i < M; i++) 
                        {
                            //
                            // Reusing "fwhich" to indicate fine cells that need massaging.
                            //
                            fwhich[i] = true;

                            if ( ! compfvalid_grown.contains(cells[i]))
                            {
                                //
                                // Go ahead and add the full correct amount to these cells.
                                // They can't touch a Fine->Crse boundary.
                                //
                                sum_fine += fracs[i];
                                //
                                // Sum up mass in first component.
                                //
                                {
                                    fab(cells[i],0) += p.rdata.arr[BL_SPACEDIM] * fracs[i];
                                }
                                //
                                // Sum up momenta in next components.
                                //
                                for (int n = 1; n < ncomp; n++) {
                                    fab(cells[i],n) += p.rdata.arr[BL_SPACEDIM+n] * p.rdata.arr[BL_SPACEDIM] * fracs[i];
				}
                                fwhich[i] = false;
                            }
                            else if (compfvalid.contains(cells[i]))
                            {
                                fwhich[i] = false;
                            }
                        }

                        const Real sum_so_far = sum_crse + sum_fine; 

                        BL_ASSERT(sum_so_far > 0);
                        BL_ASSERT(sum_so_far < 1);

                        sum_fine = 0;
                        for (int i = 0; i < M; i++) 
                        {       
                            if (fwhich[i])
                                //
                                // Got to weight cells in this direction differently.
                                //
                                sum_fine += fracs[i];
                        }

                        const Real mult = (1 - sum_so_far) / sum_fine;
                        //
                        // Now add the weighted amount to the fine cells touching the c-f interface.
                        //
                        sum_fine = 0;
                        for (int i = 0; i < M; i++)
                        {
                            if (fwhich[i])
                            {
                                //
                                // Sum up mass in first component.
                                //
                                {
                                    fab(cells[i],0) += p.rdata.arr[BL_SPACEDIM] * fracs[i] * mult;
                                }
                                //
                                // Sum up momenta in next components.
                                //
                                for (int n = 1; n < ncomp; n++) {
                                    fab(cells[i],n) += p.rdata.arr[BL_SPACEDIM+n] * p.rdata.arr[BL_SPACEDIM] * fracs[i] * mult;
				}

                                sum_fine += fracs[i] * mult;
                            }
                        }

                        BL_ASSERT(std::abs(1-(sum_fine+sum_so_far)) < 1.e-9);
                    }
                    else if (AnyCrseToFine)
                    {
                        Real sum = 0;

                        for (int i = 0; i < M; i++)
                        {
                            if (!cwhich[i])
                            {
                                // If the domain is not periodic and we want to let particles
                                //    live near the boundary but "throw away" the contribution that 
                                //    does not fall into the domain ...
                                if ( ! gm.isAllPeriodic() && allow_particles_near_boundary &&
				     ! gm.Domain().contains(ccells[i]))
				{
				  continue;
				}
                                //
                                // Sum up mass in first component.
                                //
                                {
                                    fab(cells[i],0) += p.rdata.arr[BL_SPACEDIM] * fracs[i];
                                }
                                //
                                // Sum up momenta in next components.
                                //
                                for (int n = 1; n < ncomp; n++) {
                                    fab(cells[i],n) += p.rdata.arr[BL_SPACEDIM+n] * p.rdata.arr[BL_SPACEDIM] * fracs[i];
				}

                                sum += fracs[i];
                            }
                            else
                            {
                                //
                                // We're at a Crse->Fine boundary.
                                //
                                ParticleType::FineCellsToUpdateFromCrse(p,lev,m_gdb,cells[i],cfshifts[i],fgrid,ffracs,fcells,isects);

                                for (int j = 0, nfcells = fcells.size(); j < nfcells; j++)
                                {
                                    const int who = mf[lev_index+1]->DistributionMap()[fgrid[j]];

                                    if (who == ParallelDescriptor::MyProc())
                                    {
                                        //
                                        // Sum up mass in first component.
                                        //
                                        {
                                            (*mf[lev_index+1])[fgrid[j]](fcells[j],0) += p.rdata.arr[BL_SPACEDIM] * fracs[i] * ffracs[j];
                                        }
                                        //
                                        // Sum up momenta in next components.
                                        //
                                        for (int n = 1; n < ncomp; n++) {
                                            (*mf[lev_index+1])[fgrid[j]](fcells[j],n) += p.rdata.arr[BL_SPACEDIM+n] * p.rdata.arr[BL_SPACEDIM] * fracs[i] * ffracs[j];
					}
                                    }
                                    else
                                    {
				      pb.m_lev  = lev+1;
				      pb.m_grid = fgrid[j];
				      pb.m_cell = fcells[j];

                                      //
				      // Sum up mass in first component.
				      //
				      {
					pb.m_data[0] = p.rdata.arr[BL_SPACEDIM] * fracs[i] * ffracs[j];
				      }

				      //
				      // Sum up momenta in next components.
				      //
				      for (int n = 1; n < ncomp; n++) {
					pb.m_data[0] = p.rdata.arr[BL_SPACEDIM+n] * p.rdata.arr[BL_SPACEDIM] * fracs[i] * ffracs[j];
					}

                                        data[who].push_back(pb);
                                    }

                                    sum += fracs[i] * ffracs[j];
                                }
                            }
                        }

                        BL_ASSERT(std::abs(1-sum) < 1.e-9);
                    }
                }
            }
        }
    }

    //
    // Send any needed data to other MPI processes.
    // This "may" touch ghost cells so we want to do it before
    // the SumBoundary() stuff.
    //
    AssignDensityDoit(rho_index, mf, data, ncomp, lev_min);

    for (int lev = lev_min; lev <= finest_level; lev++)
    {
        const int       lev_index = lev - lev_min;
        const Geometry& gm        = m_gdb->Geom(lev);
        const Real*     dx        = gm.CellSize();
        const Real      vol       = D_TERM(dx[0], *dx[1], *dx[2]);

        mf[lev_index]->SumBoundary(gm.periodicity());
        //
        // If ncomp > 1, first divide the momenta (component n) 
        // by the mass (component 0) in order to get velocities.
        // Be careful not to divide by zero.
        //
        for (int n = 1; n < ncomp; n++)
        {
            for (MFIter mfi(*mf[lev_index]); mfi.isValid(); ++mfi)
            {
                (*mf[lev_index])[mfi].protected_divide((*mf[lev_index])[mfi],0,n,1);
            }
        }
        //
        // Only multiply the first component by (1/vol) because this converts mass
        // to density. If there are additional components (like velocity), we don't
        // want to divide those by volume.
        //
        mf[lev_index]->mult(1/vol,0,1);
    }

    //
    // The size of the returned multifab is limited by lev_min and 
    // finest_level. In the following code, lev is the real level,  
    // lev_index is the corresponding index for mf. 
    //
    // I believe that we don't need any information in ghost cells so we don't copy those.
    //
    if ( ! all_grids_the_same)
        for (int lev = lev_min; lev <= finest_level; lev++)
        {
            const int lev_index = lev - lev_min;
            mf_to_be_filled[lev_index]->copy(*mf_part[lev_index],0,0,1);
        }
    
    if (m_verbose > 1)
    {
        Real etime = ParallelDescriptor::second() - stime;

        ParallelDescriptor::ReduceRealMax(etime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<NR,NI,C>::AssignDensity(multi-level) time: " << etime << '\n';
        }
    }
}

//
// Used by AssignDensity (Array<std::unique_ptr<MultiFab> >& mf).
//
// Passes data needed by Crse->Fine or Fine->Crse to CPU that needs it.
//
// We store the data that needs to be sent in "data". Note that m_lev is the
// real particle level, while mf may start at a fine level (e.g. lvls 1 and 2).
// Consequently, we must subtract lev_min from m_lev to get the mf lev.
//

template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::AssignDensityDoit (int               rho_index,
					       Array<std::unique_ptr<MultiFab> >&             mf,
					       std::map<int, std::deque<ParticleCommData> >& data,
					       int               ncomp,
					       int               lev_min) const
{
    if (rho_index != 0) amrex::Abort("AssignDensityDoit only works if rho_index = 0");

    BL_PROFILE("ParticleContainer<NR,NI,C>::AssignDensityDoit()");
    BL_ASSERT(NR >= ncomp);

    const int NProcs = ParallelDescriptor::NProcs();

    if (NProcs == 1)
    {
      BL_ASSERT(data.empty());
      return;
    }

#if BL_USE_MPI
    //
    // We may have data that needs to be sent to another CPU.
    //
    const int MyProc = ParallelDescriptor::MyProc();

    Array<int> Snds(NProcs,0), Rcvs(NProcs,0);

    int NumSnds = 0, NumRcvs = 0;

    for (const auto& kv : data)
    {
        NumSnds       += kv.second.size();
        Snds[kv.first] = kv.second.size();
    }

    ParallelDescriptor::ReduceIntMax(NumSnds);

    if (NumSnds == 0) {
        //
        // There's no parallel work to do.
        //
        return;
    }

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());

    BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 Rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 ParallelDescriptor::Communicator()) );
    BL_ASSERT(Rcvs[MyProc] == 0);

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    typedef std::map<int,int> IntIntMap;

    IntIntMap SndCnts, RcvCnts, rOffset;

    for (int i = 0; i < NProcs; i++) {
        if (Snds[i] > 0) {
            SndCnts[i] = Snds[i];
	}
    }

    for (int i = 0; i < NProcs; i++)
    {
        if (Rcvs[i] > 0)
        {
            RcvCnts[i] = Rcvs[i];
            rOffset[i] = NumRcvs;
            NumRcvs   += Rcvs[i];
        }
    }
    //
    // Don't need these anymore.
    //
    Array<int>().swap(Snds);
    Array<int>().swap(Rcvs);
    //
    // The data we want to receive.
    //
    const int iChunkSize = 2 + BL_SPACEDIM;
    const int rChunkSize = ncomp;

    Array<int>                    irecvdata (NumRcvs*iChunkSize);
    Array<typename ParticleType::RealType> rrecvdata (NumRcvs*rChunkSize);

    Array<int>         index(2*RcvCnts.size());
    Array<MPI_Status>  stats(2*RcvCnts.size());
    Array<MPI_Request> rreqs(2*RcvCnts.size());

    const int SeqNumI = ParallelDescriptor::SeqNum();
    const int SeqNumR = ParallelDescriptor::SeqNum();
    //
    // Post the receives.
    //
    int idx = 0;
    for (auto it = RcvCnts.cbegin(); it != RcvCnts.cend(); ++it, ++idx)
    {
        const int Who  = it->first;
        const int iCnt = it->second   * iChunkSize;
        const int rCnt = it->second   * rChunkSize;
        const int iIdx = rOffset[Who] * iChunkSize;
        const int rIdx = rOffset[Who] * rChunkSize;

        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(iCnt > 0);
        BL_ASSERT(rCnt > 0);
        BL_ASSERT(iCnt < std::numeric_limits<int>::max());
        BL_ASSERT(rCnt < std::numeric_limits<int>::max());

        rreqs[2*idx+0] = ParallelDescriptor::Arecv(&irecvdata[iIdx],iCnt,Who,SeqNumI).req();
        rreqs[2*idx+1] = ParallelDescriptor::Arecv(&rrecvdata[rIdx],rCnt,Who,SeqNumR).req();
    }
    //
    // Send the data.
    //
    Array<int>                             isenddata;
    Array<typename ParticleType::RealType> rsenddata;

    for (const auto& kv : SndCnts)
    {
        const int Who  = kv.first;
        const int iCnt = kv.second * iChunkSize;
        const int rCnt = kv.second * rChunkSize;

        BL_ASSERT(iCnt > 0);
        BL_ASSERT(rCnt > 0);
        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(iCnt < std::numeric_limits<int>::max());
        BL_ASSERT(rCnt < std::numeric_limits<int>::max());

        isenddata.resize(iCnt);
        rsenddata.resize(rCnt);

	std::deque<ParticleCommData>& pbox = data[Who];

        int ioff = 0, roff = 0;
	for (const auto& p : pbox)
        {
	  isenddata[ioff+0] = p.m_lev  - lev_min;
	  isenddata[ioff+1] = p.m_grid;

	  D_TERM(isenddata[ioff+2] = p.m_cell[0];,
		 isenddata[ioff+3] = p.m_cell[1];,
		 isenddata[ioff+4] = p.m_cell[2];);

	  ioff += iChunkSize;

	  for (int n = 0; n < ncomp; n++) {
	    rsenddata[roff+n] = p.m_data[n];
	  }

	  roff += ncomp;
        }

	std::deque<ParticleCommData>().swap(pbox);

        ParallelDescriptor::Send(isenddata.dataPtr(),iCnt,Who,SeqNumI);
        ParallelDescriptor::Send(rsenddata.dataPtr(),rCnt,Who,SeqNumR);
    }
    //
    // Receive the data.
    //
    for (int NWaits = rreqs.size(), completed; NWaits > 0; NWaits -= completed)
    {
        ParallelDescriptor::Waitsome(rreqs, completed, index, stats);
    }
    //
    // Now update "mf".
    //
    if (NumRcvs > 0)
    {
        const int*                             idata = irecvdata.dataPtr();
        const typename ParticleType::RealType* rdata = rrecvdata.dataPtr();

        for (int i = 0; i < NumRcvs; i++)
        {
            const int     lev  = idata[0];
            const int     grd  = idata[1];
            const IntVect cell (D_DECL(idata[2],idata[3],idata[4]));

            BL_ASSERT((*mf[lev]).DistributionMap()[grd] == MyProc);
	    BL_ASSERT((*mf[lev])[grd].box().contains(cell));

            for (int n = 0; n < ncomp; n++) {
                (*mf[lev])[grd](cell,n) += rdata[n];
	    }

            idata += iChunkSize;
            rdata += rChunkSize;
        }
    }

#endif /*BL_USE_MPI*/
}

//
// This is the single-level version -- it takes either cell-centered or node-centered MF's
//
template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::AssignDensitySingleLevel (int rho_index,
						    MultiFab& mf_to_be_filled,
						    int       lev,
						    int       ncomp,
						    int       particle_lvl_offset) const
{
    BL_PROFILE("ParticleContainer<NR,NI,C>::AssignDensitySingleLevel()");
    BL_ASSERT(NR >= 1);
    BL_ASSERT(ncomp == 1 || ncomp == BL_SPACEDIM+1);

    if (lev >= int(m_particles.size()))
    {
        //
        // Don't do anything if there are no particles at this level.
        //
        return;
    }

    // Keep the same external interface to the applications, but if the
    if (mf_to_be_filled.is_nodal())
    {
        NodalDepositionSingleLevel(rho_index, mf_to_be_filled,lev,ncomp,particle_lvl_offset);
    }
    else if (mf_to_be_filled.boxArray().ixType().cellCentered())
    {
        AssignCellDensitySingleLevel(rho_index, mf_to_be_filled,lev,ncomp,particle_lvl_offset);
    }
    else
    {
	amrex::Abort("AssignCellDensitySingleLevel: mixed type not supported");
    }
}

//
// This is the single-level version for cell-centered density
//
template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::AssignCellDensitySingleLevel (int rho_index,
							MultiFab& mf_to_be_filled,
							int       lev,
							int       ncomp,
							int       particle_lvl_offset) const
{
    if (rho_index != 0) amrex::Abort("AssignCellDensitySingleLevel only works if rho_index = 0");

    MultiFab* mf_pointer;

    if (OnSameGrids(lev, mf_to_be_filled))
    {
        // If we are already working with the internal mf defined on the 
        // particle_box_array, then we just work with this.
        mf_pointer = &mf_to_be_filled;
    }
    else
    {
        // If mf_to_be_filled is not defined on the particle_box_array, then we need 
        // to make a temporary here and copy into mf_to_be_filled at the end.
        mf_pointer = new MultiFab(m_gdb->ParticleBoxArray(lev), 
				  m_gdb->ParticleDistributionMap(lev),
				  ncomp, mf_to_be_filled.nGrow());
    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread its effect to an
    //    adjacent grid by first putting the value into ghost cells of its own grid.  The mf->sumBoundary call then
    //    adds the value from one grid's ghost cell to another grid's valid region.
    if (mf_pointer->nGrow() < 1) 
       amrex::Error("Must have at least one ghost cell when in AssignDensitySingleLevel");

    const Real      strttime    = ParallelDescriptor::second();
    const Geometry& gm          = m_gdb->Geom(lev);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx_particle = m_gdb->Geom(lev + particle_lvl_offset).CellSize();
    const Real*     dx          = gm.CellSize();
    const PMap&     pmap        = m_particles[lev];
    const int       ngrids      = pmap.size();

    if (gm.isAnyPeriodic() && ! gm.isAllPeriodic()) {
        amrex::Error("AssignDensity: problem must be periodic in no or all directions");
    }

    for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi) {
        (*mf_pointer)[mfi].setVal(0);
    }
    //
    // This is a little funky.  What in effect this'll do is force
    // each thread to work on a single (separate) grid at a time.  That
    // way no thread will step on any other.  If there's only one grid per CPU,
    // then oh well ....
    //
    // TODO: implement tiling with OpenMP in this grid loop.
    Array<int>         pgrd(ngrids);
    Array<const PBox*> pbxs(ngrids);

    int j = 0;
    for (const auto& kv : pmap)
    {
        pgrd[j] =   kv.first;
        pbxs[j] = &(kv.second);
	++j;
    }

    for (int j = 0; j < ngrids; j++)
    {
        const PBox& pbx = *pbxs[j];
        FArrayBox&  fab = (*mf_pointer)[pgrd[j]];
	auto N = pbx.size();

        Array<Real>    fracs;
        Array<IntVect> cells;

#ifdef _OPENMP
#pragma omp parallel for default(none) private(fracs,cells) shared(N,plo,dx,dx_particle,gm,fab,ncomp,pbx)
#endif
	for (size_t ip = 0; ip < N; ++ip)
        {
            const ParticleType& p = pbx[ip];

            if (p.idata.id <= 0) {
	      continue;
	    }

            const int M = ParticleType::CIC_Cells_Fracs(p, plo, dx, dx_particle, fracs, cells);
            //
            // If this is not fully periodic then we have to be careful that the
            // particle's support leaves the domain unless we specifically want to ignore
            // any contribution outside the boundary (i.e. if allow_particles_near_boundary = true). 
            // We test this by checking the low and high corners respectively.
            //
            if ( ! gm.isAllPeriodic() && ! allow_particles_near_boundary) {
                if ( ! gm.Domain().contains(cells[0]) || ! gm.Domain().contains(cells[M-1])) {
                    amrex::Error("AssignDensity: if not periodic, all particles must stay away from the domain boundary");
		}
	    }

            for (int i = 0; i < M; i++)
            {
                if ( ! fab.box().contains(cells[i])) {
		  continue;
		}

                // If the domain is not periodic and we want to let particles
                //    live near the boundary but "throw away" the contribution that 
                //    does not fall into the domain ...
                if ( ! gm.isAllPeriodic() && allow_particles_near_boundary && ! gm.Domain().contains(cells[i])) {
		  continue;
		}
                //
                // Sum up mass in first component.
                //
                {
#ifdef _OPENMP
#pragma omp atomic
#endif
                    fab(cells[i],0) += p.rdata.arr[BL_SPACEDIM] * fracs[i];
                }
                // 
                // Sum up momenta in next components.
                //
                for (int n = 1; n < ncomp; n++)
#ifdef _OPENMP
#pragma omp atomic
#endif
		  fab(cells[i],n) += p.rdata.arr[BL_SPACEDIM+n] * p.rdata.arr[BL_SPACEDIM] * fracs[i];
            }
        }
    }

    mf_pointer->SumBoundary(gm.periodicity());
    //
    // If ncomp > 1, first divide the momenta (component n) 
    // by the mass (component 0) in order to get velocities.
    // Be careful not to divide by zero.
    //
    for (int n = 1; n < ncomp; n++)
    {
        for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi)
        {
            (*mf_pointer)[mfi].protected_divide((*mf_pointer)[mfi],0,n,1);
        }
    }
    //
    // Only multiply the first component by (1/vol) because this converts mass
    // to density. If there are additional components (like velocity), we don't
    // want to divide those by volume.
    //
    const Real vol = D_TERM(dx[0], *dx[1], *dx[2]);

    mf_pointer->mult(1/vol,0,1);

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled.   I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled)
    {
        mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
	delete mf_pointer;
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<NR,NI,C>::AssignDensity(single-level) time: " << stoptime << '\n';
        }
    }
}

//
// This is the single-level version for nodal density
//
template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::NodalDepositionSingleLevel (int rho_index,
						      MultiFab& mf_to_be_filled,
					              int       lev,
           					      int       ncomp,
	           				      int       particle_lvl_offset) const
{
    MultiFab* mf_pointer;

    if (OnSameGrids(lev, mf_to_be_filled))
    {
        // If we are already working with the internal mf defined on the 
        // particle_box_array, then we just work with this.
        mf_pointer = &mf_to_be_filled;
    }
    else
    {
        // If mf_to_be_filled is not defined on the particle_box_array, then we need 
        // to make a temporary here and copy into mf_to_be_filled at the end.
        mf_pointer = new MultiFab(amrex::convert(m_gdb->ParticleBoxArray(lev),
						  mf_to_be_filled.boxArray().ixType()),
				  m_gdb->ParticleDistributionMap(lev),
				  ncomp, mf_to_be_filled.nGrow());
    }

    const Real      strttime    = ParallelDescriptor::second();
    const Geometry& gm          = m_gdb->Geom(lev);
    const Real*     dx          = gm.CellSize();
    const PMap&     pmap        = m_particles[lev];
    const int       ngrids      = pmap.size();

    if (gm.isAnyPeriodic() && ! gm.isAllPeriodic()) 
        amrex::Error("AssignDensity: problem must be periodic in no or all directions");

    mf_pointer->setVal(0.0);

    //
    // This is a little funky.  What in effect this'll do is force
    // each thread to work on a single (separate) grid at a time.  That
    // way no thread will step on any other.  If there's only one grid per CPU,
    // then oh well ....
    //
    // TODO: implement tiling with OpenMP in this grid loop.
    Array<int>         pgrd(ngrids);
    Array<const PBox*> pbxs(ngrids);

    int j = 0;
    for (const auto& kv : pmap)
    {
        pgrd[j] =   kv.first;
        pbxs[j] = &(kv.second);
	++j;
    }

    Array<IntVect> cells;
    cells.resize(8);

    Array<Real> fracs;
    fracs.resize(8);

#if (BL_SPACEDIM > 1)
    Array<Real> sx;
    sx.resize(2);
    Array<Real> sy;
    sy.resize(2);
#endif
#if (BL_SPACEDIM > 2)
    Array<Real> sz;
    sz.resize(2);
#endif

    for (int j = 0; j < ngrids; j++)
    {
        const PBox& pbx = *pbxs[j];
        FArrayBox&  fab = (*mf_pointer)[pgrd[j]];

	for (const auto& p : pbx)
        {
            if (p.idata.id <= 0) {
	      continue;
	    }
	    
	    ParticleLocData pld;
	    Particle<NR, NI>::Where(p, m_gdb, pld);

#if (BL_SPACEDIM == 1)
            cells[0] = pld.m_cell;
            cells[1] = pld.m_cell+IntVect(1);

            Real x = p.rdata.pos[0] / dx[0];

            int i = pld.m_cell[0];

            Real xint = x - i;

            for (int i = 0; i < 2; i++)
            {
               fab(cells[0],0) += p.rdata.arr[BL_SPACEDIM+rho_index] * (1.0 - xint);
               fab(cells[1],0) += p.rdata.arr[BL_SPACEDIM+rho_index] *        xint ;
            }
#elif (BL_SPACEDIM == 2)
            cells[0] = pld.m_cell;
            cells[1] = pld.m_cell+IntVect(1,0);
            cells[2] = pld.m_cell+IntVect(0,1);
            cells[3] = pld.m_cell+IntVect(1,1);

            Real x = p.rdata.pos[0] / dx[0];
            Real y = p.rdata.pos[1] / dx[1];

            int i = pld.m_cell[0];
            int j = pld.m_cell[1];

            Real xint = x - i;
            Real yint = y - j;

            sx[0] = 1.0-xint;
            sx[1] = xint;
            sy[0] = 1.0-yint;
            sy[1] = yint;

            fracs[0] = sx[0] * sy[0];
            fracs[1] = sx[1] * sy[0];
            fracs[2] = sx[0] * sy[1];
            fracs[3] = sx[1] * sy[1];

            for (int i = 0; i < 4; i++)
            {
               fab(cells[i],0) += p.rdata.arr[BL_SPACEDIM+rho_index] * fracs[i];
            }
#else
            cells[0] = pld.m_cell;
            cells[1] = pld.m_cell+IntVect(1,0,0);
            cells[2] = pld.m_cell+IntVect(0,1,0);
            cells[3] = pld.m_cell+IntVect(1,1,0);
            cells[4] = pld.m_cell+IntVect(0,0,1);
            cells[5] = pld.m_cell+IntVect(1,0,1);
            cells[6] = pld.m_cell+IntVect(0,1,1);
            cells[7] = pld.m_cell+IntVect(1,1,1);

            Real x = p.rdata.pos[0] / dx[0];
            Real y = p.rdata.pos[1] / dx[1];
            Real z = p.rdata.pos[2] / dx[2];

            int i = pld.m_cell[0];
            int j = pld.m_cell[1];
            int k = pld.m_cell[2];

            Real xint = x - i;
            Real yint = y - j;
            Real zint = z - k;

            sx[0] = 1.0-xint;
            sx[1] = xint;
            sy[0] = 1.0-yint;
            sy[1] = yint;
            sz[0] = 1.0-zint;
            sz[1] = zint;

            fracs[0] = sx[0] * sy[0] * sz[0];
            fracs[1] = sx[1] * sy[0] * sz[0];
            fracs[2] = sx[0] * sy[1] * sz[0];
            fracs[3] = sx[1] * sy[1] * sz[0];
            fracs[4] = sx[0] * sy[0] * sz[1];
            fracs[5] = sx[1] * sy[0] * sz[1];
            fracs[6] = sx[0] * sy[1] * sz[1];
            fracs[7] = sx[1] * sy[1] * sz[1];

            for (int i = 0; i < 8; i++)
            {
               fab(cells[i],0) += p.rdata.arr[BL_SPACEDIM+rho_index] * fracs[i];
            }
#endif
        }
    }

    mf_pointer->SumBoundary(gm.periodicity());

    //
    // Only multiply the first component by (1/vol) because this converts mass
    // to density. If there are additional components (like velocity), we don't
    // want to divide those by volume.
    //
    const Real vol = D_TERM(dx[0], *dx[1], *dx[2]);

    mf_pointer->mult(1/vol,0,1);

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled.   I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled)
    {
        mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
	delete mf_pointer;
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::NodalDepositionSingleLevel time: " << stoptime << '\n';
        }
    }
}

//
// This version takes as input the acceleration vector at cell centers, and has the option of
// returning the acceleration at the particle location in the data array, starting at
// component start_comp_for_accel
//
template <int NR, int NI, class C>
void
ParticleContainer<NR,NI,C>::moveKick (MultiFab&       acceleration,
				      int             lev,
				      Real            dt,
				      Real            a_new,
				      Real            a_half, 
				      int             start_comp_for_accel)
{
    BL_PROFILE("ParticleContainer::moveKick()");
    BL_ASSERT(NR >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < int(m_particles.size()));

    const Real strttime  = ParallelDescriptor::second();
    const Real half_dt   = Real(0.5) * dt;
    const Real a_new_inv = 1 / a_new;
    PMap&      pmap      = m_particles[lev];

    MultiFab* ac_pointer;
    if (OnSameGrids(lev,acceleration))
    {
        ac_pointer = &acceleration;
    }
    else 
    {
        ac_pointer = new MultiFab(m_gdb->ParticleBoxArray(lev),
				  m_gdb->ParticleDistributionMap(lev),
				  acceleration.nComp(),acceleration.nGrow());
        for (MFIter mfi(*ac_pointer); mfi.isValid(); ++mfi)
            ac_pointer->setVal(0.);
        ac_pointer->copy(acceleration,0,0,acceleration.nComp());
        ac_pointer->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
    }

    for (auto& kv : pmap)
    {
        const int        grid = kv.first;
        PBox&            pbox = kv.second;
        const int        n    = pbox.size();
        const FArrayBox& gfab = (*ac_pointer)[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.idata.id > 0)
            {
                //
                // Note: rdata.arr[BL_SPACEDIM] is mass, BL_SPACEDIM+1 is v_x, ...
                //
                Real grav[BL_SPACEDIM];

                ParticleType::GetGravity(gfab, m_gdb->Geom(lev), p, grav);
                //
                // Define (a u)^new = (a u)^half + dt/2 grav^new
                //
                D_TERM(p.rdata.arr[BL_SPACEDIM+1] *= a_half;,
                       p.rdata.arr[BL_SPACEDIM+2] *= a_half;,
                       p.rdata.arr[BL_SPACEDIM+3] *= a_half;);

                D_TERM(p.rdata.arr[BL_SPACEDIM+1] += half_dt * grav[0];,
                       p.rdata.arr[BL_SPACEDIM+2] += half_dt * grav[1];,
                       p.rdata.arr[BL_SPACEDIM+3] += half_dt * grav[2];);

                D_TERM(p.rdata.arr[BL_SPACEDIM+1] *= a_new_inv;,
                       p.rdata.arr[BL_SPACEDIM+2] *= a_new_inv;,
                       p.rdata.arr[BL_SPACEDIM+3] *= a_new_inv;);

                if (start_comp_for_accel > BL_SPACEDIM)
                {
                   D_TERM(p.rdata.arr[BL_SPACEDIM + start_comp_for_accel  ] = grav[0];,
                          p.rdata.arr[BL_SPACEDIM + start_comp_for_accel+1] = grav[1];,
                          p.rdata.arr[BL_SPACEDIM + start_comp_for_accel+2] = grav[2];);
                }
            }
        }
    }


    if (ac_pointer != &acceleration) delete ac_pointer;

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<NR,NI,C>::moveKick() time: " << stoptime << '\n';
        }
    }
    //
    // No need for Redistribution(), we only change the velocity.
    //
}

}

#endif /*_PARTICLES_H_*/
