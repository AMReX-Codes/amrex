#ifndef AMREX_PARTICLE_H_
#define AMREX_PARTICLE_H_

#include <string>

#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IntVect.H>

namespace amrex {

namespace
{
    std::string   aggregation_type   = "";
    int           aggregation_buffer = 1;
    constexpr int GhostParticleID    = std::numeric_limits<int>::max();
    constexpr int VirtualParticleID  = std::numeric_limits<int>::max()-1;
    constexpr int LastParticleID     = std::numeric_limits<int>::max()-2;
    constexpr int DoSplitParticleID  = std::numeric_limits<int>::max()-3;
    constexpr int NoSplitParticleID  = std::numeric_limits<int>::max()-4;
}

//
//! The struct used to store particles.
template<int NReal, int NInt=0>
struct Particle
{
    //
    //! The floating point type used for the particles.
#ifdef BL_SINGLE_PRECISION_PARTICLES
    typedef float RealType;
#else
    typedef double RealType;
#endif

    /**
    * The real data. We always have SPACEDIM position coordinates,
    * and optionally we have NReal additional real attributes.
    */
    union rm_t
    {
      RealType pos[AMREX_SPACEDIM];
      RealType arr[AMREX_SPACEDIM+NReal];
    };

    /**
    * The integer data. We always have id and cpu, and optionally we
    * have NInt additional integer attributes.
    */
    union im_t
    {
      struct {int id; int cpu;};
      int arr[2+NInt];
    };

    rm_t m_rdata;
    im_t m_idata;

    static int the_next_id;

    AMREX_GPU_HOST_DEVICE int&  id()       & {return m_idata.id;}
    AMREX_GPU_HOST_DEVICE int   id() const & {return m_idata.id;}
    AMREX_GPU_HOST_DEVICE int& cpu()       & {return m_idata.cpu;}
    AMREX_GPU_HOST_DEVICE int  cpu() const & {return m_idata.cpu;}

    AMREX_GPU_HOST_DEVICE RealType& pos(int index)       & {return m_rdata.pos[index];}
    AMREX_GPU_HOST_DEVICE RealType  pos(int index) const & {return m_rdata.pos[index];}

    AMREX_GPU_HOST_DEVICE RealType& rdata(int index)       & {return m_rdata.arr[AMREX_SPACEDIM + index];}
    AMREX_GPU_HOST_DEVICE RealType  rdata(int index) const & {return m_rdata.arr[AMREX_SPACEDIM + index];}

    AMREX_GPU_HOST_DEVICE int& idata(int index)       & {return m_idata.arr[2 + index];}
    AMREX_GPU_HOST_DEVICE int  idata(int index) const & {return m_idata.arr[2 + index];}

    static Real InterpDoit (const FArrayBox& fab, const Real* fracs, const IntVect* cells, int comp);

    static Real InterpDoit (const FArrayBox& fab, const IntVect& hi, const Real* frac, int comp);

    static void Interp (const Particle<NReal, NInt>& prt, const Geometry& geom, const FArrayBox& fab, const int* idx, Real* val, int cnt);

    static const std::string& Version ();

    static const std::string& DataPrefix ();

    static void GetGravity (const FArrayBox& gfab, const Geometry& geom, const Particle<NReal, NInt>& p, Real* grav);

    static int MaxReaders ();

    static long MaxParticlesPerRead ();


    /**
    * \brief Returns the next particle ID for this processor.
    * Particle IDs start at 1 and are never reused.
    * The pair, consisting of the ID and the CPU on which the particle is "born",
    * is a globally unique identifier for a particle.  The maximum of this value
    * across all processors must be checkpointed and then restored on restart
    * so that we don't reuse particle IDs.
    */
    static int NextID ();


    /**
    * \brief This version can only be used inside omp critical.
    */
    static int UnprotectedNextID ();


    /**
    * \brief Reset on restart.
    *
    * \param nextid
    */
    static void NextID (int nextid);

    static void CIC_Fracs (const Real* frac, Real* fracs);

    static void CIC_Cells (const IntVect& hicell, IntVect* cells);


    /**
    * \brief Old, *-based CIC for use in Interp.
    *
    * \param p
    * \param plo
    * \param dx
    * \param fracs
    * \param cells
    */
    static void CIC_Cells_Fracs_Basic (const Particle<NReal, NInt>& p,
                                       const Real*                  plo,
                                       const Real*                  dx,
                                       Real*                        fracs,
                                       IntVect*                     cells);


    /**
    * \brief Wraps the arbitrary dx function.
    *
    * \param p
    * \param plo
    * \param dx
    * \param fracs
    * \param cells
    */
    static int CIC_Cells_Fracs (const Particle<NReal, NInt>& p,
                                const Real*                  plo,
                                const Real*                  dx,
                                Vector<Real>&                 fracs,
                                Vector<IntVect>&              cells);

    /**
    * \brief Does CIC computations for arbitrary particle/grid dx's.
    *
    * \param p
    * \param plo
    * \param dx_geom
    * \param dx_part
    * \param fracs
    * \param cells
    */
    static int CIC_Cells_Fracs (const Particle<NReal, NInt>& p,
                                const Real*                  plo,
                                const Real*                  dx_geom,
                                const Real*                  dx_part,
                                Vector<Real>&                 fracs,
                                Vector<IntVect>&              cells);
};
template <int NReal, int NInt> int Particle<NReal, NInt>::the_next_id = 1;

template<int NReal, int NInt>
inline
void
Particle<NReal, NInt>::CIC_Fracs (const Real* frac, Real* fracs)
{
    //
    // "frac"  should be dimensioned: Real frac[AMREX_SPACEDIM]
    //
    // "fracs" should be dimensioned: Real fracs[AMREX_D_TERM(2,+2,+4)]
    //
#if (AMREX_SPACEDIM == 1)

    //! High
    fracs[0] = frac[0];

    // Low
    fracs[1] = (1-frac[0]);

#elif (AMREX_SPACEDIM == 2)

    //! HH
    fracs[0] = frac[0] * frac[1] ;

    // LH
    fracs[1] = (1-frac[0]) * frac[1];

    // LL
    fracs[2] = (1-frac[0]) * (1-frac[1]);

    // HL
    fracs[3] = frac[0] * (1-frac[1]);

#elif (AMREX_SPACEDIM == 3)

    //! HHH
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

template<int NReal, int NInt>
inline
void
Particle<NReal, NInt>::CIC_Cells (const IntVect& hicell, IntVect* cells)
{

    //! "cells" should be dimensioned: IntVect cells[AMREX_D_TERM(2,+2,+4)]
    IntVect cell = hicell;

#if (AMREX_SPACEDIM == 1)

    //! High
    cells[0] = cell;

    // Low
    cell[0]  = cell[0] - 1;
    cells[1] = cell;

#elif (AMREX_SPACEDIM == 2)

    //! HH
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

#elif (AMREX_SPACEDIM == 3)

    //! HHH
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

template<int NReal, int NInt>
inline
int
Particle<NReal, NInt>::CIC_Cells_Fracs (const Particle<NReal, NInt>& p,
                                        const Real*         plo,
                                        const Real*         dx,
                                        Vector<Real>&        fracs,
                                        Vector<IntVect>&     cells)
{
    return Particle<NReal, NInt>::CIC_Cells_Fracs(p,plo,dx,dx,fracs,cells);
}

template <int NReal, int NInt>
void
Particle<NReal, NInt>::CIC_Cells_Fracs_Basic (const Particle<NReal, NInt>& p,
                                              const Real*         plo,
                                              const Real*         dx,
                                              Real*               fracs,
                                              IntVect*            cells)
{
    BL_PROFILE("Particle<NReal, NInt>::CIC_Cells_Fracs_B()");
    //
    // "fracs" should be dimensioned: Real    fracs[AMREX_D_TERM(2,+2,+4)]
    //
    // "cells" should be dimensioned: IntVect cells[AMREX_D_TERM(2,+2,+4)]
    //
    const Real len[AMREX_SPACEDIM] = { AMREX_D_DECL((p.m_rdata.pos[0]-plo[0])/dx[0] + Real(0.5),
                                           (p.m_rdata.pos[1]-plo[1])/dx[1] + Real(0.5),
                                           (p.m_rdata.pos[2]-plo[2])/dx[2] + Real(0.5)) };

    const IntVect cell(AMREX_D_DECL(floor(len[0]), floor(len[1]), floor(len[2])));

    const Real frac[AMREX_SPACEDIM] = { AMREX_D_DECL(len[0]-cell[0], len[1]-cell[1], len[2]-cell[2]) };

    Particle<NReal, NInt>::CIC_Fracs(frac, fracs);
    Particle<NReal, NInt>::CIC_Cells(cell, cells);
}

template <int NReal, int NInt>
int
Particle<NReal, NInt>::CIC_Cells_Fracs (const Particle<NReal, NInt>& p,
                                        const Real*         plo,
                                        const Real*         dx_geom,
                                        const Real*         dx_part,
                                        Vector<Real>&        fracs,
                                        Vector<IntVect>&     cells)
{
    BL_PROFILE("Particle<NReal, NInt>::CIC_Cells_Fracs()");

    if (dx_geom == dx_part)
    {
        const int M = AMREX_D_TERM(2,+2,+4);
        fracs.resize(M);
        cells.resize(M);
        Particle<NReal, NInt>::CIC_Cells_Fracs_Basic(p,plo,dx_geom,fracs.dataPtr(),cells.dataPtr());
        return M;
    }
    //
    // The first element in fracs and cells is the lowest corner, the last is the highest.
    //
    const Real hilen[AMREX_SPACEDIM] = { AMREX_D_DECL((p.m_rdata.pos[0]-plo[0]+dx_part[0]/2)/dx_geom[0],
                                             (p.m_rdata.pos[1]-plo[1]+dx_part[1]/2)/dx_geom[1],
                                             (p.m_rdata.pos[2]-plo[2]+dx_part[2]/2)/dx_geom[2]) };

    const Real lolen[AMREX_SPACEDIM] = { AMREX_D_DECL((p.m_rdata.pos[0]-plo[0]-dx_part[0]/2)/dx_geom[0],
                                             (p.m_rdata.pos[1]-plo[1]-dx_part[1]/2)/dx_geom[1],
                                             (p.m_rdata.pos[2]-plo[2]-dx_part[2]/2)/dx_geom[2]) };

    const IntVect hicell(AMREX_D_DECL(floor(hilen[0]), floor(hilen[1]), floor(hilen[2])));

    const IntVect locell(AMREX_D_DECL(floor(lolen[0]), floor(lolen[1]), floor(lolen[2])));

    const Real cell_density = AMREX_D_TERM(dx_geom[0]/dx_part[0],*dx_geom[1]/dx_part[1],*dx_geom[2]/dx_part[2]);

    const int M = AMREX_D_TERM((hicell[0]-locell[0]+1),*(hicell[1]-locell[1]+1),*(hicell[2]-locell[2]+1));

    fracs.resize(M);
    cells.resize(M);

    //! This portion might be slightly inefficient. Feel free to redo it if need be.
    int i = 0;
#if (AMREX_SPACEDIM == 1)
    for (int xi = locell[0]; xi <= hicell[0]; xi++)
    {
        cells[i][0] = xi;
        fracs[i] = (std::min(hilen[0]-xi,Real(1))-std::max(lolen[0]-xi,Real(0)))*cell_density;
        i++;
    }
#elif (AMREX_SPACEDIM == 2)
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
#elif (AMREX_SPACEDIM == 3)
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

template <int NReal, int NInt>
int
Particle<NReal, NInt>::MaxReaders ()
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

template <int NReal, int NInt>
long
Particle<NReal, NInt>::MaxParticlesPerRead ()
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

template <int NReal, int NInt>
const std::string&
Particle<NReal, NInt>::DataPrefix ()
{
    //
    // The actual particle data is stored in files of the form: DATA_nnnn.
    //
    static const std::string data("DATA_");

    return data;
}

template <int NReal, int NInt>
const std::string&
Particle<NReal, NInt>::Version ()
{
    //
    // If we change the Checkpoint/Restart format we should increment this.
    //
    // Previous version strings:
    //
    //    "Version_One_Dot_Zero"
    //    "Version_One_Dot_One"
    //
    static const std::string version("Version_Two_Dot_Zero");

    return version;
}

template <int NReal, int NInt>
int
Particle<NReal, NInt>::NextID ()
{
    int next;

#ifdef _OPENMP
#if (__GNUC__ < 5)
#pragma omp critical (amrex_particle_nextid)
#else
#pragma omp atomic capture
#endif
#endif
    next = the_next_id++;

    if (next > LastParticleID)
	amrex::Abort("Particle<NReal, NInt>::NextID() -- too many particles");

    return next;
}

template <int NReal, int NInt>
int
Particle<NReal, NInt>::UnprotectedNextID ()
{
    int next = the_next_id++;
    if (next > LastParticleID)
	amrex::Abort("Particle<NReal, NInt>::NextID() -- too many particles");
    return next;
}

template <int NReal, int NInt>
void
Particle<NReal, NInt>::NextID (int nextid)
{
    the_next_id = nextid;
}

template <int NReal, int NInt>
Real
Particle<NReal, NInt>::InterpDoit (const FArrayBox& fab,
			      const IntVect&   cell,
			      const Real*      frac,
			      int              comp)
{
    BL_PROFILE("Particle<NReal, NInt>::InterpDoit(fcfc)");
    const int M = AMREX_D_TERM(2,+2,+4);

    Real    fracs[M];
    IntVect cells[M];

    Particle<NReal, NInt>::CIC_Fracs(frac, fracs);
    Particle<NReal, NInt>::CIC_Cells(cell, cells);

    Real val = Particle<NReal, NInt>::InterpDoit(fab,fracs,cells,comp);

    return val;
}

template <int NReal, int NInt>
Real
Particle<NReal, NInt>::InterpDoit (const FArrayBox& fab,
			      const Real*      fracs,
			      const IntVect*   cells,
			      int              comp)
{
    BL_PROFILE("Particle<NReal, NInt>::InterpDoit(ffcc)");
    const int M = AMREX_D_TERM(2,+2,+4);

    Real val = 0;

    for (int i = 0; i < M; i++)
    {
        val += fab(cells[i],comp) * fracs[i];
    }

    return val;
}

template <int NReal, int NInt>
void
Particle<NReal, NInt>::Interp (const Particle<NReal, NInt>& prt,
			  const Geometry&     geom,
			  const FArrayBox&    fab,
			  const int*          idx,
			  Real*               val,
			  int                 cnt)
{
    BL_PROFILE("Particle<NReal, NInt>::Interp()");
    BL_ASSERT(idx != 0);
    BL_ASSERT(val != 0);

    const int       M   = AMREX_D_TERM(2,+2,+4);
    const Real*     plo = geom.ProbLo();
    const Real*     dx  = geom.CellSize();

    Real    fracs[M];
    IntVect cells[M];
    //
    // Get "fracs" and "cells".
    //
    Particle<NReal, NInt>::CIC_Cells_Fracs_Basic(prt, plo, dx, fracs, cells);

    for (int i = 0; i < cnt; i++)
    {
        BL_ASSERT(idx[i] >= 0 && idx[i] < fab.nComp());

        val[i] = Particle<NReal, NInt>::InterpDoit(fab,fracs,cells,idx[i]);
    }
}

template <int NReal, int NInt>
void
Particle<NReal, NInt>::GetGravity (const FArrayBox&    gfab,
			      const Geometry&     geom,
			      const Particle<NReal, NInt>& p,
			      Real*               grav)
{
    BL_PROFILE("Particle<NReal, NInt>::GetGravity()");
    BL_ASSERT(grav != 0);

    int idx[AMREX_SPACEDIM] = { AMREX_D_DECL(0,1,2) };

    Particle<NReal, NInt>::Interp(p,geom,gfab,idx,grav,AMREX_SPACEDIM);
}

template <int NReal, int NInt>
std::ostream&
operator<< (std::ostream& os, const Particle<NReal, NInt>& p)
{
    os << p.m_idata.id   << ' '
       << p.m_idata.cpu  << ' ';

    for (int i = 0; i < AMREX_SPACEDIM + NReal; i++)
        os << p.m_rdata.arr[i] << ' ';

    for (int i = 2; i < 2 + NInt; i++)
        os << p.m_idata.arr[i] << ' ';

    if (!os.good())
        amrex::Error("operator<<(ostream&,Particle<NReal, NInt>&) failed");

    return os;
}

} // namespace amrex

#endif // AMREX_PARTICLE_H_
