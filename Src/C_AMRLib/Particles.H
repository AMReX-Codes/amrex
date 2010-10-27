#ifndef _PARTICLES_H_
#define _PARTICLES_H_ 

#include <deque>
#include <vector>
#include <fstream>
#include <iostream>

#include "REAL.H"
#include "IntVect.H"
#include "Array.H"
#include "Amr.H"
#include "AmrLevel.H"
#include "Utility.H"
#include "Geometry.H"
#include "VisMF.H"

struct ParticleBase
{
    int     m_id;
    int     m_cpu;
    int     m_lev;
    int     m_grid;
    IntVect m_cell;
    Real    m_pos[BL_SPACEDIM];

    ParticleBase ()
        :
        m_id(-1),
        m_cpu(-1),
        m_lev(-1),
        m_grid(-1)
        {}

    static IntVect Index (const ParticleBase& prt, int lev, const Amr* amr);

    static bool Where (ParticleBase& prt, const Amr* amr, bool update = false);

    static void PeriodicShift (ParticleBase& prt, const Amr* amr);

    static const std::string& Version ();

    static const std::string& DataPrefix ();

    static void GetGravity (const FArrayBox& gfab, const Amr* amr, int lev, const ParticleBase& p, Real* grav);
};

std::ostream& operator<< (std::ostream& os, const ParticleBase& p);

template <int N>
struct Particle
    :
    public ParticleBase
{
    //
    // The amount of Real data we hold.
    //
    // In some cases this is:
    // 0 - particle mass
    // 1 - x-velocity
    // 2 - y-velocity
    // 3 - z-velocity
    Real m_data[N];
};

template <int N>
class ParticleContainer
{
public:
    //
    // The type of Particles we hold.
    //
    typedef Particle<N> ParticleType;
    //
    // Types for iterating over the particles.
    //
    typedef typename std::deque<ParticleType>::iterator iterator;

    typedef typename std::deque<ParticleType>::const_iterator const_iterator;

    ParticleContainer (Amr* amr)
        :
        m_verbose(1), m_amr(amr) { BL_ASSERT(amr != 0); }

    void Init (const std::string& file);

    void InitRandom (long icount, unsigned long iseed, Real particleMass);

    Real sumParticleMass (const MultiFab& mf, int level) const;
    void sumParticleMomentum (const MultiFab& mf, int lev, Real* mom) const;

    void Increment (MultiFab& mf, int level) const;

    void AssignDensity (MultiFab& mf, int level) const;
 
    Real estTimestep (const MultiFab& grav_vector, int level) const;

    void Redistribute (bool where_already_called = false);

    bool OK () const;

    void ByteSpread () const;

    void MoveRandom ();

    // **************************************************************************************************************** 
    //
    // If the particles move only with self-gravity from themselves and the gas, then 
    // we can move them according to the schemes below.
    // The gravitational force must be computed between the calls of the parts of the integration scheme.
    //
    // The following two functions form a PREDICTOR CORRECTOR scheme for integrating the motion of the particles
    //
    void movePredict (const MultiFab& grav_vector, int level, Real timestep);
    void moveCorrect (const MultiFab& grav_vector_old, const MultiFab& grav_vector, int level, Real timestep);
    //
    // TODO: the methods should return a constraint on the timestep...
    //
    // The following two functions form a KICK DRIFT KICK scheme for integrating the motion of the particles
    //
    void moveKickDrift (const MultiFab& grav_vector, int level, Real timestep);
    void moveKick      (const MultiFab& grav_vector, int level, Real timestep);
    //
    // after the moveKickDrift step the positions of the particles are advanced for a full timestep,
    // so this scheme should work in the overall algorithm...
    //
    // **************************************************************************************************************** 

    void Checkpoint (const std::string& dir, const std::string& name) const;

    void Restart (const std::string& dir, const std::string& file);

    int Verbose () { return m_verbose; }

    void Verbose (int verbose) { m_verbose = verbose; }

protected:
    //
    // Helper function for Checkpoint().
    //
    void WriteParticles (int            level,
                         std::ofstream& ofs,
                         int            fnum,
                         Array<int>&    which,
                         Array<int>&    count,
                         Array<long>&   where) const;
    //
    // Helper functions for Restart().
    //
    void Restart_OneDotZero (const std::string& fullname,
                             std::ifstream&     HdrFile);

    void ReadParticles_OneDotZero (int            cnt,
                                   int            grd,
                                   int            lev,
                                   std::ifstream& ifs);
    //
    // The data.
    //
    int                      m_verbose;
    Amr*                     m_amr;
    std::deque<ParticleType> m_particles;
};

template <int N>
void
ParticleContainer<N>::ByteSpread () const
{
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    long mn = sizeof(ParticleType) * m_particles.size(), mx = mn;

    ParallelDescriptor::ReduceLongMin(mn,IOProc);
    ParallelDescriptor::ReduceLongMax(mx,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "\nParticleContainer<N> byte spread across MPI nodes: ["
                  << mn
                  << " ... "
                  << mx
                  << "]\n\n";
    }
}

template <int N>
void
ParticleContainer<N>::Init (const std::string& file)
{
    BL_ASSERT(!file.empty());

    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();

    std::ifstream ifs;

    ifs.open(file.c_str(), std::ios::in);

    if (!ifs.good())
        BoxLib::FileOpenFailed(file);

    int id  = 1;
    int cnt = 0;

    ifs >> cnt;

    ParticleType p;

    for (int i = 0; i < cnt; i++)
    {
        //
        // We don't read in m_id or m_cpu.  We'll set those later
        // in a manner to guarantee the global uniqueness of the pair.
        //
        for (int i = 0; i < BL_SPACEDIM; i++)
            ifs >> p.m_pos[i];

        for (int i = 0; i < N; i++)
            ifs >> p.m_data[i];

        if (!ifs.good())
        {
            std::string msg("ParticleContainer::Init(");
            msg += file;
            msg += ") failed";
            BoxLib::Error(msg.c_str());
        }

        if (!ParticleBase::Where(p,m_amr))
        {
            if (ParallelDescriptor::IOProcessor())
            {
                std::cout << "Ignoring particle @ pos:";
                for (int i = 0; i < BL_SPACEDIM; i++)
                    std::cout << p.m_pos[i] << ' ';
                std::cout << "-- not in domain !!!\n";
            }
        }
        else
        {
            //
            // Only the MPI process that owns the needed velocity data gets "p".
            //
            const DistributionMapping& dm = m_amr->getLevel(p.m_lev).get_new_data(0).DistributionMap();

            if (dm[p.m_grid] == MyProc)
            {
                //
                // Note that m_id always starts at 1.
                //
                p.m_id  = id++;
                p.m_cpu = MyProc;

                m_particles.push_back(p);
            }
        }
    }

    for (int i = 0; i < NProcs; i++)
    {
        if (MyProc == i)
        {
            std::cout << "Processor "
                      << i
                      << " has "
                      << m_particles.size()
                      << " particles in Init()" << std::endl;
        }
        ParallelDescriptor::Barrier();
    }
    if (ParallelDescriptor::IOProcessor())
        std::cout << '\n';

    BL_ASSERT(OK());

    ByteSpread();
}

template <int N>
void
ParticleContainer<N>::InitRandom (long          icount,
                                  unsigned long iseed,
                                  Real          mass)
{
    BL_ASSERT(iseed  > 0);
    BL_ASSERT(icount > 0);

    BL_ASSERT(m_amr != 0);

    const int  MyProc = ParallelDescriptor::MyProc();
    const int  NProcs = ParallelDescriptor::NProcs();
    //
    // We'll generate the particles in parallel.
    //
    // Each CPU will key off the given seed to get independent streams of random numbers.
    //
    BoxLib::InitRandom(iseed+MyProc);

    long M = icount / NProcs;
    //
    // Processor 0 will get the slop.
    //
    if (MyProc == 0)
        M += (icount % NProcs);

    Real len[BL_SPACEDIM], r;

    const Geometry& geom = m_amr->Geom(0);

    D_TERM(len[0]=geom.ProbLength(0);,
           len[1]=geom.ProbLength(1);,
           len[2]=geom.ProbLength(2););

    for (long icnt = 0; icnt < M; icnt++)
    {
        ParticleType p;

        for (int i = 0; i < BL_SPACEDIM; i++)
        {
            do
            {
                r = BoxLib::Random();
            }
            while (r == 0 || r == 1);

            p.m_pos[i] = geom.ProbLo(i) + (r * len[i]);

            BL_ASSERT(p.m_pos[i] < geom.ProbHi(i));
        }

        p.m_data[0] = mass;
        for (int i = 1; i < N; i++)
            //
            // Just zero out the rest of the data for lack of a better value.
            //
            p.m_data[i] = 0;

        p.m_id  = icnt+1;
        p.m_cpu = MyProc;

        m_particles.push_back(p);

        BL_ASSERT(ParticleBase::Where(p,m_amr));
    }
    //
    // Let Redistribute() sort out where the particles belong.
    //
    Redistribute();
}

template <int N>
void
ParticleContainer<N>::MoveRandom ()
{
    //
    // Move particles up to 0.25*CellSize distance in each coordinate direction.
    //
    BL_ASSERT(OK());

    BL_ASSERT(m_amr != 0);

    for (iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        ParticleType& p = *it;

        const Geometry& geom = m_amr->Geom(p.m_lev);

        for (int i = 0; i < BL_SPACEDIM; i++)
        {
            int sgn = 1;

            if (BoxLib::Random() >= 0.5)
                sgn = -1;

            p.m_pos[i] +=  sgn * 0.25 * geom.CellSize(i) * BoxLib::Random();
        }

        if (!ParticleBase::Where(p,m_amr,true))
        {
            //
            // Here's where we need to deal with boundary conditions.
            //
            // For the moment we assume all boundaries are periodic.
            //
            // Shift the particle back into the domain.
            //
            ParticleBase::PeriodicShift(p,m_amr);
            //
            // The particle should now be back in the problem domain.
            //
        }
    }

    if (m_verbose > 1 && ParallelDescriptor::IOProcessor())
        std::cout << "Calling Redistribute() from ParticleContainer<N>::MoveRandom () ...\n";

    Redistribute(true);
}

template <int N>
void
ParticleContainer<N>::Increment (MultiFab& mf,
                                 int       lev) const
{
    BL_ASSERT(OK());

    for (const_iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        const ParticleType& p = *it;

        if (p.m_lev == lev)
        {
            mf[p.m_grid](p.m_cell) += 1;
        }
    }
}

template <int N>
Real
ParticleContainer<N>::estTimestep (const MultiFab& gv,
                                   int             lev) const
{
    BL_ASSERT(N >= BL_SPACEDIM+1);

    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_amr->Geom(0);
    const Real*     dx       = geom.CellSize();
    const int       n        = m_particles.size();
    Real            dt       = 1e50;
    Real            dt_local;

#ifdef BL_USE_OMP
#pragma omp parallel private(dt_local)
#endif
    {
        dt_local = 1e50;

#ifdef BL_USE_OMP
#pragma omp for
#endif
        for (int i = 0; i < n; i++)
        {
            const ParticleType& p = m_particles[i];

            if (p.m_lev == lev)
            {
                const Real mag_vel = sqrt(D_TERM(p.m_data[1]*p.m_data[1], + p.m_data[2]*p.m_data[2], + p.m_data[3]*p.m_data[3]));

                Real dt_part = (mag_vel > 0) ? 0.1 * dx[0] / mag_vel : 1e50;

                const FArrayBox& gfab = gv[p.m_grid];

                const Real gval[BL_SPACEDIM] = { D_DECL(gfab(p.m_cell,0), gfab(p.m_cell,1), gfab(p.m_cell,2)) };

                const Real mag_grav = sqrt(D_TERM(gval[0]*gval[0], + gval[1]*gval[1], + gval[2]*gval[2]));

                if (mag_grav > 0)
                    dt_part = std::min( dt_part, 1/sqrt(mag_grav/dx[0]) );

                dt_local = std::min(dt_part, dt_local);
            }
        }

#ifdef BL_USE_OMP
#pragma omp critical(estTimestep_lock)
#endif
        {
            dt = std::min(dt, dt_local);
        }
    }

    ParallelDescriptor::ReduceRealMin(dt);

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      stoptime = ParallelDescriptor::second() - strttime;

    ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "ParticleContainer<N>::estTimestep() Time: " << stoptime << "\n\n";
    }

    return dt;
}

//
// Assumes mass is in m_data[0]!
//

template <int N>
Real
ParticleContainer<N>::sumParticleMass (const MultiFab& mf,
                                       int             lev) const
{
    BL_ASSERT(N >= 1);

    Real massSum = 0.0;

    for (const_iterator it = m_particles.begin(); it != m_particles.end(); ++it)
        if (it->m_lev == lev)
            massSum += it->m_data[0];

    ParallelDescriptor::ReduceRealSum(massSum);

    return massSum;
}

//
// Assumes mass is in m_data[0], vx in m_dat[1], ...!
// dim defines the cartesian direction in which the momentum is summed, x is 0, y is 1, ...
//

template <int N>
void
ParticleContainer<N>::sumParticleMomentum (const MultiFab& mf,
                                           int             lev,
                                           Real*           mom) const
{
    BL_ASSERT(N >= BL_SPACEDIM+1);

    Real D_DECL(mom_0, mom_1, mom_2);

    const int n = m_particles.size();

#ifdef BL_USE_OMP
#pragma omp parallel for reduction(+:D_DECL(mom_0, mom_1, mom_2))
#endif
    for (int i = 0; i < n; i++)
    {
        const ParticleType& p = m_particles[i];

        if (p.m_lev == lev)
        {
            D_TERM(mom_0 += p.m_data[0] * p.m_data[1];,
                   mom_1 += p.m_data[0] * p.m_data[2];,
                   mom_2 += p.m_data[0] * p.m_data[3];);
        }
    }

    D_TERM(mom[0] = mom_0;, mom[1] = mom_1;, mom[2] = mom_2;);

    ParallelDescriptor::ReduceRealSum(mom,BL_SPACEDIM);
}

//
// Some cic-scheme.
// Will get ugly as soon as there is more than just the rootgrid
//

template <int N>
void
ParticleContainer<N>::AssignDensity (MultiFab& mf,
                                     int       lev) const
{
    BL_ASSERT(N >= 1);

    const Real strttime = ParallelDescriptor::second();

    const Geometry& geom = m_amr->Geom(lev);

    if (!geom.isAllPeriodic())
        BoxLib::Error("AssignDensity: Problem has to be periodic...");

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        mf[mfi].setVal(0);

    const Real* dx = geom.CellSize();
    //
    // This is the # of grids we own in the MultiFab.
    //
    const int n = mf.IndexMap().size();
    //
    // This is a little funky.  What in effect this'll do is force
    // each thread to work on a single (separate) grid at a time.  That
    // way no thread will step on any other.  If there's only one grid,
    // then oh well ....
    //
#ifdef BL_USE_OMP
#pragma omp parallel for schedule(dynamic,1) if (n>1)
#endif
    for (int i = 0; i < n; i++)
    {
        const int grid = mf.IndexMap()[i];

        for (const_iterator it = m_particles.begin(); it != m_particles.end(); ++it)
        {
            const ParticleType& p = *it;

            if (p.m_grid == grid && p.m_lev == lev)
            {
                const IntVect csect(D_DECL(floor((p.m_pos[0]-geom.ProbLo(0))/dx[0] + 0.5),
                                           floor((p.m_pos[1]-geom.ProbLo(1))/dx[1] + 0.5),
                                           floor((p.m_pos[2]-geom.ProbLo(2))/dx[2] + 0.5)));

                const Real frac[BL_SPACEDIM] = { D_DECL(-csect[0] + p.m_pos[0]/dx[0] + 0.5,
                                                        -csect[1] + p.m_pos[1]/dx[1] + 0.5,
                                                        -csect[2] + p.m_pos[2]/dx[2] + 0.5) };
                FArrayBox& fab = mf[p.m_grid];

                IntVect cell = csect;

#if (BL_SPACEDIM == 1)
                // High
                fab(cell) += p.m_data[0] * frac[0];

                // Low
                cell[0]     = cell[0] - 1;
                fab(cell)  += p.m_data[0] * (1-frac[0]);
#endif
            
#if (BL_SPACEDIM == 2)
                // HH
                fab(cell) += p.m_data[0] * frac[0] * frac[1] ;
    
                // LH
                cell[0]    = cell[0] - 1;
                fab(cell) += p.m_data[0] * (1-frac[0]) * frac[1] ;
    
                // LL
                cell[1]    = cell[1] - 1;
                fab(cell) += p.m_data[0] * (1-frac[0]) * (1-frac[1]);
    
                // HL
                cell[0]    = cell[0] + 1;
                fab(cell) += p.m_data[0] * frac[0] * (1-frac[1]);
#endif

#if (BL_SPACEDIM == 3)

                // HHH
                fab(cell) += p.m_data[0] * frac[0] * frac[1] * frac[2] ;

                // LHH
                cell[0]    = cell[0] - 1;
                fab(cell) += p.m_data[0] * (1-frac[0]) * frac[1] * frac[2] ;

                // LLH
                cell[1]    = cell[1] - 1;
                fab(cell) += p.m_data[0] * frac[0] * (1-frac[1]) * frac[2] ;
    
                // HLH
                cell[0]    = cell[0] + 1;
                fab(cell) += p.m_data[0] * (1-frac[0]) * (1-frac[1]) *    frac[2] ;

                cell     = csect;
                cell[2]  = cell[2] - 1;

                // HHL
                fab(cell) += p.m_data[0] * frac[0] * frac[1] * (1-frac[2]);
    
                // LHL
                cell[0]    = cell[0] - 1;
                fab(cell) += p.m_data[0] * (1-frac[0]) * frac[1] * (1-frac[2]);

                // LLL
                cell[1]    = cell[1] - 1;
                fab(cell) += p.m_data[0] * frac[0] * (1-frac[1]) * (1-frac[2]);
    
                // HLL
                cell[0]    = cell[0] + 1;
                fab(cell) += p.m_data[0] * (1-frac[0]) * (1-frac[1]) * (1-frac[2]);
#endif
            }
        }
    }

    mf.SumBoundary();
    geom.SumPeriodicBoundary(mf);

    const Real cell_vol = D_TERM(dx[0], * dx[1], * dx[2]);

    mf.mult(1/cell_vol);

    if (m_verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::AssignDensity() time: " << stoptime << "\n\n";
        }
    }
}

template <int N>
void
ParticleContainer<N>::movePredict (const MultiFab& gv,
                                   int             lev,
                                   Real            dt)
{
    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);

    const Real strttime = ParallelDescriptor::second();

    const int n = m_particles.size();

#ifdef BL_USE_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < n; i++)
    {
        ParticleType& p = m_particles[i];

        if (p.m_lev == lev)
        {
            //
            // Note: m_data[0] is mass, 1 is v_x, ...
            //
            const FArrayBox& gfab = gv[p.m_grid];

            D_TERM(p.m_data[1] += dt * gfab(p.m_cell,0);,
                   p.m_data[2] += dt * gfab(p.m_cell,1);,
                   p.m_data[3] += dt * gfab(p.m_cell,2););

            D_TERM(p.m_pos[0]  += dt * p.m_data[1];,
                   p.m_pos[1]  += dt * p.m_data[2];,
                   p.m_pos[2]  += dt * p.m_data[3];);

            if (!ParticleBase::Where(p,m_amr,true))
            {
                //
                // Here's where we need to deal with boundary conditions.
                //
                // For the moment we assume all boundaries are periodic.
                //
                // Shift the particle back into the domain.
                //
                ParticleBase::PeriodicShift(p,m_amr);
                //
                // The particle should now be back in the problem domain.
                //
            }
        }
    }

    if (m_verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::movePredict() time: " << stoptime << "\n\n";
        }
    }

    Redistribute(true);
}

template <int N>
void
ParticleContainer<N>::moveCorrect (const MultiFab& gv_old,
                                   const MultiFab& gv,
                                   int             lev,
                                   Real            dt)
{
    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);

    const Real strttime = ParallelDescriptor::second();

    const int n = m_particles.size();

#ifdef BL_USE_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < n; i++)
    {
        ParticleType& p = m_particles[i];

        if (p.m_lev == lev)
        {
            //
            // Note: m_data[0] is mass, 1 is v_x, ...
            //
            const FArrayBox& gfab     = gv[p.m_grid];
            const FArrayBox& gfab_old = gv_old[p.m_grid];

            D_TERM(p.m_pos[0]  -= 0.5 * dt * p.m_data[1];,
                   p.m_pos[1]  -= 0.5 * dt * p.m_data[2];,
                   p.m_pos[2]  -= 0.5 * dt * p.m_data[3];);
            
            D_TERM(p.m_data[1] += 0.5 * dt * ( gfab(p.m_cell,0) - gfab_old(p.m_cell,0) );,
                   p.m_data[2] += 0.5 * dt * ( gfab(p.m_cell,1) - gfab_old(p.m_cell,1) );,
                   p.m_data[3] += 0.5 * dt * ( gfab(p.m_cell,2) - gfab_old(p.m_cell,2) ););

            D_TERM(p.m_pos[0]  += 0.5 * dt * p.m_data[1];,
                   p.m_pos[1]  += 0.5 * dt * p.m_data[2];,
                   p.m_pos[2]  += 0.5 * dt * p.m_data[3];);

            if (!ParticleBase::Where(p,m_amr,true))
            {
                //
                // Here's where we need to deal with boundary conditions.
                //
                // For the moment we assume all boundaries are periodic.
                //
                // Shift the particle back into the domain.
                //
                ParticleBase::PeriodicShift(p,m_amr);
                //
                // The particle should now be back in the problem domain.
                //
            }
        }
    }

    if (m_verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::moveCorrect() time: " << stoptime << "\n\n";
        }
    }

    Redistribute(true);
}

template <int N>
void
ParticleContainer<N>::moveKickDrift (const MultiFab& gv,
                                     int             lev,
                                     Real            dt)
{
    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);

    const Real strttime = ParallelDescriptor::second();

    const int n = m_particles.size();

#ifdef BL_USE_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < n; i++)
    {
        ParticleType& p = m_particles[i];

        if (p.m_lev == lev)
        {
            //
            // note: m_data[0] is mass, 1 is v_x, ...
            //
            const FArrayBox& gfab = gv[p.m_grid];

            Real grav[BL_SPACEDIM];

            ParticleBase::GetGravity(gfab, m_amr, lev, p, grav);

            D_TERM(p.m_data[1] += 0.5 * dt * grav[0];,
                   p.m_data[2] += 0.5 * dt * grav[1];,
                   p.m_data[3] += 0.5 * dt * grav[2];);

            D_TERM(p.m_pos[0] += dt * p.m_data[1];,
                   p.m_pos[1] += dt * p.m_data[2];,
                   p.m_pos[2] += dt * p.m_data[3];);

            if (!ParticleBase::Where(p,m_amr,true))
            {
                //
                // Here's where we need to deal with boundary conditions.
                //
                // For the moment we assume all boundaries are periodic.
                //
                // Shift the particle back into the domain.
                //
                ParticleBase::PeriodicShift(p,m_amr);
                //
                // Modulo any arithmetic funny business the particle should in the problem domain.
                //
            }
        }
    }

    if (m_verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::moveKickDrift() time: " << stoptime << "\n\n";
        }
    }

    Redistribute(true);
}

template <int N>
void
ParticleContainer<N>::moveKick (const MultiFab& gv,
                                int             lev,
                                Real            dt)
{
    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);

    const Real strttime = ParallelDescriptor::second();

    const int n = m_particles.size();

#ifdef BL_USE_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < n; i++)
    {
        ParticleType& p = m_particles[i];

        if (p.m_lev == lev)
        {
            //
            // Note: m_data[0] is mass, 1 is v_x, ...
            //
            const FArrayBox& gfab = gv[p.m_grid];

            Real grav[BL_SPACEDIM];

            ParticleBase::GetGravity(gfab, m_amr, lev, p, grav);

            D_TERM(p.m_data[1] += 0.5 * dt * grav[0];,
                   p.m_data[2] += 0.5 * dt * grav[1];,
                   p.m_data[3] += 0.5 * dt * grav[2];);
        }
    }

    if (m_verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::moveKick() time: " << stoptime << "\n\n";
        }
    }
    //
    // No need for Redistribution(), we only change the velocity.
    //
}

template <int N>
void
ParticleContainer<N>::Redistribute (bool where_already_called)
{
    const int  MyProc   = ParallelDescriptor::MyProc();
    const Real strttime = ParallelDescriptor::second();

    std::vector<int>      owner;
    std::vector<iterator> index;

    for (iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        ParticleType& p = *it;

        BL_ASSERT(p.m_id > 0);

        if (!where_already_called)
        {
            bool valid = ParticleBase::Where(p,m_amr);

            BL_ASSERT(valid);
        }
        //
        // The owner of the particle is the CPU owning the finest grid
        // in state data that contains the particle.
        //
        const int who = m_amr->getLevel(p.m_lev).get_new_data(0).DistributionMap()[p.m_grid];

        if (who != MyProc)
        {
            //
            // We're not the rightful owner of this particle.
            //
            index.push_back(it);
            owner.push_back(who);
        }
    }

    const int NProcs = ParallelDescriptor::NProcs();

    if (NProcs == 1)
    {
        BL_ASSERT(index.empty());
        BL_ASSERT(owner.empty());

        return;
    }

#if BL_USE_MPI
    //
    // We may now have particles that are rightfully owned by another CPU.
    //
    Array<int> Snds(NProcs,0);
    Array<int> Rcvs(NProcs,0);

    for (int i = 0; i < owner.size(); i++)
        Snds[owner[i]]++;

    BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 Rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 ParallelDescriptor::Communicator()) );
    BL_ASSERT(Rcvs[MyProc] == 0);

    int NumRcvs = 0;
    for (int i = 0; i < NProcs; i++)
        NumRcvs += Rcvs[i];

    int NumSnds = 0;
    for (int i = 0; i < NProcs; i++)
        NumSnds += Snds[i];

    int maxsendrecv = 0;

    Array<ParticleType> nparticles(NumRcvs);

    if (m_verbose > 1)
    {
        maxsendrecv = std::max(NumRcvs,NumSnds);

        ParallelDescriptor::ReduceIntMax(maxsendrecv);

        if (maxsendrecv > 0)
        {
            for (int i = 0; i < NProcs; i++)
            {
                if (MyProc == i)
                {
                    std::cout << "Processor "
                              << i
                              << " : NumRcvs: "
                              << NumRcvs
                              << " NumSnds: "
                              << NumSnds
                              << " in Redistribute()" << std::endl;
                }
                ParallelDescriptor::Barrier();
            }

            if (ParallelDescriptor::IOProcessor())
                std::cout << '\n';
        }
    }

    {
        //
        // The int data.
        //
        const int iChunkSize = BL_SPACEDIM+4;

        Array<int> recvdata (NumRcvs*iChunkSize);
        Array<int> senddata (NumSnds*iChunkSize);

        Array<int> sendcnts(NProcs,0), sdispls(NProcs,0);
        Array<int> recvcnts(NProcs,0), rdispls(NProcs,0), offset(NProcs,0);

        for (int i = 0; i < NProcs; i++)
        {
            recvcnts[i] = Rcvs[i] * iChunkSize;
            sendcnts[i] = Snds[i] * iChunkSize;
        }

        for (int i = 1; i < NProcs; i++)
        {
            offset [i] = offset [i-1] + sendcnts[i-1];
            rdispls[i] = rdispls[i-1] + recvcnts[i-1];
            sdispls[i] = sdispls[i-1] + sendcnts[i-1];
        }

        for (int i = 0; i < index.size(); i++)
        {
            int idx = offset[owner[i]];

            ParticleType& p = *index[i];

            senddata[idx+0] = p.m_id;
            senddata[idx+1] = p.m_cpu;
            senddata[idx+2] = p.m_lev;
            senddata[idx+3] = p.m_grid;

            idx += 4;

            D_TERM(senddata[idx+0]=p.m_cell[0];,
                   senddata[idx+1]=p.m_cell[1];,
                   senddata[idx+2]=p.m_cell[2];);

            offset[owner[i]] += iChunkSize;
        }

        BL_MPI_REQUIRE( MPI_Alltoallv(NumSnds == 0 ? 0 : senddata.dataPtr(),
                                      sendcnts.dataPtr(),
                                      sdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<int>::type(),
                                      NumRcvs == 0 ? 0 : recvdata.dataPtr(),
                                      recvcnts.dataPtr(),
                                      rdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<int>::type(),
                                      ParallelDescriptor::Communicator()) );
        //
        // Unpack data into the new particles.
        //
        for (int i = 0, idx = 0; i < NumRcvs; i++)
        {
            ParticleType& p = nparticles[i];

            p.m_id   = recvdata[idx+0];
            p.m_cpu  = recvdata[idx+1];
            p.m_lev  = recvdata[idx+2];
            p.m_grid = recvdata[idx+3];

            idx += 4;

            D_TERM(p.m_cell[0]=recvdata[idx+0];,
                   p.m_cell[1]=recvdata[idx+1];,
                   p.m_cell[2]=recvdata[idx+2];);

            idx += BL_SPACEDIM;

            BL_ASSERT(m_amr->getLevel(p.m_lev).get_new_data(0).DistributionMap()[p.m_grid] == MyProc);
        }
    }

    {
        //
        // The Real data (m_pos & m_data).
        //
        const int rChunkSize = BL_SPACEDIM+N;

        Array<Real> recvdata (NumRcvs*rChunkSize);
        Array<Real> senddata (NumSnds*rChunkSize);

        Array<int> sendcnts(NProcs,0), sdispls(NProcs,0);
        Array<int> recvcnts(NProcs,0), rdispls(NProcs,0), offset(NProcs,0);

        for (int i = 0; i < NProcs; i++)
        {
            recvcnts[i] = Rcvs[i] * rChunkSize;
            sendcnts[i] = Snds[i] * rChunkSize;
        }

        for (int i = 1; i < NProcs; i++)
        {
            offset [i] = offset [i-1] + sendcnts[i-1];
            rdispls[i] = rdispls[i-1] + recvcnts[i-1];
            sdispls[i] = sdispls[i-1] + sendcnts[i-1];
        }

        for (int i = 0; i < index.size(); i++)
        {
            int idx = offset[owner[i]];

            ParticleType& p = *index[i];

            D_TERM(senddata[idx+0]=p.m_pos[0];,
                   senddata[idx+1]=p.m_pos[1];,
                   senddata[idx+2]=p.m_pos[2];);

            idx += BL_SPACEDIM;

            for (int j = 0; j < N; j++)
                senddata[idx+j] = p.m_data[j];

            offset[owner[i]] += rChunkSize;
        }

        BL_MPI_REQUIRE( MPI_Alltoallv(NumSnds == 0 ? 0 : senddata.dataPtr(),
                                      sendcnts.dataPtr(),
                                      sdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<Real>::type(),
                                      NumRcvs == 0 ? 0 : recvdata.dataPtr(),
                                      recvcnts.dataPtr(),
                                      rdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<Real>::type(),
                                      ParallelDescriptor::Communicator()) );
        //
        // Remove particles we no longer own.
        //
        for (int i = 0; i < index.size(); i++)
        {
            ParticleType& p = *index[i];

            BL_ASSERT(p.m_id > 0);

            p.m_id = -p.m_id;
        }

        std::deque<ParticleType> tmp;

        while (!m_particles.empty())
        {
            if (m_particles.front().m_id > 0)
                tmp.push_back(m_particles.front());
            m_particles.pop_front();
        }

        m_particles.swap(tmp);
        //
        // Unpack data into the new particles.
        //
        for (int i = 0, idx = 0; i < NumRcvs; i++)
        {
            ParticleType& p = nparticles[i];

            D_TERM(p.m_pos[0]=recvdata[idx+0];,
                   p.m_pos[1]=recvdata[idx+1];,
                   p.m_pos[2]=recvdata[idx+2];);

            idx += BL_SPACEDIM;

            for (int j = 0; j < N; j++)
                p.m_data[j] = recvdata[idx+j];

            idx += N;

            m_particles.push_back(p);
        }
    }

    BL_ASSERT(OK());

    if (m_verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        if (m_verbose > 1 && maxsendrecv > 0)
        {
            for (int i = 0; i < NProcs; i++)
            {
                if (MyProc == i)
                {
                    std::cout << "Processor "
                              << i
                              << " has "
                              << m_particles.size()
                              << " particles after Redistribute()" << std::endl;
                }
                ParallelDescriptor::Barrier();
            }

            if (ParallelDescriptor::IOProcessor())
                std::cout << '\n';
        }

        ByteSpread();

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::Redistribute() time: " << stoptime << "\n\n";
        }
    }
#endif
}

template <int N>
bool
ParticleContainer<N>::OK () const
{
    //
    // Check that the integer data in each valid particle is what it should be.
    //
    for (const_iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        ParticleType p = *it;

        BL_ASSERT(p.m_id > 0);

        const int     lev  = p.m_lev;
        const int     grid = p.m_grid;
        const IntVect cell = p.m_cell;

        if (!ParticleBase::Where(p,m_amr))
            return false;

        if (lev != p.m_lev || grid != p.m_grid || cell != p.m_cell)
            return false;
    }
    return true;
}

template <int N>
void
ParticleContainer<N>::Checkpoint (const std::string& dir,
                                  const std::string& name) const
{
    BL_ASSERT(OK());

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
        if (!BoxLib::UtilCreateDirectory(pdir, 0755))
            BoxLib::CreateDirectoryFailed(pdir);
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

    long nparticles = m_particles.size();

    ParallelDescriptor::ReduceLongSum(nparticles,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        std::string HdrFileName = pdir;

        if (!HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/')
            HdrFileName += '/';

        HdrFileName += "Header";

        HdrFile.open(HdrFileName.c_str(), std::ios::out|std::ios::trunc);

        if (!HdrFile.good())
            BoxLib::FileOpenFailed(HdrFileName);
        //
        // First thing written is our Checkpoint/Restart version string.
        //
        HdrFile << ParticleBase::Version() << '\n';
        //
        // BL_SPACEDIM and N for sanity checking.
        //
        HdrFile << BL_SPACEDIM << '\n';

        HdrFile << N << '\n';
        //
        // The total number of particles.
        //
        HdrFile << nparticles << '\n';
        //
        // Then the finest level of the AMR hierarchy.
        //
        HdrFile << m_amr->finestLevel() << '\n';
        //
        // Then the number of grids at each level.
        //
        for (int lev = 0; lev <= m_amr->finestLevel(); lev++)
        {
            HdrFile << m_amr->boxArray(lev).size() << '\n';
        }
    }
    //
    // We want to write the data out in parallel.
    //
    // We'll allow up to nOutFiles active writers at a time.
    //
    const int nOutFiles = std::min(64,NProcs);

    for (int lev = 0; lev <= m_amr->finestLevel(); lev++)
    {
        //
        // We store the particles at each level in their own subdirectory.
        //
        std::string LevelDir = pdir;

        if (!LevelDir.empty() && LevelDir[LevelDir.size()-1] != '/')
            LevelDir += '/';
        
        char buf[32];
        sprintf(buf, "Level_%d", lev);

        LevelDir += buf;

        if (ParallelDescriptor::IOProcessor())
            if (!BoxLib::UtilCreateDirectory(LevelDir, 0755))
                BoxLib::CreateDirectoryFailed(LevelDir);
        //
        // Force other processors to wait till directory is built.
        //
        ParallelDescriptor::Barrier();

        const MultiFab& state = m_amr->getLevel(lev).get_new_data(0);
        //
        // We eventually want to write out the file name and the offset
        // into that file into which each grid of particles is written.
        //
        Array<int>  which(state.size(),0);
        Array<int > count(state.size(),0);
        Array<long> where(state.size(),0);

        const int FileNumber = MyProc % nOutFiles;

        sprintf(buf, "%04d", FileNumber);

        std::string FullFileName = LevelDir;

        FullFileName += '/';
        FullFileName += ParticleBase::DataPrefix();
        FullFileName += buf;

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
                    BoxLib::FileOpenFailed(FullFileName);
                //
                // Write out all the particles we own at the specified level.
                // Do it grid block by grid block remembering the seek offset
                // for the start of writing of each block of data.
                //
                WriteParticles(lev, ParticleFile, FileNumber, which, count, where);

                ParticleFile.flush();

#ifdef BL_USECLOSE
                ParticleFile.close();
#endif
                if (!ParticleFile.good())
                    BoxLib::Abort("ParticleContainer<N>::Checkpoint(): problem writing ParticleFile");

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
            //
            // Unlink any zero-length data files.
            //
            Array<long> cnt(nOutFiles,0);

            for (int i = 0; i < count.size(); i++)
                cnt[which[i]] += count[i];

            for (int i = 0; i < cnt.size(); i++)
            {
                if (cnt[i] == 0)
                {
                    char buf[32];

                    sprintf(buf, "%04d", i);

                    std::string FullFileName = LevelDir;

                    FullFileName += '/';
                    FullFileName += ParticleBase::DataPrefix();
                    FullFileName += buf;

                    BoxLib::UnlinkFile(FullFileName.c_str());
                }
            }
        }
    }

    Real stoptime = ParallelDescriptor::second() - strttime;

    ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        HdrFile.flush();

#ifdef BL_USECLOSE
        HdrFile.close();
#endif
        if (!HdrFile.good())
            BoxLib::Abort("ParticleContainer<N>::Checkpoint(): problem writing HdrFile");

        std::cout << "\nParticleContainer<N>::Checkpoint() time: " << stoptime << "\n\n";
    }
}

template <int N>
void
ParticleContainer<N>::WriteParticles (int            lev,
                                      std::ofstream& ofs,
                                      int            fnum,
                                      Array<int>&    which,
                                      Array<int>&    count,
                                      Array<long>&   where) const
{
    BL_ASSERT(m_amr->finestLevel() <= lev+1);

    const MultiFab& state = m_amr->getLevel(lev).get_new_data(0);

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const int grid = mfi.index();

        int cnt = 0;

        for (const_iterator it = m_particles.begin(); it != m_particles.end(); ++it)
            if (it->m_lev == lev && it->m_grid == grid)
                cnt++;

        which[grid] = fnum;
        count[grid] = cnt;
        where[grid] = VisMF::FileOffset(ofs);

        if (cnt == 0) continue;

        {
            //
            // First write out the integer data in binary.
            //
            const int iChunkSize = 4+BL_SPACEDIM;

            Array<int> istuff(cnt*iChunkSize);

            int* iptr = istuff.dataPtr();

            for (const_iterator it = m_particles.begin(); it != m_particles.end(); ++it)
            {
                const ParticleType& p = *it;

                if (p.m_lev == lev && p.m_grid == grid)
                {
                    iptr[0] = p.m_id;
                    iptr[1] = p.m_cpu;
                    iptr[2] = p.m_lev;
                    iptr[3] = p.m_grid;

                    for (int i = 0; i < BL_SPACEDIM; i++)
                        iptr[4+i] = p.m_cell[i];

                    iptr += iChunkSize;
                }
            }

            ofs.write((char*)istuff.dataPtr(),istuff.size()*sizeof(int));
        }

        {
            //
            // Then the Real data in binary.
            //
            const int rChunkSize = BL_SPACEDIM+N;

            Array<Real> rstuff(cnt*rChunkSize);

            Real* rptr = rstuff.dataPtr();

            for (const_iterator it = m_particles.begin(); it != m_particles.end(); ++it)
            {
                const ParticleType& p = *it;

                if (p.m_lev == lev && p.m_grid == grid)
                {
                    for (int i = 0; i < BL_SPACEDIM; i++)
                        rptr[i] = p.m_pos[i];

                    for (int i = 0; i < N; i++)
                        rptr[BL_SPACEDIM+i] = p.m_data[i];

                    rptr += rChunkSize;
                }
            }

            ofs.write((char*)rstuff.dataPtr(),rstuff.size()*sizeof(Real));
        }
    }
}

template <int N>
void
ParticleContainer<N>::Restart (const std::string& dir,
                               const std::string& file)
{
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
        BoxLib::FileOpenFailed(HdrFileName);
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

        for (int i = 0; i < version.size(); i++)
            vbuf[i] = version[i];

        vbuf[version.size()] = '\0';
    }
    ParallelDescriptor::Bcast(vbuf.dataPtr(), vbuf.size(), IOProc);

    version = vbuf.dataPtr();

    if (version == ParticleBase::Version())
    {
        Restart_OneDotZero(fullname,HdrFile);
    }
    else
    {
        std::string msg("ParticleContainer<N>::Restart(): unknown version string: ");
        msg += version;
        BoxLib::Abort(msg.c_str());
    }

    Real stoptime = ParallelDescriptor::second() - strttime;

    ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "\nParticleContainer<N>::Restart() time: " << stoptime << "\n\n";
    }
}

template <int N>
void
ParticleContainer<N>::Restart_OneDotZero (const std::string& fullname,
                                          std::ifstream&     HdrFile)
{
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
            BoxLib::Abort("ParticleContainer<N>::Restart(): dm != BL_SPACEDIM");
    }
    ParallelDescriptor::Bcast(&dm, 1, IOProc);
    //
    // Next value should be our "N".
    //
    int n;

    if (ParallelDescriptor::IOProcessor())
    {
        HdrFile >> n;

        if (n != N)
            BoxLib::Abort("ParticleContainer<N>::Restart(): n != N");
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

    BL_ASSERT(finest_level == m_amr->finestLevel());

    if (ParallelDescriptor::IOProcessor())
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            HdrFile >> ngrids[lev];

            BL_ASSERT(ngrids[lev] > 0);
            BL_ASSERT(ngrids[lev] == m_amr->boxArray(lev).size());
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

        const MultiFab& state = m_amr->getLevel(lev).get_new_data(0);

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

            char buf[64];

            sprintf(buf,
                    "Level_%d/%s%04d",
                    lev,
                    ParticleBase::DataPrefix().c_str(),
                    which[grid]);

            name += buf;

            std::ifstream ParticleFile;

            ParticleFile.open(name.c_str(), std::ios::in);

            if (!ParticleFile.good())
                BoxLib::FileOpenFailed(name);

            ParticleFile.seekg(where[grid], std::ios::beg);

            ReadParticles_OneDotZero(count[grid],grid,lev,ParticleFile);

#ifdef BL_USECLOSE
            ParticleFile.close();
#endif
            if (!ParticleFile.good())
                BoxLib::Abort("ParticleContainer<N>::Restart(): problem reading particles");
        }
    }

    BL_ASSERT(OK());
}

template <int N>
void
ParticleContainer<N>::ReadParticles_OneDotZero (int            cnt,
                                                int            grd,
                                                int            lev,
                                                std::ifstream& ifs)
{
    BL_ASSERT(cnt > 0);
    BL_ASSERT(lev >= 0 && lev <= m_amr->finestLevel());
    BL_ASSERT(grd >= 0 && grd < m_amr->boxArray(lev).size());

    ParticleType p;
    //
    // First read in the integer data in binary.
    //
    const int iChunkSize = 4+BL_SPACEDIM;

    Array<int> istuff(cnt*iChunkSize);

    ifs.read((char*)istuff.dataPtr(),istuff.size()*sizeof(int));
    //
    // Then the Real data in binary.
    //
    const int rChunkSize = BL_SPACEDIM+N;

    Array<Real> rstuff(cnt*rChunkSize);

    ifs.read((char*)rstuff.dataPtr(),rstuff.size()*sizeof(Real));
    //
    // Now reassemble the particles.
    //
    int*  iptr = istuff.dataPtr();
    Real* rptr = rstuff.dataPtr();

    for (int i = 0; i < cnt; i++)
    {
        p.m_id   = iptr[0];
        p.m_cpu  = iptr[1];
        p.m_lev  = iptr[2];
        p.m_grid = iptr[3];

        BL_ASSERT(p.m_lev  == lev);
        BL_ASSERT(p.m_grid == grd);

        for (int i = 0; i < BL_SPACEDIM; i++)
            p.m_cell[i] = iptr[4+i];

        iptr += iChunkSize;

        for (int i = 0; i < BL_SPACEDIM; i++)
            p.m_pos[i] = rptr[i];

        for (int i = 0; i < N; i++)
            p.m_data[i] = rptr[BL_SPACEDIM+i];

        rptr += rChunkSize;

        m_particles.push_back(p);
    }
}

#endif /*_PARTICLES_H_*/
