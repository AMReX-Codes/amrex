#ifndef _PARTICLES_H_
#define _PARTICLES_H_ 

#include <list>
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
    typedef typename std::list<ParticleType>::iterator iterator;

    typedef typename std::list<ParticleType>::const_iterator const_iterator;

    ParticleContainer (Amr* amr)
        :
        m_verbose(1), m_amr(amr) { BL_ASSERT(amr != 0); }

    void Init (const std::string& file);

    void InitRandom (long icount, unsigned long iseed, Real particleMass);

    Real sumParticleMass (const MultiFab& mf, int level) const;

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
    void movePredict(const MultiFab& grav_vector, int level, Real timestep);
    void moveCorrect(const MultiFab& grav_vector_old, const MultiFab& grav_vector, int level, Real timestep);
    //
    // TODO: the methods should return a constraint on the timestep...
    //
    // The following two functions form a KICK DRIFT KICK scheme for integrating the motion of the particles
    //
    void moveKickDrift(const MultiFab& grav_vector, int level, Real timestep);
    void moveKick     (const MultiFab& grav_vector, int level, Real timestep);
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
                         Array<long>&   where) const;
    //
    // Helper function for Restart().
    //
    void ReadParticles (int            cnt,
                        int            grd,
                        int            lev,
                        std::ifstream& ifs);
    //
    // The data.
    //
    int                     m_verbose;
    Amr*                    m_amr;
    std::list<ParticleType> m_particles;
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
    const int MyProc = ParallelDescriptor::MyProc();

    std::ifstream ifs;

    ifs.open(file.c_str(), std::ios::in);

    if (!ifs.good())
        BoxLib::FileOpenFailed(file);

    int id = 1;

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

    const int NProcs = ParallelDescriptor::NProcs();

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
ParticleContainer<N>::InitRandom (long icount, unsigned long iseed, Real particleMass)
{
    BL_ASSERT(iseed  > 0);
    BL_ASSERT(icount > 0);

    BL_ASSERT(m_amr != 0);

    const int  MyProc = ParallelDescriptor::MyProc();
    const int  NProcs = ParallelDescriptor::NProcs();
    const int  IOProc = ParallelDescriptor::IOProcessorNumber();
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

        p.m_data[0] = particleMass;
        for (int i = 1; i < N; i++)
            //
            // Just zero out the data for lack of a better value.
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

    const int  MyProc = ParallelDescriptor::MyProc();
    const int  NProcs = ParallelDescriptor::NProcs();

    for (iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        ParticleType& p = *it;

        BL_ASSERT(p.m_id > 0);

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
ParticleContainer<N>::Increment (MultiFab& mf, int level) const
{
    const int MyProc = ParallelDescriptor::MyProc();

    const DistributionMapping& dm = mf.DistributionMap();

    for (const_iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        const ParticleType& p = *it;

        if (p.m_lev == level && dm[p.m_grid] == MyProc)
        {
            mf[p.m_grid](p.m_cell) += 1;
        }
    }
}

template <int N>
Real
ParticleContainer<N>::estTimestep (const MultiFab& grav_vector, int level) const
{
    BL_ASSERT(N >= BL_SPACEDIM+1);

    const int MyProc = ParallelDescriptor::MyProc();

    const DistributionMapping& dm = grav_vector.DistributionMap();

    Real dt = 1e200, mag_vel, mag_grav, dt_part;

    const Geometry& geom = m_amr->Geom(0);
    const Real*     dx   = geom.CellSize();

    for (const_iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        const ParticleType& p = *it;

        if (p.m_lev == level && dm[p.m_grid] == MyProc)
        {
            D_TERM(mag_vel  = p.m_data[1] * p.m_data[1];,
                   mag_vel += p.m_data[2] * p.m_data[2];,
                   mag_vel += p.m_data[3] * p.m_data[3];);

            mag_vel = sqrt(mag_vel);

            dt_part = 1e200;

            if (mag_vel > 0)
                dt_part = 0.1 * dx[0] / mag_vel;

            const FArrayBox& gfab = grav_vector[p.m_grid];

            Real gval[BL_SPACEDIM] = { D_DECL(gfab(p.m_cell,0), gfab(p.m_cell,1), gfab(p.m_cell,2)) };

            D_TERM(mag_grav  = gval[0]*gval[0];,
                   mag_grav += gval[1]*gval[1];,
                   mag_grav += gval[2]*gval[2];);

            mag_grav = sqrt(mag_grav);

            if (mag_grav > 0)
                dt_part = std::min( dt_part, 1./sqrt( mag_grav / dx[0] ) );

            dt = std::min(dt_part, dt); 
        }
    }

    ParallelDescriptor::ReduceRealMin(dt);

    return dt;
}

//
// Assumes mass is in m_data[0]!
//

template <int N>
Real
ParticleContainer<N>::sumParticleMass (const MultiFab& mf, int level) const
{
    BL_ASSERT(N >= 1);

    const int MyProc = ParallelDescriptor::MyProc();

    const DistributionMapping& dm = mf.DistributionMap();

    Real massSum = 0.0;

    for (const_iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        const ParticleType& p = *it;

        if (p.m_lev == level && dm[p.m_grid] == MyProc)
        {
            massSum += p.m_data[0];
        }
    }

    ParallelDescriptor::ReduceRealSum(massSum);

    return massSum;
}

//
// Some "pseudo" cic-scheme (as in enzo)
// assumes grid has the size of 8x8 ;) -> here be dragons
// Will get ugly as soon as there is more than just the rootgrid
// NOTE: does NOT return density but mass assigned in a proper way!
// might once be ported to fortran.
//

template <int N>
void
ParticleContainer<N>::AssignDensity (MultiFab& mf, int level) const
{
    BL_ASSERT(N >= 1);

    const int MyProc = ParallelDescriptor::MyProc();

    const DistributionMapping& dm = mf.DistributionMap();

    const Geometry& geom = m_amr->Geom(level);

    if (!geom.isAllPeriodic())
        BoxLib::Error("AssignDensity: Problem has to be periodic...");

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        mf[mfi].setVal(0);
    }

    //what the cave? does this really work?
    //grow box, it shrinks again when the ghost cells are summed in the end.
    //    BoxLib::grow(mf.boxArray()[0],2);
    //    for (MFIter mfi(mf); mfi.isValid(); ++mfi) 
    //    { 
    //       Box bx = BoxLib::grow(mfi.validbox(),1);
    //    }

    const Real* dx = geom.CellSize();
    Real fracOfCell[BL_SPACEDIM], cell_vol;

    D_TERM(cell_vol  = dx[0];,
           cell_vol *= dx[1];,
           cell_vol *= dx[2];);

    for (const_iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        const ParticleType& p = *it;

        if (p.m_lev == level && dm[p.m_grid] == MyProc)
        {
            IntVect centralIntersection;

//FIXME: m_pos seems to be relative to geom.ProbLo?!?
            D_TERM(centralIntersection [0] = floor( p.m_pos[0]/dx[0] + 0.5);,
                   centralIntersection [1] = floor( p.m_pos[1]/dx[1] + 0.5);,
                   centralIntersection [2] = floor( p.m_pos[2]/dx[2] + 0.5););

            IntVect cell = centralIntersection;
            
            D_TERM(fracOfCell[0] = -centralIntersection[0] + p.m_pos[0]/dx[0] + 0.5;,
                   fracOfCell[1] = -centralIntersection[1] + p.m_pos[1]/dx[1] + 0.5;,
                   fracOfCell[2] = -centralIntersection[2] + p.m_pos[2]/dx[2] + 0.5;);

            FArrayBox& fab = mf[p.m_grid];

#if (BL_SPACEDIM == 1)
            fab(centralIntersection) += p.m_data[0] *    fracOfCell[0];
    
            cell[0] = centralIntersection[0] - 1;
            fab(cell)                += p.m_data[0] * (1-fracOfCell[0]);
#endif
            
#if (BL_SPACEDIM == 2)
            fab(centralIntersection) += p.m_data[0] *    fracOfCell[0]  *    fracOfCell[1] ;
    
            cell[0] = centralIntersection[0] - 1;
            fab(cell)                += p.m_data[0] * (1-fracOfCell[0]) *    fracOfCell[1] ;
    
            cell = centralIntersection;
            cell[1] = centralIntersection[1] - 1;
            fab(cell)                += p.m_data[0] *    fracOfCell[0]  * (1-fracOfCell[1]);
    
            cell[0] = centralIntersection[0] - 1;
            fab(cell)                += p.m_data[0] * (1-fracOfCell[0]) * (1-fracOfCell[1]);
#endif
 

#if (BL_SPACEDIM == 3)
            fab(centralIntersection) += p.m_data[0] *    fracOfCell[0]  *    fracOfCell[1]  *    fracOfCell[2] ;
    
            cell[0] = centralIntersection[0] - 1;
            fab(cell)                += p.m_data[0] * (1-fracOfCell[0]) *    fracOfCell[1]  *    fracOfCell[2] ;
    
            cell = centralIntersection;
            cell[1] = centralIntersection[1] - 1;
            fab(cell)                += p.m_data[0] *    fracOfCell[0]  * (1-fracOfCell[1]) *    fracOfCell[2] ;
    
            cell[0] = centralIntersection[0] - 1;
            fab(cell)                += p.m_data[0] * (1-fracOfCell[0]) * (1-fracOfCell[1]) *    fracOfCell[2] ;

            
            centralIntersection[3] = centralIntersection[3] - 1; //ugly hack - but this way I can "reuse" the code fr2m 2d...
            fab(centralIntersection) += p.m_data[0] *    fracOfCell[0]  *    fracOfCell[1]  * (1-fracOfCell[2]);
    
            cell[0] = centralIntersection[0] - 1;
            fab(cell)                += p.m_data[0] * (1-fracOfCell[0]) *    fracOfCell[1]  * (1-fracOfCell[2]);
    
            cell = centralIntersection;
            cell[1] = centralIntersection[1] - 1;
            fab(cell)                += p.m_data[0] *    fracOfCell[0]  * (1-fracOfCell[1]) * (1-fracOfCell[2]);
    
            cell[0] = centralIntersection[0] - 1;
            fab(cell)                += p.m_data[0] * (1-fracOfCell[0]) * (1-fracOfCell[1]) * (1-fracOfCell[2]);
#endif
        }
    }
    mf.SumBoundary();
    geom.SumPeriodicBoundary(mf);
    mf.mult(1./cell_vol);
}

template <int N>
void
ParticleContainer<N>::movePredict(const MultiFab& grav_vector, int level, Real timestep)
{
    BL_ASSERT(N >= BL_SPACEDIM+1);

    const int MyProc = ParallelDescriptor::MyProc();

    const DistributionMapping& dm = grav_vector.DistributionMap();

    for (iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        ParticleType& p = *it;

        if (p.m_lev == level && dm[p.m_grid] == MyProc)
        {
            //
            // Note: m_data[0] is mass, 1 is v_x, ...
            //
            const FArrayBox& gfab = grav_vector[p.m_grid];

            D_TERM(p.m_data[1] += timestep * gfab(p.m_cell,0);,
                   p.m_data[2] += timestep * gfab(p.m_cell,1);,
                   p.m_data[3] += timestep * gfab(p.m_cell,2););

            D_TERM(p.m_pos[0]  += timestep * p.m_data[1];,
                   p.m_pos[1]  += timestep * p.m_data[2];,
                   p.m_pos[2]  += timestep * p.m_data[3];);
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
        std::cout << "Calling Redistribute() from ParticleContainer<N>::MovePredict () ...\n";

    Redistribute();
}

template <int N>
void
ParticleContainer<N>::moveCorrect(const MultiFab& grav_vector_old, const MultiFab& grav_vector, int level, Real timestep)
{
    BL_ASSERT(N >= BL_SPACEDIM+1);

    const int MyProc = ParallelDescriptor::MyProc();

    const DistributionMapping& dm = grav_vector.DistributionMap();

    for (iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        ParticleType& p = *it;

        if (p.m_lev == level && dm[p.m_grid] == MyProc)
        {
            D_TERM(p.m_pos[0]  -= 0.5 * timestep * p.m_data[1];,
                   p.m_pos[1]  -= 0.5 * timestep * p.m_data[2];,
                   p.m_pos[2]  -= 0.5 * timestep * p.m_data[3];);
            //
            // Note: m_data[0] is mass, 1 is v_x, ...
            //
            const FArrayBox& gfab     = grav_vector[p.m_grid];
            const FArrayBox& gfab_old = grav_vector_old[p.m_grid];
            
            D_TERM(p.m_data[1] += 0.5 * timestep * ( gfab(p.m_cell,0) - gfab_old(p.m_cell,0) );, 
                   p.m_data[2] += 0.5 * timestep * ( gfab(p.m_cell,1) - gfab_old(p.m_cell,1) );,
                   p.m_data[3] += 0.5 * timestep * ( gfab(p.m_cell,2) - gfab_old(p.m_cell,2) ););

            D_TERM(p.m_pos[0]  += 0.5 * timestep * p.m_data[1];,
                   p.m_pos[1]  += 0.5 * timestep * p.m_data[2];,
                   p.m_pos[2]  += 0.5 * timestep * p.m_data[3];);
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
        std::cout << "Calling Redistribute() from ParticleContainer<N>::MoveCorrect () ...\n";

    Redistribute();
}

template <int N>
void
ParticleContainer<N>::moveKickDrift(const MultiFab& grav_vector, int level, Real timestep)
{
    BL_ASSERT(N >= BL_SPACEDIM+1);

    const int MyProc = ParallelDescriptor::MyProc();

    const DistributionMapping& dm = grav_vector.DistributionMap();

    for (iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        ParticleType& p = *it;

        if (p.m_lev == level && dm[p.m_grid] == MyProc)
        {
            //
            // note: m_data[0] is mass, 1 is v_x, ...
            //
            const FArrayBox& gfab = grav_vector[p.m_grid];

            D_TERM(p.m_data[1] += 0.5 * timestep * gfab(p.m_cell,0);,
                   p.m_data[2] += 0.5 * timestep * gfab(p.m_cell,1);,
                   p.m_data[3] += 0.5 * timestep * gfab(p.m_cell,2););

            D_TERM(p.m_pos[0]  +=       timestep * p.m_data[1];,
                   p.m_pos[1]  +=       timestep * p.m_data[2];,
                   p.m_pos[2]  +=       timestep * p.m_data[3];);
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
            // Modulo any arithmetic funny business the particle should in the problem domain.
            //
        }
    }

    if (m_verbose > 1 && ParallelDescriptor::IOProcessor())
        std::cout << "Calling Redistribute() from ParticleContainer<N>::moveKickDrift () ...\n";

    Redistribute();
}

template <int N>
void
ParticleContainer<N>::moveKick(const MultiFab& grav_vector, int level, Real timestep)
{
    BL_ASSERT(N >= BL_SPACEDIM+1);

    const int MyProc = ParallelDescriptor::MyProc();

    const DistributionMapping& dm = grav_vector.DistributionMap();

    for (iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        ParticleType& p = *it;

        if (p.m_lev == level && dm[p.m_grid] == MyProc)
        {
            //
            // Note: m_data[0] is mass, 1 is v_x, ...
            //
            const FArrayBox& gfab = grav_vector[p.m_grid];

            D_TERM(p.m_data[1] += 0.5 * timestep * gfab(p.m_cell,0);, //how do I use grav_vec correctly? this seems to be some kind of s...
                   p.m_data[2] += 0.5 * timestep * gfab(p.m_cell,1);,
                   p.m_data[3] += 0.5 * timestep * gfab(p.m_cell,2););
        }
    }
    //
    // No need for redistribution, we only change the velocity.
    //
}

template <int N>
void
ParticleContainer<N>::Redistribute (bool where_already_called)
{
    const int  MyProc   = ParallelDescriptor::MyProc();
    const Real strttime = ParallelDescriptor::second();

    std::vector<int>      owner;
    std::vector<iterator> which;

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
            which.push_back(it);
            owner.push_back(who);
        }
    }

    const int NProcs = ParallelDescriptor::NProcs();

    if (NProcs == 1)
    {
        BL_ASSERT(which.empty());
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
        const int ChunkSize = BL_SPACEDIM+4;

        Array<int> recvdata (NumRcvs*ChunkSize);
        Array<int> senddata (NumSnds*ChunkSize);

        Array<int> sendcnts(NProcs,0), sdispls(NProcs,0);
        Array<int> recvcnts(NProcs,0), rdispls(NProcs,0), offset(NProcs,0);

        for (int i = 0; i < NProcs; i++)
        {
            recvcnts[i] = Rcvs[i] * ChunkSize;
            sendcnts[i] = Snds[i] * ChunkSize;
        }

        for (int i = 1; i < NProcs; i++)
        {
            offset [i] = offset [i-1] + sendcnts[i-1];
            rdispls[i] = rdispls[i-1] + recvcnts[i-1];
            sdispls[i] = sdispls[i-1] + sendcnts[i-1];
        }

        for (int i = 0; i < which.size(); i++)
        {
            int idx = offset[owner[i]];

            ParticleType& p = *which[i];

            senddata[idx+0] = p.m_id;
            senddata[idx+1] = p.m_cpu;
            senddata[idx+2] = p.m_lev;
            senddata[idx+3] = p.m_grid;

            idx += 4;

            D_TERM(senddata[idx+0]=p.m_cell[0];,
                   senddata[idx+1]=p.m_cell[1];,
                   senddata[idx+2]=p.m_cell[2];);

            offset[owner[i]] += ChunkSize;
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
        const int ChunkSize = BL_SPACEDIM+N;

        Array<Real> recvdata (NumRcvs*ChunkSize);
        Array<Real> senddata (NumSnds*ChunkSize);

        Array<int> sendcnts(NProcs,0), sdispls(NProcs,0);
        Array<int> recvcnts(NProcs,0), rdispls(NProcs,0), offset(NProcs,0);

        for (int i = 0; i < NProcs; i++)
        {
            recvcnts[i] = Rcvs[i] * ChunkSize;
            sendcnts[i] = Snds[i] * ChunkSize;
        }

        for (int i = 1; i < NProcs; i++)
        {
            offset [i] = offset [i-1] + sendcnts[i-1];
            rdispls[i] = rdispls[i-1] + recvcnts[i-1];
            sdispls[i] = sdispls[i-1] + sendcnts[i-1];
        }

        for (int i = 0; i < which.size(); i++)
        {
            int idx = offset[owner[i]];

            ParticleType& p = *which[i];

            D_TERM(senddata[idx+0]=p.m_pos[0];,
                   senddata[idx+1]=p.m_pos[1];,
                   senddata[idx+2]=p.m_pos[2];);

            idx += BL_SPACEDIM;

            for (int j = 0; j < N; j++)
                senddata[idx+j] = p.m_data[j];

            offset[owner[i]] += ChunkSize;
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
        for (int i = 0; i < which.size(); i++)
        {
            m_particles.erase(which[i]);
        }
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

    if (m_verbose > 0 && NProcs > 1)
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
ParticleContainer<N>::Checkpoint (const std::string& dir, const std::string& name) const
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
    // The actual particle data is stored in files of the form: DATA_nnnn.
    //
    static const std::string DataPrefix = "DATA_";
    //
    // We want to write the data in parallel.
    // We'll allow up to nOutFiles active writers at a time.
    //
    int nOutFiles = std::min(64,NProcs);

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
        Array<long> where(state.size(),0);

        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

        std::ofstream ParticleFile;

        ParticleFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        const int FileNumber = MyProc % nOutFiles;

        sprintf(buf, "%04d", FileNumber);

        std::string FullFileName = LevelDir;

        FullFileName += '/';
        FullFileName += DataPrefix;
        FullFileName += buf;

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
                {
                    //
                    // First set.
                    //
                    ParticleFile.open(FullFileName.c_str(),
                                      std::ios::out|std::ios::trunc|std::ios::binary);
                }
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
                WriteParticles(lev, ParticleFile, FileNumber, which, where);

                ParticleFile.flush();
#ifdef BL_USECLOSE
                ParticleFile.close();
#endif
                if (!ParticleFile.good())
                    BoxLib::Abort("ParticleContainer<N>::Checkpoint(): problem writing ParticleFile");

                int iBuff     = 0;
                int wakeUpPID = (MyProc + nOutFiles);
                int tag       = (MyProc % nOutFiles);

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
                int iBuff;
                int waitForPID = (MyProc - nOutFiles);
                int tag        = (MyProc % nOutFiles);
                ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
            }
        }

        ParallelDescriptor::ReduceIntSum(which.dataPtr(), which.size(), IOProc);

        ParallelDescriptor::ReduceLongSum(where.dataPtr(), where.size(), IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            char buf[64];

            for (int j = 0; j < state.size(); j++)
            {
                //
                // We now write the file name & file offset into which
                // the data for each grid was written to the header file.
                // We only write the filename relative to the main
                // directory containing the particles.
                //
                sprintf(buf,
                        "Level_%d/%s%04d",
                        lev,
                        DataPrefix.c_str(),
                        which[j]);

                HdrFile << buf << ' ' << where[j] << '\n';
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
                                      Array<long>&   where) const
{
    BL_ASSERT(m_amr->finestLevel() <= lev);

    const MultiFab& state = m_amr->getLevel(lev).get_new_data(0);

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const int grid = mfi.index();

        which[grid] = fnum;
        where[grid] = VisMF::FileOffset(ofs);

        int cnt = 0;

        for (const_iterator it = m_particles.begin(); it != m_particles.end(); ++it)
            if (it->m_lev == lev && it->m_grid == grid)
                cnt++;
        
        ofs << cnt << '\n';

        if (cnt == 0) continue;

        {
            //
            // First write out the integer data in binary.
            //
            const int ChunkSize = 4+BL_SPACEDIM;

            Array<int> istuff(cnt*ChunkSize);

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

                    iptr += ChunkSize;
                }
            }

            ofs.write((char*)istuff.dataPtr(),istuff.size()*sizeof(int));
        }

        {
            //
            // Then the Real data in binary.
            //
            const int ChunkSize = BL_SPACEDIM+N;

            Array<Real> rstuff(cnt*ChunkSize);

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

                    rptr += ChunkSize;
                }
            }

            ofs.write((char*)rstuff.dataPtr(),rstuff.size()*sizeof(Real));
        }
    }
}

template <int N>
void
ParticleContainer<N>::Restart (const std::string& dir,const std::string& file)
{
    const int  IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real strttime = ParallelDescriptor::second();

    std::string fullname = dir;

    if (!fullname.empty() && fullname[fullname.size()-1] != '/')
        fullname += '/';

    fullname += file;
    //
    // The header contains the info we need to read back in the particles.
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
    // First value should be BL_SPACEDIM;
    //
    int dm;
        
    HdrFile >> dm;

    if (dm != BL_SPACEDIM)
        BoxLib::Abort("ParticleContainer<N>::Restart(): dm != BL_SPACEDIM");
    //
    // Next value should be our "N".
    //
    int n;

    HdrFile >> n;

    if (n != N)
        BoxLib::Abort("ParticleContainer<N>::Restart(): n != N");

    long nparticles;
    //
    // The total number of particles.
    //
    HdrFile >> nparticles;

    BL_ASSERT(nparticles >= 0);
    //
    // Then the finest level of the AMR hierarchy.
    //
    int finest_level;

    HdrFile >> finest_level;

    BL_ASSERT(finest_level >= 0);
    //
    // Then the number of grids at each level.
    //
    Array<int> ngrids(finest_level+1);

    BL_ASSERT(finest_level == m_amr->finestLevel());

    for (int lev = 0; lev <= finest_level; lev++)
    {
        HdrFile >> ngrids[lev];

        BL_ASSERT(ngrids[lev] > 0);
        BL_ASSERT(ngrids[lev] == m_amr->boxArray(lev).size());
    }
    //
    // The rest of HdrFile consists of pairs of the form:
    //
    //   filename offset
    //
    // One for each grid at each level from 0 -> finest_level.
    //
    for (int lev = 0; lev <= finest_level; lev++)
    {
        typedef std::pair< std::string, long > line;
        //
        // Read in the grid file/offset info for this level.
        //
        Array<line> lines(ngrids[lev]);

        for (int i = 0; i < ngrids[lev]; i++)
        {
            HdrFile >> lines[i].first >> lines[i].second;
        }

        const MultiFab& state = m_amr->getLevel(lev).get_new_data(0);

        for (MFIter mfi(state); mfi.isValid(); ++mfi)
        {
            const int grid = mfi.index();
            //
            // The file names in the header file are relative.
            //
            std::string name = fullname;

            if (!name.empty() && name[name.size()-1] != '/')
                name += '/';

            name += lines[grid].first;

            std::ifstream ParticleFile;

            ParticleFile.open(name.c_str(), std::ios::in);

            if (!ParticleFile.good())
                BoxLib::FileOpenFailed(name);

            ParticleFile.seekg(lines[grid].second, std::ios::beg);
            //
            // First thing is an ASCII line of the form:
            //
            // cnt\n
            //
            int cnt;

            ParticleFile >> cnt;
            //
            // Skip past the newline
            //
            int ignore = ParticleFile.get();

            if (cnt > 0)
            {
                ReadParticles(cnt,grid,lev,ParticleFile);

                if (!ParticleFile.good())
                    BoxLib::Abort("ParticleContainer<N>::Restart(): problem reading particles");
            }

#ifdef BL_USECLOSE
            ParticleFile.close();
#endif
        }
    }

    BL_ASSERT(OK());

#ifndef NDEBUG
    long sz = m_particles.size();

    ParallelDescriptor::ReduceLongSum(sz,IOProc);

    if (ParallelDescriptor::IOProcessor())
        BL_ASSERT(sz == nparticles);
#endif

    Real stoptime = ParallelDescriptor::second() - strttime;

    ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "\nParticleContainer<N>::Restart() time: " << stoptime << "\n\n";
    }
}

template <int N>
void
ParticleContainer<N>::ReadParticles (int            cnt,
                                     int            grd,
                                     int            lev,
                                     std::ifstream& ifs)
{
    BL_ASSERT(cnt > 0);
    BL_ASSERT(lev >= 0 && lev <= m_amr->finestLevel());
    BL_ASSERT(grd >= 0 && grd < m_amr->boxArray(lev).size());

    ParticleType p;

    std::list<ParticleType> newparticles;
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
