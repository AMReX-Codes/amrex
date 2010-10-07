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

    void InitRandom (long icount, unsigned long iseed);

    Real sumParticleMass (MultiFab& mf, int level) const;

    void Increment (MultiFab& mf, int level) const;

    void AssignDensity (MultiFab& mf, int level) const;

    void Redistribute (bool where_already_called = false);

    bool OK () const;

    void ByteSpread () const;

    void MoveRandom ();

    // **************************************************************************************************************** 

    // If the particles move only with self-gravity from themselves and the gas, then 
    //   we can move them according to the schemes below.
    // The gravitational force must be computed between the calls of the parts of the integration scheme.

    // The following two functions form a PREDICTOR CORRECTOR scheme for integrating the motion of the particles
    void movePredict(MultiFab& grav_vector, int level, Real timestep);
    void moveCorrect(MultiFab& grav_vector_old, MultiFab& grav_vector, int level, Real timestep);
    //TODO: the methods should return a constraint on the timestep...

    // The following two functions form a KICK DRIFT KICK scheme for integrating the motion of the particles
    void moveKickDrift(MultiFab& grav_vector, int level, Real timestep);
    void moveKick( MultiFab& grav_vector, int level, Real timestep);
    // after the moveKickDrift step the positions of the particles are advanced for a full timestep,
    // so this scheme should work in the overall algorithm...

    // **************************************************************************************************************** 

    int Verbose () { return m_verbose; }

    void Verbose (int verbose) { m_verbose = verbose; }

protected:
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
                p.m_id  = m_particles.size() + 1;
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
ParticleContainer<N>::InitRandom (long icount, unsigned long iseed)
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
        
        p.m_data[0] = 1.;
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

    if (ParallelDescriptor::IOProcessor())
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

        if (p.m_id > 0 && p.m_lev == level && dm[p.m_grid] == MyProc)
        {
            mf[p.m_grid](p.m_cell) += 1;
        }
    }
}

//assumes mass is in m_data[0]!
//FIXME: might not work with MPI!!!
template <int N>
Real
ParticleContainer<N>::sumParticleMass (MultiFab& mf, int level) const
{
    const int MyProc = ParallelDescriptor::MyProc();

    const DistributionMapping& dm = mf.DistributionMap();

    Real massSum = 0.0;

    for (const_iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        const ParticleType& p = *it;

        if (p.m_id > 0 && p.m_lev == level && dm[p.m_grid] == MyProc)
        {
            massSum += p.m_data[0];
        }
    }
    return massSum;
}

//some "pseudo" cic-scheme (as in enzo)
//assumes grid has the size of 8x8 ;) -> here be dragons
//Will get ugly as soon as there is more than just the rootgrid
//NOTE: does NOT return density but mass assigned in a proper way!
//might once be ported to fortran
template <int N>
void
ParticleContainer<N>::AssignDensity (MultiFab& mf, int level) const
{
    const int MyProc = ParallelDescriptor::MyProc();

    const DistributionMapping& dm = mf.DistributionMap();

    const Geometry& geom = m_amr->Geom(level);
    if(!geom.isAllPeriodic()){
        BoxLib::Error("AssignDensity: Problem has to be periodic...");
    }
    mf.setVal(0.0,0,1,1);

//what the cave? does this really work?
//grow box, it shrinks again when the ghost cells are summed in the end.
//    BoxLib::grow(mf.boxArray()[0],2);
//    for (MFIter mfi(mf); mfi.isValid(); ++mfi) 
//    { 
//       Box bx = BoxLib::grow(mfi.validbox(),1);
//    }

    const Real* dx = geom.CellSize();
    Real fracOfCell[BL_SPACEDIM];
    Real cell_vol;
    D_TERM(cell_vol  = dx[0];,
           cell_vol *= dx[1];,
           cell_vol *= dx[2];);

    for (const_iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        const ParticleType& p = *it;


        if (p.m_id > 0 && p.m_lev == level && dm[p.m_grid] == MyProc)
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

#if (BL_SPACEDIM == 1)
            mf[p.m_grid](centralIntersection) += p.m_data[0] *    fracOfCell[0];
    
            cell[0] = centralIntersection[0] - 1;
            mf[p.m_grid](cell)                += p.m_data[0] * (1-fracOfCell[0]);
#endif
            
#if (BL_SPACEDIM == 2)
            mf[p.m_grid](centralIntersection) += p.m_data[0] *    fracOfCell[0]  *    fracOfCell[1] ;
    
            cell[0] = centralIntersection[0] - 1;
            mf[p.m_grid](cell)                += p.m_data[0] * (1-fracOfCell[0]) *    fracOfCell[1] ;
    
            cell = centralIntersection;
            cell[1] = centralIntersection[1] - 1;
            mf[p.m_grid](cell)                += p.m_data[0] *    fracOfCell[0]  * (1-fracOfCell[1]);
    
            cell[0] = centralIntersection[0] - 1;
            mf[p.m_grid](cell)                += p.m_data[0] * (1-fracOfCell[0]) * (1-fracOfCell[1]);
#endif


#if (BL_SPACEDIM == 3)
            mf[p.m_grid](centralIntersection) += p.m_data[0] *    fracOfCell[0]  *    fracOfCell[1]  *    fracOfCell[2] ;
    
            cell[0] = centralIntersection[0] - 1;
            mf[p.m_grid](cell)                += p.m_data[0] * (1-fracOfCell[0]) *    fracOfCell[1]  *    fracOfCell[2] ;
    
            cell = centralIntersection;
            cell[1] = centralIntersection[1] - 1;
            mf[p.m_grid](cell)                += p.m_data[0] *    fracOfCell[0]  * (1-fracOfCell[1]) *    fracOfCell[2] ;
    
            cell[0] = centralIntersection[0] - 1;
            mf[p.m_grid](cell)                += p.m_data[0] * (1-fracOfCell[0]) * (1-fracOfCell[1]) *    fracOfCell[2] ;

            
            centralIntersection[3] = centralIntersection[3] - 1; //ugly hack - but this way I can "reuse" the code fr2m 2d...
            mf[p.m_grid](centralIntersection) += p.m_data[0] *    fracOfCell[0]  *    fracOfCell[1]  * (1-fracOfCell[2]);
    
            cell[0] = centralIntersection[0] - 1;
            mf[p.m_grid](cell)                += p.m_data[0] * (1-fracOfCell[0]) *    fracOfCell[1]  * (1-fracOfCell[2]);
    
            cell = centralIntersection;
            cell[1] = centralIntersection[1] - 1;
            mf[p.m_grid](cell)                += p.m_data[0] *    fracOfCell[0]  * (1-fracOfCell[1]) * (1-fracOfCell[2]);
    
            cell[0] = centralIntersection[0] - 1;
            mf[p.m_grid](cell)                += p.m_data[0] * (1-fracOfCell[0]) * (1-fracOfCell[1]) * (1-fracOfCell[2]);
#endif
        }
    }
    mf.SumBoundary();
    geom.SumPeriodicBoundary(mf);
    mf.mult(1./cell_vol);
}

template <int N>
void
ParticleContainer<N>::movePredict(MultiFab& grav_vector, int level, Real timestep)
{
    const int MyProc = ParallelDescriptor::MyProc();

    const DistributionMapping& dm = grav_vector.DistributionMap();


    for (iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        ParticleType& p = *it;

        if (p.m_id > 0 && p.m_lev == level && dm[p.m_grid] == MyProc)
        {
            //note: m_data[0] is mass, 1 is v_x, ...
            D_TERM(p.m_data[1] += timestep * grav_vector[p.m_grid](p.m_cell,0);,
                   p.m_data[2] += timestep * grav_vector[p.m_grid](p.m_cell,1);,
                   p.m_data[3] += timestep * grav_vector[p.m_grid](p.m_cell,2););
            D_TERM(p.m_pos[0]  += timestep * p.m_data[1];,
                   p.m_pos[1]  += timestep * p.m_data[2];,
                   p.m_pos[2]  += timestep * p.m_data[3];);
        }
    }
//NO redistribution here, this is done in the corrector-step (which is hopefully sufficient).
}


template <int N>
void
ParticleContainer<N>::moveCorrect(MultiFab& grav_vector_old, MultiFab& grav_vector, int level, Real timestep)
{
    const int MyProc = ParallelDescriptor::MyProc();

    const DistributionMapping& dm = grav_vector.DistributionMap();


    for (iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        ParticleType& p = *it;

        if (p.m_id > 0 && p.m_lev == level && dm[p.m_grid] == MyProc)
        {
            D_TERM(p.m_pos[0]  -= 0.5 * timestep * p.m_data[1];,
                   p.m_pos[1]  -= 0.5 * timestep * p.m_data[2];,
                   p.m_pos[2]  -= 0.5 * timestep * p.m_data[3];);
            //note: m_data[0] is mass, 1 is v_x, ...
            D_TERM(p.m_data[1] += 0.5 * timestep * ( grav_vector[p.m_grid](p.m_cell,0) - grav_vector_old[p.m_grid](p.m_cell,0) );, 
                   p.m_data[2] += 0.5 * timestep * ( grav_vector[p.m_grid](p.m_cell,1) - grav_vector_old[p.m_grid](p.m_cell,1) );,
                   p.m_data[3] += 0.5 * timestep * ( grav_vector[p.m_grid](p.m_cell,2) - grav_vector_old[p.m_grid](p.m_cell,2) ););
            D_TERM(p.m_pos[0]  += 0.5 * timestep * p.m_data[1];,
                   p.m_pos[1]  += 0.5 * timestep * p.m_data[2];,
                   p.m_pos[2]  += 0.5 * timestep * p.m_data[3];);
        }
        //
        // Here's where we need to deal with boundary conditions.
        //
        // For the moment we assume all boundaries are periodic.
        //
        // Shift the particle back into the domain.
        //
        bool shifted = ParticleBase::PeriodicShift(p,m_amr);
        //
        // Modulo any arithmetic funny business the particle should in the problem domain.
        //
        BL_ASSERT(shifted == false || ParticleBase::Where(p,m_amr));
    }
    Redistribute();
}

template <int N>
void
ParticleContainer<N>::moveKickDrift(MultiFab& grav_vector, int level, Real timestep)
{
    const int MyProc = ParallelDescriptor::MyProc();

    const DistributionMapping& dm = grav_vector.DistributionMap();


    for (iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        ParticleType& p = *it;

        if (p.m_id > 0 && p.m_lev == level && dm[p.m_grid] == MyProc)
        {
            //note: m_data[0] is mass, 1 is v_x, ...
            D_TERM(p.m_data[1] += 0.5 * timestep * grav_vector[p.m_grid](p.m_cell,0);, //how do I use grav_vec correctly? this seems to be some kind of s...
                   p.m_data[2] += 0.5 * timestep * grav_vector[p.m_grid](p.m_cell,1);,
                   p.m_data[3] += 0.5 * timestep * grav_vector[p.m_grid](p.m_cell,2););
            D_TERM(p.m_pos[0]  +=       timestep * p.m_data[1];,
                   p.m_pos[1]  +=       timestep * p.m_data[2];,
                   p.m_pos[2]  +=       timestep * p.m_data[3];);
        }
        //
        // Here's where we need to deal with boundary conditions.
        //
        // For the moment we assume all boundaries are periodic.
        //
        // Shift the particle back into the domain.
        //
        bool shifted = ParticleBase::PeriodicShift(p,m_amr);
        //
        // Modulo any arithmetic funny business the particle should in the problem domain.
        //
        BL_ASSERT(shifted == false || ParticleBase::Where(p,m_amr));
    }
    Redistribute();
}

template <int N>
void
ParticleContainer<N>::moveKick(MultiFab& grav_vector, int level, Real timestep)
{
    const int MyProc = ParallelDescriptor::MyProc();

    const DistributionMapping& dm = grav_vector.DistributionMap();


    for (iterator it = m_particles.begin(); it != m_particles.end(); ++it)
    {
        ParticleType& p = *it;

        if (p.m_id > 0 && p.m_lev == level && dm[p.m_grid] == MyProc)
        {
            //note: m_data[0] is mass, 1 is v_x, ...
            D_TERM(p.m_data[1] += 0.5 * timestep * grav_vector[p.m_grid](p.m_cell,0);, //how do I use grav_vec correctly? this seems to be some kind of s...
                   p.m_data[2] += 0.5 * timestep * grav_vector[p.m_grid](p.m_cell,1);,
                   p.m_data[3] += 0.5 * timestep * grav_vector[p.m_grid](p.m_cell,2););
        }
    }
    // No need for redistribution, we only change the velocity.
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

        if (p.m_id <= 0) continue;

        if (!ParticleBase::Where(p,m_amr))
            return false;
    }
    return true;
}

#endif /*_PARTICLES_H_*/
