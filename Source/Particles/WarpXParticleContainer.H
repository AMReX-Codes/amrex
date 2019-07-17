#ifndef WARPX_WarpXParticleContainer_H_
#define WARPX_WarpXParticleContainer_H_

#include <memory>

#include <AMReX_Particles.H>
#include <AMReX_AmrCore.H>

enum struct ConvertDirection{WarpX_to_SI, SI_to_WarpX};

struct PIdx
{
    enum { // Particle Attributes stored in amrex::ParticleContainer's struct of array
	w = 0,  // weight
	ux, uy, uz, Ex, Ey, Ez, Bx, By, Bz,
#ifdef WARPX_RZ
        theta, // RZ needs all three position components
#endif
	nattribs
    };
};

struct DiagIdx
{
    enum {
        w = 0,
        x, y, z, ux, uy, uz,
        nattribs
    };
};

namespace ParticleStringNames
{
    const std::map<std::string, int> to_index = {
        {"w",     PIdx::w    },
        {"ux",    PIdx::ux   },
        {"uy",    PIdx::uy   },
        {"uz",    PIdx::uz   },
        {"Ex",    PIdx::Ex   },
        {"Ey",    PIdx::Ey   },
        {"Ez",    PIdx::Ez   },
        {"Bx",    PIdx::Bx   },
        {"By",    PIdx::By   },
        {"Bz",    PIdx::Bz   }
    };
}

class WarpXParIter
    : public amrex::ParIter<0,0,PIdx::nattribs>
{
public:
    using amrex::ParIter<0,0,PIdx::nattribs>::ParIter;

    WarpXParIter (ContainerType& pc, int level);

#if (AMREX_SPACEDIM == 2)
    void GetPosition (amrex::Cuda::ManagedDeviceVector<amrex::Real>& x,
                      amrex::Cuda::ManagedDeviceVector<amrex::Real>& y,
                      amrex::Cuda::ManagedDeviceVector<amrex::Real>& z) const;
    void SetPosition (const amrex::Cuda::ManagedDeviceVector<amrex::Real>& x,
                      const amrex::Cuda::ManagedDeviceVector<amrex::Real>& y,
                      const amrex::Cuda::ManagedDeviceVector<amrex::Real>& z);
#endif
    const std::array<RealVector, PIdx::nattribs>& GetAttribs () const { 
        return GetStructOfArrays().GetRealData(); 
    }
    
    std::array<RealVector, PIdx::nattribs>& GetAttribs () { 
        return GetStructOfArrays().GetRealData(); 
    }

    const RealVector& GetAttribs (int comp) const { 
        return GetStructOfArrays().GetRealData(comp);
    }
    
    RealVector& GetAttribs (int comp) { 
        return GetStructOfArrays().GetRealData(comp);
    }
};

class MultiParticleContainer;

class WarpXParticleContainer
    : public amrex::ParticleContainer<0,0,PIdx::nattribs>
{
public:
    friend MultiParticleContainer;

    // amrex::StructOfArrays with DiagIdx::nattribs amrex::Real components 
    // and 0 int components for the particle data.
    using DiagnosticParticleData = amrex::StructOfArrays<DiagIdx::nattribs, 0>;
    // DiagnosticParticles is a vector, with one element per MR level.
    // DiagnosticParticles[lev] is typically a key-value pair where the key is 
    // a pair [grid_index, tile_index], and the value is the corresponding 
    // DiagnosticParticleData (see above) on this tile.
    using DiagnosticParticles = amrex::Vector<std::map<std::pair<int, int>, DiagnosticParticleData> >;
    
    WarpXParticleContainer (amrex::AmrCore* amr_core, int ispecies);
    virtual ~WarpXParticleContainer() {}

    virtual void InitData () = 0;

    virtual void FieldGatherES (const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                                const amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > >& masks) {}

    virtual void FieldGatherFortran (int lev,
                                     const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
                                     const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz) {}

#ifdef WARPX_DO_ELECTROSTATIC    
    virtual void EvolveES (const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                                 amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho, 
                           amrex::Real t, amrex::Real dt) = 0;
#endif // WARPX_DO_ELECTROSTATIC
    
    virtual void Evolve (int lev,
			 const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
			 const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
			 amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz,
                         amrex::MultiFab* cjx, amrex::MultiFab* cjy, amrex::MultiFab* cjz,
                         amrex::MultiFab* rho, amrex::MultiFab* crho,
                         const amrex::MultiFab* cEx, const amrex::MultiFab* cEy, const amrex::MultiFab* cEz,
                         const amrex::MultiFab* cBx, const amrex::MultiFab* cBy, const amrex::MultiFab* cBz,
                         amrex::Real t, amrex::Real dt) = 0;

    virtual void PostRestart () = 0;

    virtual void GetParticleSlice(const int direction,     const amrex::Real z_old,
                                  const amrex::Real z_new, const amrex::Real t_boost, 
                                  const amrex::Real t_lab, const amrex::Real dt,
                                  DiagnosticParticles& diagnostic_particles) {}
    
    void AllocData ();

    ///
    /// This pushes the particle positions by one half time step.
    /// It is used to desynchronize the particles after initializaton
    /// or when restarting from a checkpoint.
    /// This is the electrostatic version of the particle push.
    /// 
    void PushXES (amrex::Real dt);

    ///
    /// This pushes the particle positions by one half time step.
    /// It is used to desynchronize the particles after initializaton
    /// or when restarting from a checkpoint.
    /// This is the electromagnetic version of the particle push.
    /// 
    void PushX (         amrex::Real dt);
    void PushX (int lev, amrex::Real dt);

    ///
    /// This pushes the particle momenta by dt.
    /// 
    virtual void PushP (int lev, amrex::Real dt,
                        const amrex::MultiFab& Ex,
                        const amrex::MultiFab& Ey,
                        const amrex::MultiFab& Ez,
                        const amrex::MultiFab& Bx,
                        const amrex::MultiFab& By,
                        const amrex::MultiFab& Bz) = 0;

    void DepositCharge(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                       bool local = false);
    std::unique_ptr<amrex::MultiFab> GetChargeDensity(int lev, bool local = false);

    virtual void DepositCharge(WarpXParIter& pti,
                               RealVector& wp,
                               amrex::MultiFab* rhomf,
                               amrex::MultiFab* crhomf,
                               int icomp,
                               const long np_current,
                               const long np,
                               int thread_num,
                               int lev );

    virtual void DepositCurrent(WarpXParIter& pti,
                                RealVector& wp,
                                RealVector& uxp,
                                RealVector& uyp,
                                RealVector& uzp,
                                amrex::MultiFab* jx,
                                amrex::MultiFab* jy,
                                amrex::MultiFab* jz,
                                const long offset,
                                const long np_to_depose,
                                int thread_num,
                                int lev,
                                int depos_lev,
                                amrex::Real dt);

    virtual void DepositCurrentFortran(WarpXParIter& pti,
                                       RealVector& wp,
                                       RealVector& uxp,
                                       RealVector& uyp,
                                       RealVector& uzp,
                                       amrex::MultiFab* jx,
                                       amrex::MultiFab* jy,
                                       amrex::MultiFab* jz,
                                       const long offset,
                                       const long np_to_depose,
                                       int thread_num,
                                       int lev,
                                       int depos_lev,
                                       amrex::Real dt);

    // If particles start outside of the domain, ContinuousInjection 
    // makes sure that they are initialized when they enter the domain, and 
    // NOT before. Virtual function, overriden by derived classes. 
    // Current status:
    // PhysicalParticleContainer: implemented.
    // LaserParticleContainer: implemented.
    // RigidInjectedParticleContainer: not implemented.
    virtual void ContinuousInjection(const amrex::RealBox& injection_box) {}
    // Update optional sub-class-specific injection location.
    virtual void UpdateContinuousInjectionPosition(amrex::Real dt) {}

    /// 
    /// This returns the total charge for all the particles in this ParticleContainer.
    /// This is needed when solving Poisson's equation with periodic boundary conditions.
    ///
    amrex::Real sumParticleCharge(bool local = false);

    std::array<amrex::Real, 3> meanParticleVelocity(bool local = false);

    amrex::Real maxParticleVelocity(bool local = false);

    void AddNParticles (int lev,
                        int n, const amrex::Real* x, const amrex::Real* y, const amrex::Real* z,
			const amrex::Real* vx, const amrex::Real* vy, const amrex::Real* vz,
			int nattr, const amrex::Real* attr, int uniqueparticles, int id=-1);

    void AddOneParticle (int lev, int grid, int tile,
                         amrex::Real x, amrex::Real y, amrex::Real z,
                         std::array<amrex::Real,PIdx::nattribs>& attribs);

    void AddOneParticle (ParticleTileType& particle_tile,
                         amrex::Real x, amrex::Real y, amrex::Real z,
                         std::array<amrex::Real,PIdx::nattribs>& attribs);

    virtual void ReadHeader (std::istream& is);

    virtual void WriteHeader (std::ostream& os) const;

    virtual void ConvertUnits (ConvertDirection convert_dir){};

    static void ReadParameters ();

    static int NextID () { return ParticleType::NextID(); }

    bool do_splitting = false;

    // split along axes (0) or diagonals (1)
    int split_type = 0;

    using amrex::ParticleContainer<0, 0, PIdx::nattribs>::AddRealComp; 
    using amrex::ParticleContainer<0, 0, PIdx::nattribs>::AddIntComp;
   
    void AddRealComp (const std::string& name, bool comm=true)
    {
        particle_comps[name] = NumRealComps();
        AddRealComp(comm);
    }

    void AddIntComp (const std::string& name, bool comm=true)
    {
        particle_comps[name] = NumIntComps();
        AddIntComp(comm);
    }

    int DoBoostedFrameDiags () const { return do_boosted_frame_diags; }    

protected:

    std::map<std::string, int> particle_comps;
    
    int species_id;

    amrex::Real charge;
    amrex::Real mass;

    bool deposit_on_main_grid = false;

    static int do_not_push;

    // Whether to allow particles outside of the simulation domain to be 
    // initialized when they enter the domain.
    // This is currently required because continuous injection does not 
    // support all features allowed by direct injection.
    int do_continuous_injection = 0;

    int do_boosted_frame_diags = 1;

    amrex::Vector<amrex::FArrayBox> local_rho;
    amrex::Vector<amrex::FArrayBox> local_jx;
    amrex::Vector<amrex::FArrayBox> local_jy;
    amrex::Vector<amrex::FArrayBox> local_jz;

    amrex::Vector<amrex::Cuda::ManagedDeviceVector<amrex::Real> > m_xp, m_yp, m_zp, m_giv;

    // Whether to dump particle quantities. 
    // If true, particle position is always dumped.
    int plot_species = 1;
    // For all particle attribs (execept position), whether or not 
    // to dump to file.
    amrex::Vector<int> plot_flags;
    // list of names of attributes to dump.
    amrex::Vector<std::string> plot_vars;

private:
    virtual void particlePostLocate(ParticleType& p, const amrex::ParticleLocData& pld,
                                    const int lev) override;

};

#endif
