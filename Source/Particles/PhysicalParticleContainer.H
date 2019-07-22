#ifndef WARPX_PhysicalParticleContainer_H_
#define WARPX_PhysicalParticleContainer_H_

#include <map>

#include <AMReX_IArrayBox.H>

#include <PlasmaInjector.H>
#include <WarpXParticleContainer.H>

class PhysicalParticleContainer
    : public WarpXParticleContainer
{
public:
    PhysicalParticleContainer (amrex::AmrCore* amr_core,
                               int ispecies,
                               const std::string& name);

    PhysicalParticleContainer (amrex::AmrCore* amr_core);

    virtual ~PhysicalParticleContainer () {}

    virtual void InitData () override;

#ifdef WARPX_DO_ELECTROSTATIC    
    virtual void FieldGatherES(const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                               const amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > >& masks) override;

    virtual void EvolveES (const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                                 amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                           amrex::Real t, amrex::Real dt) override;
#endif // WARPX_DO_ELECTROSTATIC
    
    virtual void FieldGather(int lev,
                             const amrex::MultiFab& Ex,
                             const amrex::MultiFab& Ey,
                             const amrex::MultiFab& Ez,
                             const amrex::MultiFab& Bx,
                             const amrex::MultiFab& By,
                             const amrex::MultiFab& Bz) final;

    virtual void Evolve (int lev,
			 const amrex::MultiFab& Ex,
                         const amrex::MultiFab& Ey,
                         const amrex::MultiFab& Ez,
			 const amrex::MultiFab& Bx,
                         const amrex::MultiFab& By,
                         const amrex::MultiFab& Bz,
			 amrex::MultiFab& jx,
                         amrex::MultiFab& jy,
                         amrex::MultiFab& jz,
                         amrex::MultiFab* cjx,
                         amrex::MultiFab* cjy,
                         amrex::MultiFab* cjz,
                         amrex::MultiFab* rho,
                         amrex::MultiFab* crho,
                         const amrex::MultiFab* cEx,
                         const amrex::MultiFab* cEy,
                         const amrex::MultiFab* cEz,
                         const amrex::MultiFab* cBx,
                         const amrex::MultiFab* cBy,
                         const amrex::MultiFab* cBz,
                         amrex::Real t,
                         amrex::Real dt) override;

    virtual void PushPX(WarpXParIter& pti,
	                amrex::Cuda::ManagedDeviceVector<amrex::Real>& xp,
                        amrex::Cuda::ManagedDeviceVector<amrex::Real>& yp,
                        amrex::Cuda::ManagedDeviceVector<amrex::Real>& zp,
                        amrex::Cuda::ManagedDeviceVector<amrex::Real>& giv,
                        amrex::Real dt);

    virtual void PushP (int lev, amrex::Real dt,
                        const amrex::MultiFab& Ex,
                        const amrex::MultiFab& Ey,
                        const amrex::MultiFab& Ez,
                        const amrex::MultiFab& Bx,
                        const amrex::MultiFab& By,
                        const amrex::MultiFab& Bz) override;
                        
    void copy_attribs(WarpXParIter& pti,const amrex::Real* xp,
                        const amrex::Real* yp, const amrex::Real* zp);

    virtual void PostRestart () final {}

    void SplitParticles(int lev);

    // Inject particles in Box 'part_box'
    virtual void AddParticles (int lev);
    void AddPlasma(int lev, amrex::RealBox part_realbox = amrex::RealBox());
    void AddPlasmaCPU (int lev, amrex::RealBox part_realbox);
#ifdef AMREX_USE_GPU
    void AddPlasmaGPU (int lev, amrex::RealBox part_realbox);
#endif

    void MapParticletoBoostedFrame(amrex::Real& x, amrex::Real& y, amrex::Real& z, std::array<amrex::Real, 3>& u);

    void AddGaussianBeam(amrex::Real x_m, amrex::Real y_m, amrex::Real z_m,
                         amrex::Real x_rms, amrex::Real y_rms, amrex::Real z_rms,
                         amrex::Real q_tot, long npart, int do_symmetrize);

    void CheckAndAddParticle(amrex::Real x, amrex::Real y, amrex::Real z,
                             std::array<amrex::Real, 3> u, amrex::Real weight);

    virtual void GetParticleSlice(const int direction, const amrex::Real z_old,
                                  const amrex::Real z_new, const amrex::Real t_boost, 
                                  const amrex::Real t_lab, const amrex::Real dt,
                                  DiagnosticParticles& diagnostic_particles) final;

    virtual void ConvertUnits (ConvertDirection convert_dir) override;

protected:

    std::string species_name;
    std::unique_ptr<PlasmaInjector> plasma_injector;

    // When true, adjust the transverse particle positions accounting
    // for the difference between the Lorentz transformed time of the
    // particle and the time of the boosted frame.
    bool boost_adjust_transverse_positions = false;
    bool do_backward_propagation = false;

    long NumParticlesToAdd (const amrex::Box& overlap_box,
			    const amrex::RealBox& overlap_realbox,
			    const amrex::RealBox& tile_real_box,
			    const amrex::RealBox& particle_real_box);
  
    int GetRefineFac(const amrex::Real x, const amrex::Real y, const amrex::Real z);
    std::unique_ptr<amrex::IArrayBox> m_refined_injection_mask = nullptr;

    // Inject particles during the whole simulation
    void ContinuousInjection(const amrex::RealBox& injection_box) override;

};

#endif
