#ifndef AMREX_INTEGRATOR_BASE_H
#define AMREX_INTEGRATOR_BASE_H
#include <AMReX_Config.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

#if defined(AMREX_PARTICLES)
#include <AMReX_Particles.H>
#endif

#include <functional>
#include <type_traits>

namespace amrex {

template<class T, typename Tv = void> struct IntegratorOps;

#if defined(AMREX_PARTICLES)
template<class T>
struct IntegratorOps<T, typename std::enable_if<std::is_base_of<amrex::ParticleContainerBase, T>::value>::type>
{

    static void CreateLike (amrex::Vector<std::unique_ptr<T> >& V, const T& Other)
    {
        // Emplace a new T in V with the same size as Other and get a reference
        V.emplace_back(std::make_unique<T>(Other.Geom(0), Other.ParticleDistributionMap(0), Other.ParticleBoxArray(0)));
        T& pc = *V[V.size()-1];

        // We want the particles to have all the same position, cpu, etc.
        // as in Other, so do a copy from Other to our new particle container.
        Copy(pc, Other);
    }

    static void Copy (T& Y, const T& Other)
    {
        // Copy the contents of Other into Y
        const bool local = true;
        Y.copyParticles(Other, local);
    }

    static void Saxpy (T& Y, const amrex::Real a, T& X)
    {
        // Calculate Y += a * X using a particle-level saxpy function supplied by the particle container T
        typedef amrex::ParIter<T::NStructReal, T::NStructInt, T::NArrayReal, T::NArrayInt> TParIter;
        typedef amrex::Particle<T::NStructReal, T::NStructInt> ParticleType;

        int lev = 0;
        TParIter pty(Y, lev);
        TParIter ptx(X, lev);

        auto checkValid = [&]() -> bool {
            bool pty_v = pty.isValid();
            bool ptx_v = ptx.isValid();
            AMREX_ASSERT(pty_v == ptx_v);
            return pty_v && ptx_v;
        };

        auto ptIncrement = [&](){ ++pty; ++ptx; };

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (; checkValid(); ptIncrement())
        {
            const int npy  = pty.numParticles();
            const int npx  = ptx.numParticles();
            AMREX_ALWAYS_ASSERT(npy == npx);

            ParticleType* psy = &(pty.GetArrayOfStructs()[0]);
            ParticleType* psx = &(ptx.GetArrayOfStructs()[0]);

            auto particle_apply_rhs = T::particle_apply_rhs;

            amrex::ParallelFor ( npy, [=] AMREX_GPU_DEVICE (int i) {
                ParticleType& py = psy[i];
                const ParticleType& px = psx[i];
                particle_apply_rhs(py, a, px);
            });
        }
    }

};
#endif

template<class T>
struct IntegratorOps<T, typename std::enable_if<std::is_same<amrex::Vector<amrex::MultiFab>, T>::value>::type>
{

    static void CreateLike (amrex::Vector<std::unique_ptr<T> >& V, const T& Other, bool Grow = false)
    {
        // Emplace a new T in V with the same size as Other
        V.emplace_back(std::make_unique<T>());
        for (auto const& other_mf : Other) {
            IntVect nGrow = Grow ? other_mf.nGrowVect() : IntVect(0);
            V.back()->push_back(amrex::MultiFab(other_mf.boxArray(), other_mf.DistributionMap(), other_mf.nComp(), nGrow));
        }
    }

    static void Copy (T& Y, const T& Other, const Vector<int>& scomp={}, const Vector<int>& ncomp={}, bool Grow = true)
    {
        // Copy the contents of Other into Y
        const int size = Y.size();
        bool specify_components = !scomp.empty() && ncomp.size() == scomp.size();
        for (int i = 0; i < size; ++i) {
            IntVect nGrow = Grow ? Other[i].nGrowVect() : IntVect(0);
            const int iscomp = specify_components ? scomp[i] : 0;
            const int incomp = specify_components ? ncomp[i] : Other[i].nComp();
            if (incomp > 0) {
                amrex::MultiFab::Copy(Y[i], Other[i], iscomp, iscomp, incomp, nGrow);
            }
        }
    }

    static void Saxpy (T& Y, const amrex::Real a, const T& X, const Vector<int>& scomp={}, const Vector<int>& ncomp={}, bool Grow = false)
    {
        // Calculate Y += a * X
        const int size = Y.size();
        bool specify_components = !scomp.empty() && ncomp.size() == scomp.size();
        for (int i = 0; i < size; ++i) {
            IntVect nGrow = Grow ? X[i].nGrowVect() : IntVect(0);
            const int iscomp = specify_components ? scomp[i] : 0;
            const int incomp = specify_components ? ncomp[i] : X[i].nComp();
            if (incomp > 0) {
                amrex::MultiFab::Saxpy(Y[i], a, X[i], iscomp, iscomp, incomp, nGrow);
            }
        }
    }

};

template<class T>
struct IntegratorOps<T, typename std::enable_if<std::is_same<amrex::MultiFab, T>::value>::type>
{

    static void CreateLike (amrex::Vector<std::unique_ptr<T> >& V, const T& Other, bool Grow = false)
    {
        // Emplace a new T in V with the same size as Other
        IntVect nGrow = Grow ? Other.nGrowVect() : IntVect(0);
        V.emplace_back(std::make_unique<T>(Other.boxArray(), Other.DistributionMap(), Other.nComp(), nGrow));
    }

    static void Copy (T& Y, const T& Other, const int scomp=0, const int ncomp=-1, bool Grow = true)
    {
        // Copy the contents of Other into Y
        IntVect nGrow = Grow ? Other.nGrowVect() : IntVect(0);
        const int mf_ncomp = ncomp > 0 ? ncomp : Other.nComp();
        amrex::MultiFab::Copy(Y, Other, scomp, scomp, mf_ncomp, nGrow);
    }

    static void Saxpy (T& Y, const amrex::Real a, const T& X, const int scomp=0, const int ncomp=-1, bool Grow = false)
    {
        // Calculate Y += a * X
        IntVect nGrow = Grow ? X.nGrowVect() : IntVect(0);
        const int mf_ncomp = ncomp > 0 ? ncomp : X.nComp();
        amrex::MultiFab::Saxpy(Y, a, X, scomp, scomp, mf_ncomp, nGrow);
    }

};

template<class T>
class IntegratorBase
{
private:
   /**
    * \brief Fun is the right-hand-side function the integrator will use.
    */
    std::function<void(T&, const T&, const amrex::Real)> Fun;

   /**
    * \brief FastFun is the fast timescale right-hand-side function for a multirate integration problem.
    */
    std::function<void(T&, T&, const T&, const amrex::Real)> FastFun;

protected:
   /**
    * \brief Integrator timestep size (Real)
    */
    amrex::Real timestep;

   /**
    * \brief For multirate problems, the ratio of slow timestep size / fast timestep size (int)
    */
    int slow_fast_timestep_ratio = 0;

   /**
    * \brief For multirate problems, the fast timestep size (Real)
    */
    Real fast_timestep = 0.0;

   /**
    * \brief The post_update function is called by the integrator on state data before using it to evaluate a right-hand side.
    */
    std::function<void (T&, amrex::Real)> post_update;

public:
    IntegratorBase () = default;

    IntegratorBase (const T& /* S_data */) {}

    virtual ~IntegratorBase () = default;

    virtual void initialize (const T& S_data) = 0;

    void set_rhs (std::function<void(T&, const T&, const amrex::Real)> F)
    {
        Fun = F;
    }

    void set_fast_rhs (std::function<void(T&, T&, const T&, const amrex::Real)> F)
    {
        FastFun = F;
    }

    void set_slow_fast_timestep_ratio (const int timestep_ratio = 1)
    {
        slow_fast_timestep_ratio = timestep_ratio;
    }

    void set_fast_timestep (const Real fast_dt = 1.0)
    {
        fast_timestep = fast_dt;
    }

    void set_post_update (std::function<void (T&, amrex::Real)> F)
    {
        post_update = F;
    }

    std::function<void (T&, amrex::Real)> get_post_update ()
    {
        return post_update;
    }

    std::function<void(T&, const T&, const amrex::Real)> get_rhs ()
    {
        return Fun;
    }

    std::function<void(T&, T&, const T&, const amrex::Real)> get_fast_rhs ()
    {
        return FastFun;
    }

    int get_slow_fast_timestep_ratio ()
    {
        return slow_fast_timestep_ratio;
    }

    Real get_fast_timestep ()
    {
        return fast_timestep;
    }

    void rhs (T& S_rhs, const T& S_data, const amrex::Real time)
    {
        Fun(S_rhs, S_data, time);
    }

    void fast_rhs (T& S_rhs, T& S_extra, const T& S_data, const amrex::Real time)
    {
        FastFun(S_rhs, S_extra, S_data, time);
    }

    virtual amrex::Real advance (T& S_old, T& S_new, amrex::Real time, amrex::Real dt) = 0;

    virtual void time_interpolate (const T& S_new, const T& S_old, amrex::Real timestep_fraction, T& data) = 0;

    virtual void map_data (std::function<void(T&)> Map) = 0;
};

}

#endif
