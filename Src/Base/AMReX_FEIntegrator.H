#ifndef AMREX_FE_INTEGRATOR_H
#define AMREX_FE_INTEGRATOR_H
#include <functional>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include "AMReX_IntegratorBase.H"

template<class T>
class FEIntegrator : public IntegratorBase<T>
{
private:
    typedef IntegratorBase<T> BaseT;

    amrex::Vector<std::unique_ptr<T> > F_nodes;

public:
    FEIntegrator(T& S_data)
    {
        IntegratorOps<T>::CreateLike(F_nodes, S_data);
    }

    virtual ~FEIntegrator() {}

    amrex::Real advance(T& S_old, T& S_new, amrex::Real time, const amrex::Real timestep)
    {
        // Assume before advance() that S_old is valid data at the current time ("time" argument)
        // So we initialize S_new by copying the old state.
        IntegratorOps<T>::Copy(S_new, S_old);

        // F = RHS(S, t)
        T& F = *F_nodes[0];
        BaseT::rhs(F, S_new, time);

        // S_new += timestep * dS/dt
        IntegratorOps<T>::Saxpy(S_new, timestep, F);

        // Call the post-update hook for S_new
        BaseT::post_update(S_new, time + timestep);

        // Return timestep
        return timestep;
    }

    void time_interpolate(const T& S_new, const T& S_old, amrex::Real timestep_fraction, T& data)
    {
        amrex::Error("Time interpolation not yet supported by forward euler integrator.");
    }

    void map_data(std::function<void(T&)> Map)
    {
        for (auto& F : F_nodes) {
            Map(*F);
        }
    }

};

#endif
