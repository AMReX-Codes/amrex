#ifndef AMREX_FE_INTEGRATOR_H
#define AMREX_FE_INTEGRATOR_H
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_IntegratorBase.H>
#include <functional>

namespace amrex {

template<class T>
class FEIntegrator : public IntegratorBase<T>
{
private:
    typedef IntegratorBase<T> BaseT;

    amrex::Vector<std::unique_ptr<T> > F_nodes;

    void initialize_stages (const T& S_data)
    {
        IntegratorOps<T>::CreateLike(F_nodes, S_data);
    }

public:
    FEIntegrator () {}

    FEIntegrator (const T& S_data)
    {
        initialize(S_data);
    }

    virtual ~FEIntegrator () {}

    void initialize (const T& S_data) override
    {
        initialize_stages(S_data);
    }

    amrex::Real advance (T& S_old, T& S_new, amrex::Real time, const amrex::Real time_step) override
    {
        BaseT::timestep = time_step;
        // Assume before advance() that S_old is valid data at the current time ("time" argument)
        // So we initialize S_new by copying the old state.
        IntegratorOps<T>::Copy(S_new, S_old);

        // F = RHS(S, t)
        T& F = *F_nodes[0];
        BaseT::rhs(F, S_new, time);

        // S_new += timestep * dS/dt
        IntegratorOps<T>::Saxpy(S_new, BaseT::timestep, F);

        // Call the post-update hook for S_new
        BaseT::post_update(S_new, time + BaseT::timestep);

        // Return timestep
        return BaseT::timestep;
    }

    virtual void time_interpolate (const T& /* S_new */, const T& /* S_old */, amrex::Real /* timestep_fraction */, T& /* data */) override
    {
        amrex::Error("Time interpolation not yet supported by forward euler integrator.");
    }

    virtual void map_data (std::function<void(T&)> Map) override
    {
        for (auto& F : F_nodes) {
            Map(*F);
        }
    }

};

}

#endif
