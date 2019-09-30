#ifndef IONIZATIONEVENT_H_
#define IONIZATIONEVENT_H_

#include "WarpXParticleContainer.H"

class copyAndTransformParticle
{
public:
    void operator () (
        int is, int ip, int pid_product, int cpuid, int do_boosted_product,
        const amrex::GpuArray<amrex::ParticleReal*,3>& runtime_uold_source,
        WarpXParticleContainer::ParticleType& p_source,
        WarpXParticleContainer::ParticleType& p_product,
        const amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs>& attribs_source,
        const amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs>& attribs_product,
        const amrex::GpuArray<amrex::ParticleReal*,6>& runtime_attribs_product) const
    {
        amrex::Print()<<"is = "<<is<<'\n';

        // Copy particle from source to product: AoS
        p_product.id() = pid_product + ip;
        p_product.cpu() = cpuid;
        p_product.pos(0) = p_source.pos(0);
        p_product.pos(1) = p_source.pos(1);
#if (AMREX_SPACEDIM == 3)
        p_product.pos(2) = p_source.pos(2);
#endif
        // Copy particle from source to product: SoA
        for (int ia = 0; ia < PIdx::nattribs; ++ia) {
            attribs_product[ia][ip] = attribs_source[ia][is];
        }
        // Update xold etc. if boosted frame diagnostics required
        // for product species. Fill runtime attribs with a copy of
        // current properties (xold = x etc.).
        if (do_boosted_product) {
            runtime_attribs_product[0][ip] = p_source.pos(0);
            runtime_attribs_product[1][ip] = p_source.pos(1);
            runtime_attribs_product[2][ip] = p_source.pos(2);
            runtime_attribs_product[3][ip] = runtime_uold_source[0][ip];
            runtime_attribs_product[4][ip] = runtime_uold_source[1][ip];
            runtime_attribs_product[5][ip] = runtime_uold_source[2][ip];
        }
    }
};

#endif // IONIZATIONEVENT_H_
