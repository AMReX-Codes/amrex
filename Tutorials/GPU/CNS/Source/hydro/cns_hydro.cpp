
#include "CNS_K.H"
#include "CNS_index_macros.H"

using namespace amrex;

AMREX_GPU_DEVICE
void cns_ctoprim (Box const& bx, FArrayBox const& u, FArrayBox & q)
{
    const auto up0 = u.stridedPtr(bx);
    const auto qp0 = q.stridedPtr(bx);
    const auto len = bx.length3d();
    const Real smallr = 1.e-19;
    const Real smallp = 1.e-10;
    const Real gamma = 1.4;

    for         (int k = 0; k < len[2]; ++k) {
        for     (int j = 0; j < len[1]; ++j) {
            for (int i = 0; i < len[0]; ++i) {
                Real rho  = *up0(i,j,k,URHO);
                Real mx   = *up0(i,j,k,UMX);
                Real my   = *up0(i,j,k,UMY);
                Real mz   = *up0(i,j,k,UMY);
                Real etot = *up0(i,j,k,UEDEN);
                Real ei   = *up0(i,j,k,UEINT);

                Real* AMREX_RESTRICT qrho  = qp0(i,j,k,QRHO);
                Real* AMREX_RESTRICT qu    = qp0(i,j,k,QU);
                Real* AMREX_RESTRICT qv    = qp0(i,j,k,QV);
                Real* AMREX_RESTRICT qw    = qp0(i,j,k,QW);
                Real* AMREX_RESTRICT qei   = qp0(i,j,k,QEINT);
                Real* AMREX_RESTRICT qpres = qp0(i,j,k,QPRES);
                Real* AMREX_RESTRICT qcs   = qp0(i,j,k,QCS);
                Real* AMREX_RESTRICT qtemp = qp0(i,j,k,QTEMP);

                (*qrho) = (rho > smallr) ? rho : smallr;
                Real rhoinv = 1.0/(*qrho);
                (*qu) = mx*rhoinv;
                (*qv) = my*rhoinv;
                (*qw) = mz*rhoinv;
                Real kineng = 0.5*(*qrho)*((*qu)*(*qu)+(*qv)*(*qv)+(*qw)*(*qw));
                (*qei) = (etot <= kineng) ? ei*rhoinv : (etot-kineng)*rhoinv;
                (*qpres) = (gamma-1.0)*ei;
                (*qpres) = ((*qpres) > smallp) ? (*qpres) : smallp;
                (*qcs) = std::sqrt(gamma*(*qpres)*rhoinv);
                (*qtemp) = 0.0;
            }
        }
    }

}

