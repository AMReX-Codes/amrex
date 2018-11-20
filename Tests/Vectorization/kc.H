
void SIMDORNOT(flux_to_dudt_c)
    (Box const& bx,
     FArrayBox& dudtfab,
     FArrayBox const& fxfab,
     FArrayBox const& fyfab,
     FArrayBox const& fzfab,
     Array<Real,AMREX_SPACEDIM> const& a_dxinv)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const int ncomp = dudtfab.nComp();
    const auto dudt = dudtfab.view(lo);
    const auto fx   =   fxfab.view(lo);
    const auto fy   =   fyfab.view(lo);
    const auto fz   =   fzfab.view(lo);
    const Real dxinv = a_dxinv[0];
    const Real dyinv = a_dxinv[1];
    const Real dzinv = a_dxinv[2];

    for (int n = 0; n < 7; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    dudt(i,j,k,n) = dxinv * (fx(i,j,k,n) - fx(i+1,j,k,n))
                        +           dyinv * (fy(i,j,k,n) - fy(i,j+1,k,n))
                        +           dzinv * (fz(i,j,k,n) - fz(i,j,k+1,n));
                }
            }
        }
    }
}

void SIMDORNOT(ctoprim_c)
    (Box const& bx, FArrayBox const& ufab, FArrayBox & qfab)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto u = ufab.view(lo);
    const auto q = qfab.view(lo);
    const Real smallr = 1.e-19;
    const Real smallp = 1.e-10;
    const Real gamma = 1.4;

    const int URHO  = 0;
    const int UMX   = 1;
    const int UMY   = 2;
    const int UMZ   = 3;
    const int UEDEN = 4;
    const int UEINT = 5;
    const int UTEMP = 6;

    const int QRHO   = 0;
    const int QU     = 1;
    const int QV     = 2;
    const int QW     = 3;
    const int QPRES  = 4;
    const int QCS    = 5;
    const int QEINT  = 6;
    const int QTEMP  = 7;

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                Real rho = u(i,j,k,URHO);
                rho = (rho > smallr) ? rho : smallr;
                Real rhoinv = 1.0/rho;
                Real ux = u(i,j,k,UMX)*rhoinv;
                Real uy = u(i,j,k,UMY)*rhoinv;
                Real uz = u(i,j,k,UMZ)*rhoinv;
                Real kineng = 0.5*rho*(ux*ux+uy*uy+uz*uz);
                Real ei = u(i,j,k,UEDEN) - kineng;
                if (ei <= 0.0) ei = u(i,j,k,UEINT);
                Real p = (gamma-1.0)*ei;
                p = (p > smallp) ? p : smallp;
                ei *= rhoinv;

                q(i,j,k,QRHO) = rho;
                q(i,j,k,QU) = ux;
                q(i,j,k,QV) = uy;
                q(i,j,k,QW) = uz;
                q(i,j,k,QEINT) = ei;
                q(i,j,k,QPRES) = p;
                q(i,j,k,QCS) = std::sqrt(gamma*p*rhoinv);
                q(i,j,k,QTEMP) = 0.0;
            }
        }
    }
}

