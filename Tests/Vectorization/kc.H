
void SIMDORNOT(flux_to_dudt_c)
    (Box const& bx,
     FArrayBox& dudtfab,
     FArrayBox const& fxfab,
     FArrayBox const& fyfab,
     FArrayBox const& fzfab,
     Array<Real,AMREX_SPACEDIM> const& dxinv)
{
    const auto len = bx.length3d();
    const int ncomp = dudtfab.nComp();
    const auto dudt = dudtfab.view(bx);
    const auto fx   =   fxfab.view(bx);
    const auto fy   =   fyfab.view(bx);
    const auto fz   =   fzfab.view(bx);

    for (int n = 0; n < 5; ++n) {
        for         (int k = 0; k < len[2]; ++k) {
            for     (int j = 0; j < len[1]; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len[0]; ++i) {
                    dudt(i,j,k,n) = dxinv[0] * (fx(i,j,k,n) - fx(i+1,j,k,n))
                        +           dxinv[1] * (fy(i,j,k,n) - fy(i,j+1,k,n))
                        +           dxinv[2] * (fz(i,j,k,n) - fz(i,j,k+1,n));
                }
            }
        }
    }
}

void SIMDORNOT(ctoprim_c)
    (Box const& bx, FArrayBox const& ufab, FArrayBox & qfab)
{
    const auto len = bx.length3d();
    const auto u = ufab.view(bx);
    const auto q = qfab.view(bx);
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

    for         (int k = 0; k < len[2]; ++k) {
        for     (int j = 0; j < len[1]; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len[0]; ++i) {
                q(i,j,k,QRHO) = (u(i,j,k,URHO) > smallr)
                               ? u(i,j,k,URHO) : smallr;
                Real rhoinv = 1.0/q(i,j,k,QRHO);
                q(i,j,k,QU) = u(i,j,k,UMX)*rhoinv;
                q(i,j,k,QV) = u(i,j,k,UMY)*rhoinv;
                q(i,j,k,QW) = u(i,j,k,UMZ)*rhoinv;
                Real kineng = 0.5*q(i,j,k,QRHO)*(q(i,j,k,QU)*q(i,j,k,QU)
                                               + q(i,j,k,QV)*q(i,j,k,QV)
                                               + q(i,j,k,QW)*q(i,j,k,QW));
                q(i,j,k,QEINT) = (u(i,j,k,UEDEN) <= kineng)
                                ? u(i,j,k,UEINT)*rhoinv
                                : (u(i,j,k,UEDEN)-kineng)*rhoinv;
                q(i,j,k,QPRES) = (gamma-1.0)*u(i,j,k,UEINT);
                q(i,j,k,QPRES) = (q(i,j,k,QPRES) > smallp)
                                ? q(i,j,k,QPRES) : smallp;
                q(i,j,k,QCS) = std::sqrt(gamma*q(i,j,k,QPRES)*rhoinv);
                q(i,j,k,QTEMP) = 0.0;
            }
        }
    }
}

