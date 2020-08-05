#ifdef AMREX_USE_LAPACKE
#ifndef AMREX_GP_3D_H
#define AMREX_GP_3D_H
#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H> 
#include <vector>
#include <array> 
#include <cmath> 

// Contains the GP class for Steven Reeves' research into Gaussian Process based Models for AMR prolongation 

class GP
{
    public: 
    GP(const amrex::IntVect Ratio, const amrex::Real *del);
    ~GP(){}    

    // Member data
    amrex::IntVect r;
    amrex::Real dx[3]; 
    //
    // Eigen Values of Covariance Matrix 
    //
    amrex::Real lam[7] = {}; 
    //
    //  Eigen Vectors of Covariance Matrix
    //
    amrex::Real V[7][7] = {};
    //
    //  Weights to be applied for interpolation
    //
    std::vector<std::array<std::array<amrex::Real, 7>, 7>> ks;
    //
    //  Gammas needed for smoothness indication
    //
    std::vector<std::array<amrex::Real, 7>> gam;

    //Pointers for device functions
    amrex::Real *lamd; 
    amrex::Real *Vd; 
    amrex::Real *ksd; 
    amrex::Real *gamd; 
    amrex::Real l;
    amrex::Real sig;  

//Fill Managed Arrays  
    void h2mfill()
    {
        for(int i = 0; i < gam.size(); ++i){
            for(int j = 0; j < 7; ++j){
                gamd[i*7 + j] = gam[i][j];
                for(int k = 0; k < 7; ++k){
                    ksd[(i*7 + j)*7 + k] = ks[i][j][k]; 
                }
            }
        }
       for(int i = 0; i < 7; ++i){
            lamd[i] = lam[i]; 
            for(int j = 0; j < 7; ++j){
                Vd[i*7 + j] = V[j][i]; //Transpose so we load with fast index in the interpolater                
            }
        }
    }

//Free Temporary Mem
    void GP_finalize(){
#ifdef AMREX_USE_CUDA     
        cudaFree(gamd);
        cudaFree(ksd);
        cudaFree(lamd); 
        cudaFree(Vd);    
#else
        delete gamd;
        delete ksd; 
        delete lamd; 
        delete Vd;   
#endif 
    } 


// Linear Algebra Functions
    template<int n>
    AMREX_GPU_HOST_DEVICE
    inline
    static amrex::Real inner_prod(const amrex::Real x[n], const amrex::Real y[n])
    {
        amrex::Real result = 0.e0; 
        for(int i = 0; i < n; ++i) result += x[i]*y[i];
        return result;  
    }

    template<int n>
    inline
    static amrex::Real inner_prod(std::array<amrex::Real, n> &x, const amrex::Real y[n])
    {
        amrex::Real result = 0.e0; 
        for(int i = 0; i < n; ++i) result += x[i]*y[i];
        return result;  
    }


    template<int n> void CholeskyDecomp(amrex::Real (&K)[n][n]);
 

    template<int n> 
    void
    cholesky(amrex::Real (&b)[n], amrex::Real const K[n][n]); 
    
    template<int n> 
    void
    cholesky(std::array<amrex::Real, n> &b, amrex::Real const K[n][n]);

    void
    Decomp(amrex::Real (&K)[7][7], amrex::Real (&Kt)[25][25]); 

//GP functions! 
    template<class T> 
    amrex::Real sqrexp(const T x[3],const T y[3]);

    amrex::Real sqrexp(const std::array<amrex::Real, 3> x, const amrex::Real y[3]); 
    
    amrex::Real sqrexp2(const amrex::Real x[3], const amrex::Real y[3]); 
 
    inline
    amrex::Real cov1(const amrex::Real xc[3],
                     const amrex::Real yc[3],
                     const amrex::Real par)
    {
        amrex::Real x[3] = {xc[0]*dx[0], xc[1]*dx[1], xc[2]*dx[2]};
        amrex::Real y[3] = {yc[0]*dx[0], yc[1]*dx[1], yc[2]*dx[2]};

        amrex::Real rt2 = std::sqrt(2.e0);
        amrex::Real pi  = std::atan(1.e0)*4;
        amrex::Real dkh = (x[0] - y[0])/dx[0];
        amrex::Real arg1 = (dkh + 1.)/(rt2*par/dx[0]);
        amrex::Real arg2 = (dkh)/(rt2*par/dx[0]);
        amrex::Real arg3 = (dkh - 1.)/(rt2*par/dx[0]);
        amrex::Real cx = std::sqrt(pi)*(par*par/(dx[0]*dx[0]))*((arg1*std::erf(arg1)
                       + arg3*std::erf(arg3)) + 1.0/std::sqrt(pi)*
                         (std::exp(-arg1*arg1) + std::exp(-arg3*arg3))
                       - 2.0*(arg2*std::erf(arg2) + 1./std::sqrt(pi)*
                         std::exp(-arg2*arg2)));


        dkh = (x[1] - y[1])/dx[1];
        arg1 = (dkh + 1.)/(rt2*par/dx[1]);
        arg2 = (dkh)/(rt2*par/dx[1]);
        arg3 = (dkh - 1.)/(rt2*par/dx[1]);
        amrex::Real cy = std::sqrt(pi)*(par*par/(dx[1]*dx[1]))*((arg1*std::erf(arg1)
                       + arg3*std::erf(arg3)) + 1.0/std::sqrt(pi)*
                         (std::exp(-arg1*arg1) + std::exp(-arg3*arg3))
                       - 2.0*(arg2*std::erf(arg2) + 1./std::sqrt(pi)*
                         std::exp(-arg2*arg2)));

        dkh = (x[2] - y[2])/dx[2];
        arg1 = (dkh + 1.)/(rt2*par/dx[2]);
        arg2 = (dkh)/(rt2*par/dx[2]);
        arg3 = (dkh - 1.)/(rt2*par/dx[2]);
        amrex::Real cz = std::sqrt(pi)*(par*par/(dx[2]*dx[2]))*((arg1*std::erf(arg1)
                       + arg3*std::erf(arg3)) + 1.0/std::sqrt(pi)*
                         (std::exp(-arg1*arg1) + std::exp(-arg3*arg3))
                       - 2.0*(arg2*std::erf(arg2) + 1./std::sqrt(pi)*
                         std::exp(-arg2*arg2)));

        return cx*cy*cz;
    }


    inline
    amrex::Real cov2(const std::array<amrex::Real, 3> xc,
                     const amrex::Real yc[3])
    {
        amrex::Real x[3] = {xc[0]*dx[0], xc[1]*dx[1], xc[2]*dx[2]};
        amrex::Real y[3] = {yc[0]*dx[0], yc[1]*dx[1], yc[2]*dx[2]};
        amrex::Real dks = (x[0] - y[0])/dx[0];
        amrex::Real arg[4] = {static_cast<amrex::Real>( dks + (r[0] -1.0)/(2.0*r[0]) ),
                              static_cast<amrex::Real>( dks + (r[0] +1.0)/(2.0*r[0]) ),
                              static_cast<amrex::Real>( dks - (r[0] -1.0)/(2.0*r[0]) ),
                              static_cast<amrex::Real>( dks - (r[0] +1.0)/(2.0*r[0]) )};
        amrex::Real pi = std::atan(1.0)*4.0;
        amrex::Real cx = 0.0;
        amrex::Real iarg;
        for(int i = 0; i < 4; i++){
            iarg = arg[i]/(std::sqrt(2)*(l/dx[0]));
            cx += pow(-1,i+1)*(iarg*std::erf(iarg) + 1./(std::sqrt(pi))*std::exp(-iarg*iarg));
        }
        cx *= (l*l/(dx[0]*dx[0])); 

        dks = (x[1] - y[1])/dx[1];
        arg[0] = dks +(r[1] -1.)/(2.*r[1]);
        arg[1] = dks +(r[1] +1.)/(2.*r[1]);
        arg[2] = dks -(r[1] -1.)/(2.*r[1]);
        arg[3] = dks -(r[1] +1.)/(2.*r[1]);
        amrex::Real cy = 0.e0;
        for(int i = 0; i < 4; i++){
            iarg = arg[i]/(std::sqrt(2)*(l/dx[1])); 
            cy += pow(-1,i+1)*(iarg*std::erf(iarg) + 1./(std::sqrt(pi))*std::exp(-iarg*iarg));
        }
        cy *= (l*l/(dx[1]*dx[1])); 

        dks = (x[2] - y[2])/dx[2];
        arg[0] = dks +(r[2] -1.)/(2.*r[2]);
        arg[1] = dks +(r[2] +1.)/(2.*r[2]);
        arg[2] = dks -(r[2] -1.)/(2.*r[2]);
        arg[3] = dks -(r[2] +1.)/(2.*r[2]);
        amrex::Real cz = 0.e0;
        for(int i = 0; i < 4; i++){
            iarg = arg[i]/(std::sqrt(2)*(l/dx[2])); 
            cz += pow(-1,i+1)*(iarg*std::erf(iarg) + 1./(std::sqrt(pi))*std::exp(-iarg*iarg));
        }
        cz *= (l*l/(dx[2]*dx[2])); 
    
        return pi*std::sqrt(pi)*cx*cy*cz*r[0]*r[1]*r[2]; 
    }


    // Set up for the multi-sampled Weighted GP interpolation 
    // Build K makes the Coviarance Kernel Matrices for each Samples 
    // And for Total Stencil
    void GetK(amrex::Real (&K)[7][7], amrex::Real (&Ktot)[25][25]); 
    //
    // Get Weights builds k*Kinv for each stencil
    //
    void GetKs(const amrex::Real K[7][7]);
    //
    //  Get Weights for the LSQ RHS
    //
    void GetKtotks(const amrex::Real K1[25][25], std::vector<std::array<amrex::Real, 25>> &kt); 
    //
    // Get Gamma by solving a LSQ problem only need this once. 
    //
    void GetGamma(std::array<std::array<amrex::Real, 7>, 7> const& k,
                  std::array<amrex::Real,25> const &kt, std::array<amrex::Real,7> &ga); 
    //
    //  Get EigenVecs and EigenValues for smoothness indicators. 
    //  Will use the Shifted QR algorithm with deflation 
    //
    void GetEigen();
};


#endif
#endif 
