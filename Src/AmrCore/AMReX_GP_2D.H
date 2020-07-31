#ifndef AMREX_GP_2D_H
#define AMREX_GP_2D_H
#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Array.H>
#include <AMReX_IntVect.H> 
#include <vector>
#include <array> 
#include <cmath> 

// Contains the GP class for Steven Reeves' research into Gaussian Process based Models for AMR prolongation 

class GP
{
    public: 
    GP(const amrex::IntVect Ratio, const amrex::Real *del);
    ~GP(){};      
    
    // Member data
    amrex::IntVect r; 
    amrex::Real dx[2];  
    //
    // Eigen Values of Covariance Matrix 
    //
    amrex::Real lam[5] = {}; 
    //
    //  Eigen Vectors of Covariance Matrix
    //
    amrex::Real V[5][5] = {};
    amrex::Real *lamd; 
    amrex::Real *Vd; 
    amrex::Real *ksd; 
    amrex::Real *gamd; 
    //
    //  Weights to be applied for interpolation
    //
    std::vector<std::array<std::array<amrex::Real, 5>, 5>> ks; 
    //
    //  Gammas needed for smoothness indication
    //
    std::vector<std::array<amrex::Real, 5>> gam;
    amrex::Real l;
    amrex::Real sig;  



//Fill Managed Arrays  
    void h2mfill()
    {
        for(int i = 0; i < gam.size(); ++i) 
            for(int j = 0; j < 5; ++j){
                gamd[i*5 + j] = gam[i][j]; 
                for(int k = 0; k < 5; ++k)
                    ksd[(i*5 + j)*5 + k] = ks[i][j][k]; 
            }
        for(int i = 0; i < 5; ++i){
            lamd[i] = lam[i]; 
            for(int j = 0; j < 5; ++j){
                Vd[i*5 + j] = V[j][i];                
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
    void
    cholesky(amrex::Real (&b)[n], amrex::Real const K[n][n]); 
    
    template<int n> 
    void
    cholesky(std::array<amrex::Real, n> &b, amrex::Real const K[n][n]);

    void
    Decomp(amrex::Real (&K)[5][5], amrex::Real (&Kt)[13][13]); 

//GP functions! 
    inline
    amrex::Real cov1(const amrex::Real xc[2], 
                     const amrex::Real yc[2], 
                     const amrex::Real par)
    {
        amrex::Real x[2] = {xc[0]*dx[0], xc[1]*dx[1]}; 
        amrex::Real y[2] = {yc[0]*dx[0], yc[1]*dx[1]}; 
        
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
        
        return cx*cy; 
    }

    inline
    amrex::Real cov2(const std::array<amrex::Real, 2> xc, 
                     const amrex::Real yc[2])
    {
        amrex::Real x[2] = {xc[0]*dx[0], xc[1]*dx[1]}; 
        amrex::Real y[2] = {yc[0]*dx[0], yc[1]*dx[1]}; 
        amrex::Real dks = (x[0] - y[0])/dx[0]; 
        amrex::Real arg[4] = { dks + (r[0] -1.0)/(2.0*r[0]), 
                               dks + (r[0] +1.0)/(2.0*r[0]), 
                               dks - (r[0] -1.0)/(2.0*r[0]), 
                               dks - (r[0] +1.0)/(2.0*r[0])};
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

        return pi*cx*cy*r[0]*r[1]; 
    }

    // Set up for the multi-sampled Weighted GP interpolation 
    // Build K makes the Coviarance Kernel Matrices for each Samples 
    // And for Total Stencil
    void GetK(amrex::Real (&K)[5][5], amrex::Real (&Ktot)[13][13]); 
    //
    // Get Weights builds k*Kinv for each stencil
    //
    void GetKs(const amrex::Real K[5][5]);
    //
    //  Get Weights for the LSQ RHS
    //
    void GetKtotks(const amrex::Real K1[13][13], std::vector<std::array<amrex::Real, 13>> &kt); 
    //
    // Get Gamma by solving a LSQ problem only need this once. 
    //
    void GetGamma(std::array<std::array<amrex::Real, 5>, 5> const& k,
                  std::array<amrex::Real,13> const &kt, std::array<amrex::Real,5> &ga); 
    //
    //  Get EigenVecs and EigenValues for smoothness indicators. 
    //  Will use the Shifted QR algorithm with deflation 
    //
    void GetEigen();
};



#endif 
