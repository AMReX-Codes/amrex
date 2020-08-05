#ifdef AMREX_USE_LAPACKE
#include <AMReX_GP_3D.H>
#include <AMReX_Gpu.H>
#include <lapacke.h> 

template<class T>
amrex::Real
GP::sqrexp(const T x[3], const T y[3])
{
    amrex::Real result = std::exp(-0.5*((x[0] - y[0])*(x[0] - y[0])*dx[0]*dx[0] + 
                                        (x[1] - y[1])*(x[1] - y[1])*dx[1]*dx[1] + 
                                        (x[2] - y[2])*(x[2] - y[2])*dx[2]*dx[2])/(l*l));
    return result;    
}
 
amrex::Real
GP::sqrexp(const std::array<amrex::Real, 3> x, const amrex::Real y[3])
{
    amrex::Real result = std::exp(-0.5*((x[0] - y[0])*(x[0] - y[0])*dx[0]*dx[0] + 
                                        (x[1] - y[1])*(x[1] - y[1])*dx[1]*dx[1] + 
                                        (x[2] - y[2])*(x[2] - y[2])*dx[2]*dx[2])/(l*l));
    return result;    
} 

amrex::Real
GP::sqrexp2(const amrex::Real x[3], const amrex::Real y[3])
{
    amrex::Real result = std::exp(-0.5*((x[0] - y[0])*(x[0] - y[0])*dx[0]*dx[0] + 
                                        (x[1] - y[1])*(x[1] - y[1])*dx[1]*dx[1] + 
                                        (x[2] - y[2])*(x[2] - y[2])*dx[2]*dx[2])/(sig*sig));
    return result;    
} 

    //Perfroms Cholesky Decomposition on covariance matrix K

GP::GP (const amrex::IntVect Ratio, const amrex::Real *del)
{
    BL_PROFILE_VAR("GP::GP()", gp_ctor); 
    D_DECL(dx[0] = del[0], dx[1] = del[1], dx[2] = del[2]); 
    r = Ratio;
//    if(dx[0] > 1./128.) l = 0.1; 
//    else 
    l = 12.*std::min(dx[0], std::min(dx[1], dx[2]));      
    sig = 3*std::min(dx[0], std::min(dx[1], dx[2])); 
    const int expfactor = D_TERM(r[0],*r[1],*r[2]);  
#ifdef AMREX_USE_CUDA     
    cudaMallocManaged(&gamd, expfactor*7*sizeof(amrex::Real));
    cudaMallocManaged(&ksd, expfactor*49*sizeof(amrex::Real));
    cudaMallocManaged(&lamd, 7*sizeof(amrex::Real)); 
    cudaMallocManaged(&Vd, 49*sizeof(amrex::Real)); 
#else
    gamd = new amrex::Real[expfactor*7];
    ksd = new amrex::Real[expfactor*49]; 
    Vd = new amrex::Real[49]; 
    lamd = new amrex::Real[7]; 
#endif

    amrex::Real K[7][7] = {}; //The same for every ratio;  
    amrex::Real Ktot[25][25] = {}; // The same for every ratio; 
    std::vector<std::array<amrex::Real, 25>> kt(r[0]*r[1]*r[2], std::array<amrex::Real, 25>{{0}});
         // First dim is rx*ry; 
    GetK(K, Ktot); // Builds Covariance Matrices of Base Sample and Extended Samples/stencils  
    GetEigen(); //Gets Eigenvalues and Vectors from K for use in the interpolation 
    Decomp(K, Ktot); //Decomposes K and Ktot into their Cholesky Versions
    ks.resize(r[0]*r[1]*r[2], std::array<std::array<amrex::Real, 7>, 7>()); 
    GetKs(K); 
    // K and Ktot are not actually necessary for the rest of the GP interpolation 
    // They are only used to construct the weights w = ks^T Kinv 
    // and gam = Rinv Q^T kt; 
    // ks, gam, lam and V are part of the class and will be used in the main interpolation routine. 
    gam.resize(r[0]*r[1]*r[2], std::array<amrex::Real, 7>()); 
    GetKtotks(Ktot, kt); 
    for(int i = 0; i < r[0]*r[1]*r[2]; ++i){
        GetGamma(ks[i], kt[i], gam[i]); //Gets the gamma's
    }
    h2mfill(); 
#ifdef AMREX_USE_CUDA 
    AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(lamd, 7*sizeof(amrex::Real),
                                              amrex::Gpu::Device::deviceId(),
                                              amrex::Gpu::gpuStream()));
    AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(Vd, 49*sizeof(amrex::Real),
                                              amrex::Gpu::Device::deviceId(),
                                              amrex::Gpu::gpuStream()));
    AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(ksd, expfactor*49*sizeof(amrex::Real),
                                              amrex::Gpu::Device::deviceId(),
                                              amrex::Gpu::gpuStream()));
    AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(gamd, expfactor*7*sizeof(amrex::Real),
                                              amrex::Gpu::Device::deviceId(),
                                              amrex::Gpu::gpuStream()));
#endif     
    BL_PROFILE_VAR_STOP(gp_ctor); 
}

template<int n>
void
GP::CholeskyDecomp(amrex::Real (&K)[n][n])
{
     for(int j = 0; j < n; ++j){
        for( int k = 0; k < j; ++k){
            K[j][j] -= (K[j][k]*K[j][k]);
        }
        K[j][j] = std::sqrt(K[j][j]);
        for(int i = j+1; i < n; ++i){
            for(int k = 0; k < j; ++k){
                K[i][j] -= K[i][k]*K[j][k];
            }
            K[i][j] /= K[j][j];
        }
    }
}

//Performs Cholesky Backsubstitution
template<int n> 
void
GP::cholesky(amrex::Real (&b)[n], amrex::Real const K[n][n])
{
    /* Forward sub Ly = b */
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < i; ++j) b[i] -= b[j]*K[i][j];
        b[i] /= K[i][i];
    }

    /* Back sub Ux = y */
    for(int i = n-1; i >= 0; --i){
        for(int j = i+1; j < n; ++j) b[i] -= K[j][i]*b[j];
        b[i] /= K[i][i];
    }
}

//Performs Cholesky Backsubstitution
template<int n> 
void
GP::cholesky(std::array<amrex::Real, n> &b, amrex::Real const K[n][n])
{
    /* Forward sub Ly = b */
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < i; ++j) b[i] -= b[j]*K[i][j];
        b[i] /= K[i][i];
    }

    /* Back sub Ux = y */
    for(int i = n-1; i >= 0; --i){
        for(int j = i+1; j < n; ++j) b[i] -= K[j][i]*b[j];
        b[i] /= K[i][i];
    }
}

//Builds the Covariance matrix K if uninitialized --> if(!init) GetK, weights etc.
//Four K totals to make the gammas.  
void
GP::GetK(amrex::Real (&K)[7][7], amrex::Real (&Ktot)[25][25])
{

    amrex::Real pnt[7][3] = {{ 0,  0, -1},
                             { 0, -1,  0}, 
                             {-1,  0,  0}, 
                             { 0,  0,  0}, 
                             { 1,  0,  0}, 
                             { 0,  1,  0},
                             { 0,  0,  1}}; 

//Small K
    for(int i = 0; i < 7; ++i)
        for(int j = i; j < 7; ++j){
            K[i][j] = cov1(pnt[i], pnt[j],l); 
            K[j][i] = K[i][j]; 
        }

    //Super K positions 
    amrex::Real spnt[25][3] =  {{ 0,  0, -2}, 
                                { 0, -1, -1},  
                                {-1,  0, -1}, 
                                { 0,  0, -1}, 
                                { 1,  0, -1},
                                { 0,  1, -1}, 
                                { 0, -2,  0}, 
                                {-1, -1,  0}, 
                                { 0, -1,  0}, 
                                { 1, -1,  0}, 
                                {-2,  0,  0}, 
                                {-1,  0,  0}, 
                                { 0,  0,  0}, 
                                { 1,  0,  0}, 
                                { 2,  0,  0}, 
                                {-1,  1,  0}, 
                                { 0,  1,  0}, 
                                { 1,  1,  0}, 
                                { 0,  2,  0}, 
                                { 0, -1,  1}, 
                                {-1,  0,  1}, 
                                { 0,  0,  1}, 
                                { 1,  0,  1}, 
                                { 0,  1,  1}, 
                                { 0,  0,  2}}; 

    for(int i = 0; i < 25; ++i)
        for(int j = i; j < 25; ++j){
            Ktot[i][j] = cov1(spnt[i], spnt[j],l); 
            Ktot[j][i] = Ktot[i][j]; 
        }
}

//We need to to the decomposition outside of the GetK routine so we can use K to get the 
//EigenVectors and Values. 

void
GP::Decomp(amrex::Real (&K)[7][7], amrex::Real (&Kt)[25][25])
{
    CholeskyDecomp<7>(K); 
    CholeskyDecomp<25>(Kt); 
}

//Use a Cholesky Decomposition to solve for k*K^-1 
//Inputs: K, outputs w = k*K^-1. 
//We need weights for each stencil. Therefore we'll have 5 arrays of 16 X 5 each. 

void 
GP::GetKs(const amrex::Real K[7][7])
{
    //Locations of new points relative to i,j 
    std::vector<std::array<amrex::Real,3>> pnt(r[0]*r[1]*r[2], std::array<amrex::Real, 3>()); 
    int id; 
    amrex::Real strt = -0.25;

    if(r[0] == 4 && r[1] == 4 && r[2] == 4){
        strt = -0.375; 
    }

    for(int k = 0; k < r[2]; ++k){
        for(int j = 0; j < r[1]; ++j){
            for(int i = 0; i < r[0]; ++i){ 
                id = i + r[0]*(j + r[1]*k); 
                pnt[id][0] = strt + 1.0/double(r[0])*i; 
                pnt[id][1] = strt + 1.0/double(r[1])*j; 
                pnt[id][2] = strt + 1.0/double(r[2])*k;
            }
        }
    }
    amrex::Real spnt[7][3] = {{ 0,  0, -1},
                              { 0, -1,  0}, 
                              {-1,  0,  0}, 
                              { 0,  0,  0}, 
                              { 1,  0,  0}, 
                              { 0,  1,  0},
                              { 0,  0,  1}}; 

    amrex::Real temp[3]; 
    //Build covariance vector between interpolant points and stencil 
     for(int i = 0; i < r[0]*r[1]*r[2]; ++i){
        for(int j = 0; j < 7; ++j){
            temp[0] = spnt[j][0], temp[1] = spnt[j][1], temp[2] = spnt[j][2] - 1.0; //sten_km
            ks[i][0][j] = cov2(pnt[i], temp);

            temp[0] = spnt[j][0], temp[1] = spnt[j][1] - 1.0, temp[2] = spnt[j][2]; //sten_jm
            ks[i][1][j] = cov2(pnt[i], temp);

            temp[0] = spnt[j][0] - 1.0, temp[1] = spnt[j][1], temp[2] = spnt[j][2]; //sten_im
            ks[i][2][j] = cov2(pnt[i], temp);

            ks[i][3][j] = cov2(pnt[i], spnt[j]); //sten_cen
    
            temp[0] = spnt[j][0] + 1.0, temp[1] = spnt[j][1], temp[2] = spnt[j][2]; //sten_ip
            ks[i][4][j] = cov2(pnt[i], temp);

            temp[0] = spnt[j][0], temp[1] = spnt[j][1] + 1.0, temp[2] = spnt[j][2];  //sten_jp
            ks[i][5][j] = cov2(pnt[i], temp);

            temp[0] = spnt[j][0], temp[1] = spnt[j][1], temp[2] = spnt[j][2] + 1.0;  //sten_kp
            ks[i][6][j] = cov2(pnt[i], temp);
        }
     //Backsubstitutes for k^TK^{-1} 
        for(int k = 0; k < 7; ++k)
            cholesky<7>(ks[i][k], K); 
   }
}

// Here we are using Kt to get the weights for the overdetermined  
// In this case, we will have 16 new points
// Therefore, we will need 16 b =  k*^T Ktot^(-1)
// K1 is already Choleskied  
void 
GP::GetKtotks(const amrex::Real K1[25][25], std::vector<std::array<amrex::Real, 25>> &kt)
{
    //Locations of new points relative to i,j 
    std::vector<std::array<amrex::Real,3>> pnt(r[0]*r[1]*r[2], std::array<amrex::Real,3>()); 
    int id; 
    amrex::Real strt = -0.25;

    if(r[0] == 4 && r[1] == 4 && r[2] == 4){
        strt = -0.375; 
    }

    for(int k = 0; k < r[2]; ++k){
        for(int j = 0; j < r[1]; ++j){
            for(int i = 0; i < r[0]; ++i){ 
                id = i + r[0]*(j + r[1]*k); 
                pnt[id][0] = strt + 1.0/double(r[0])*i; 
                pnt[id][1] = strt + 1.0/double(r[1])*j; 
                pnt[id][2] = strt + 1.0/double(r[2])*k;
            }
        }
    }
   
    //Super K positions 
    amrex::Real spnt[25][3] =  {{ 0,  0, -2}, 
                                { 0, -1, -1},  
                                {-1,  0, -1}, 
                                { 0,  0, -1}, 
                                { 1,  0, -1},
                                { 0,  1, -1}, 
                                { 0, -2,  0}, 
                                {-1, -1,  0}, 
                                { 0, -1,  0}, 
                                { 1, -1,  0}, 
                                {-2,  0,  0}, 
                                {-1,  0,  0}, 
                                { 0,  0,  0}, 
                                { 1,  0,  0}, 
                                { 2,  0,  0}, 
                                {-1,  1,  0}, 
                                { 0,  1,  0}, 
                                { 1,  1,  0}, 
                                { 0,  2,  0}, 
                                { 0, -1,  1}, 
                                {-1,  0,  1}, 
                                { 0,  0,  1}, 
                                { 1,  0,  1}, 
                                { 0,  1,  1}, 
                                { 0,  0,  2}}; 

    for(int i = 0; i < r[0]*r[1]*r[2]; i++){
       for (int j = 0; j < 25; j++){
            kt[i][j] = cov2(pnt[i], spnt[j]); 
       }
       cholesky<25>(kt[i], K1); 
    } 
}

//Each point will have its
//own set of gammas. 
//Use x = R^-1Q'b 
void
GP::GetGamma(std::array<std::array<amrex::Real, 7>, 7> const& k,
             std::array<amrex::Real, 25> const& kt, 
             std::array<amrex::Real, 7> &ga)
{
//Extended matrix Each column contains the vector of coviarances corresponding 
//to each sample (weno-like stencil)

                        //  km       jm       im       cen      ip       jp       kp 
    amrex::Real A[25*7] = { k[0][0], 0.e0   , 0.e0   , 0.e0   , 0.e0   , 0.e0   , 0.e0   ,  //i   j   k-2
                            k[0][1], k[1][0], 0.e0   , 0.e0   , 0.e0   , 0.e0   , 0.e0   ,  //i   j-1 k-1
                            k[0][2], 0.e0   , k[2][0], 0.e0   , 0.e0   , 0.e0   , 0.e0   ,  //i-1 j   k-1
                            k[0][3], 0.e0   , 0.e0   , k[3][0], 0.e0   , 0.e0   , 0.e0   ,  //i   j   k-1
                            k[0][4], 0.e0   , 0.e0   , 0.e0   , k[4][0], 0.e0   , 0.e0   ,  //i+1 j   k-1
                            k[0][5], 0.e0   , 0.e0   , 0.e0   , 0.e0   , k[5][0], 0.e0   ,  //i   j+1 k-1 
                            0.e0   , k[1][1], 0.e0   , 0.e0   , 0.e0   , 0.e0   , 0.e0   ,  //i   j-2 k
                            0.e0   , k[1][2], k[2][1], 0.e0   , 0.e0   , 0.e0   , 0.e0   ,  //i-1 j-1 k 
                            0.e0   , k[1][3], 0.e0   , k[3][1], 0.e0   , 0.e0   , 0.e0   ,  //i   j-1 k
                            0.e0   , k[1][4], 0.e0   , 0.e0   , k[4][1], 0.e0   , 0.e0   ,  //i+1 j-1 k 
                            0.e0   , 0.e0   , k[2][2], 0.e0   , 0.e0   , 0.e0   , 0.e0   ,  //i-2 j   k 
                            0.e0   , 0.e0   , k[2][3], k[3][2], 0.e0   , 0.e0   , 0.e0   ,  //i-1 j   k 
                            k[0][6], k[1][5], k[2][4], k[3][3], k[4][2], k[5][1], k[6][0],  //i   j   k 
                            0.e0   , 0.e0   , 0.e0   , k[3][4], k[4][3], 0.e0   , 0.e0   ,  //i+1 j   k 
                            0.e0   , 0.e0   , 0.e0   , 0.e0   , k[4][4], 0.e0   , 0.e0   ,  //i+2 j   k 
                            0.e0   , 0.e0   , k[2][5], 0.e0   , 0.e0   , k[5][2], 0.e0   ,  //i-1 j+1 k 
                            0.e0   , 0.e0   , 0.e0   , k[3][5], 0.e0   , k[5][3], 0.e0   ,  //i   j+1 k 
                            0.e0   , 0.e0   , 0.e0   , 0.e0   , k[4][5], k[5][4], 0.e0   ,  //i+1 j+1 k 
                            0.e0   , 0.e0   , 0.e0   , 0.e0   , 0.e0   , k[5][5], 0.e0   ,  //i   j+2 k 
                            0.e0   , k[1][6], 0.e0   , 0.e0   , 0.e0   , 0.e0   , k[6][1],  //i   j-1 k+1 
                            0.e0   , 0.e0   , k[2][6], 0.e0   , 0.e0   , 0.e0   , k[6][2],  //i-1 j   k+1
                            0.e0   , 0.e0   , 0.e0   , k[3][6], 0.e0   , 0.e0   , k[6][3],  //i   j   k+1
                            0.e0   , 0.e0   , 0.e0   , 0.e0   , k[4][6], 0.e0   , k[6][4],  //i+1 j   k+1 
                            0.e0   , 0.e0   , 0.e0   , 0.e0   , 0.e0   , k[5][6], k[6][5],  //i   j+1 k+1
                            0.e0   , 0.e0   , 0.e0   , 0.e0   , 0.e0   , 0.e0   , k[6][6]}; //i   j   k+2  
 
    int m = 25, n = 7, nrhs = 1; 
    int lda = 7, ldb = 1, lwork = -1, info; 
    double temp[25]; 
    for(int i = 0; i < 25; i++) temp[i] = kt[i]; 
    double workt;
    info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', m, n, nrhs, A, lda, temp, ldb);
    for(int i = 0; i < 7; ++i) ga[i] = temp[i];
}
 
extern "C"
{
 void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,
                double* w, double* work, int* lwork, int* info );
}
 

void 
GP::GetEigen()
{
    double A[49];
    amrex::Real pnt[7][3] = {{ 0,  0, -1}, // i   j    k-1
                             { 0, -1,  0}, // i   j-1  k
                             {-1,  0,  0}, // i-1 j    k
                             { 0,  0,  0}, // i   j    k
                             { 1,  0,  0}, // i+1 j    k
                             { 0,  1,  0}, // i   j+1  k 
                             { 0,  0,  1}};// i   j    k+1 

    for (int j = 0; j < 7; ++j){
        for(int i = j; i < 7; ++i){
             A[i + j*7] = cov1(pnt[i], pnt[j], sig); //this is K_sig
             A[j + 7*i] = A[i + j*7]; 
        }
    }

    int N = 7, lda = 7, info, lwork;
    double wkopt;
    double* work;
    lwork = -1;
    dsyev_( "Vectors", "Upper", &N, A, &lda, lam, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work = (double*)malloc( lwork*sizeof(double) );
    dsyev_( "Vectors", "Upper", &N, A, &lda, lam, work, &lwork, &info );
    for (int j = 0; j < 7; ++j){
        for(int i = 0; i < 7; ++i) V[i][j] =  A[i + j*7];
    }
    free(work);
}
#endif 
