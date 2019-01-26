//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
// Lu = a*alpha[]*u[] - b*divergence( beta[]*gradient(u[]) )
//------------------------------------------------------------------------------------------------------------------------------
#ifndef DEFINES_H
#define DEFINES_H
//------------------------------------------------------------------------------------------------------------------------------
#define  VECTOR_TEMP         0 // 
#define  VECTOR_E            1 // error used in residual correction FMG
#define  VECTOR_R            2 // cell centered residual (f-Av)
//------------------------------------------------------------------------------------------------------------------------------
#define  VECTOR_F            3 // original right-hand side (Au=f), cell centered
#define  VECTOR_U            4 // numerical solution
#define  VECTOR_ALPHA        5 // cell centered coefficient
#define  VECTOR_BETA_I       6 // face centered coefficient (n.b. element 0 is the left face of the ghost zone element)
#define  VECTOR_BETA_J       7 // face centered coefficient (n.b. element 0 is the back face of the ghost zone element)
#define  VECTOR_BETA_K       8 // face centered coefficient (n.b. element 0 is the bottom face of the ghost zone element)
//------------------------------------------------------------------------------------------------------------------
#define  VECTOR_DINV         9 // cell centered relaxation parameter (e.g. inverse of the diagonal)
#define  VECTOR_L1INV       10 // cell centered relaxation parameter (e.g. inverse of the L1 norm of each row)
//------------------------------------------------------------------------------------------------------------------
#define VECTORS_RESERVED    11 // total number of vectors and the starting location for any auxillary bottom solver vectors
//------------------------------------------------------------------------------------------------------------------------------
#endif
