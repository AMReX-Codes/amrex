//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
//------------------------------------------------------------------------------------------------------------------------------
#ifndef    CA_KRYLOV_S
#define    CA_KRYLOV_S     4
#endif
//------------------------------------------------------------------------------------------------------------------------------
#include "matmul.c"
//------------------------------------------------------------------------------------------------------------------------------
// z[r] = alpha*A[r][c]*x[c]+beta*y[r]   // [row][col]
// z[r] = alpha*A[r][c]*x[c]+beta*y[r]   // [row][col]
#define gemv(z,alpha,A,x,beta,y,rows,cols)  {int r,c;double sum;for(r=0;r<(rows);r++){sum=0.0;for(c=0;c<(cols);c++){sum+=(A)[r][c]*(x)[c];}(z)[r]=(alpha)*sum+(beta)*(y)[r];}}
static inline void axpy(double * z, double alpha, double * x, double beta, double * y, int n){ // z[n] = alpha*x[n]+beta*y[n]
  int nn;
  for(nn=0;nn<n;nn++){
    z[nn] = alpha*x[nn] + beta*y[nn];
  }
}
static inline double vdotv(double * x, double * y, int n){ // x[n].y[n]
  int nn;
  double sum = 0.0;
  for(nn=0;nn<n;nn++){
    sum += x[nn]*y[nn];
  }
  return(sum);
}
static inline void zero(double * z, int n){ // z[n] = 0.0
  int nn;
  for(nn=0;nn<n;nn++){
    z[nn] = 0.0;
  }
}


//------------------------------------------------------------------------------------------------------------------------------
void CACG(level_type * level, int e_id, int R_id, double a, double b, double desired_reduction_in_norm){
  // based on Lauren Goodfriend, Yinghui Huang, and David Thorman's derivation in their Spring 2013 CS267 Report
  int    r0_id = VECTORS_RESERVED+0;
  int     r_id = VECTORS_RESERVED+1;
  int     p_id = VECTORS_RESERVED+2;
  int  PRrt_id = VECTORS_RESERVED+3;

  double  temp1[2*CA_KRYLOV_S+1];                                                             //
  double  temp2[2*CA_KRYLOV_S+1];                                                             //
  double  temp3[2*CA_KRYLOV_S+1];                                                             //
  double     aj[2*CA_KRYLOV_S+1];                                                             //
  double     cj[2*CA_KRYLOV_S+1];                                                             //
  double     ej[2*CA_KRYLOV_S+1];                                                             //
  double   Tpaj[2*CA_KRYLOV_S+1];                                                             //
  double     Tp[2*CA_KRYLOV_S+1][2*CA_KRYLOV_S+1];                                          // T'  indexed as [row][col]
  double      G[2*CA_KRYLOV_S+1][2*CA_KRYLOV_S+1];                                          // extracted from first 2*CA_KRYLOV_S+1 columns of Gg[][].  indexed as [row][col]
  double   Gbuf[(2*CA_KRYLOV_S+1)*(2*CA_KRYLOV_S+1)];                                       // buffer to hold the Gram-like matrix produced by matmul().  indexed as [row*(2*CA_KRYLOV_S+1) + col]
  int      PR[2*CA_KRYLOV_S+1];                                                               // vector_id's of the concatenation of the S+1 matrix powers of P, and the S matrix powers of R
  int *P = PR+              0;                                                                  // vector_id's of the S+1 Matrix Powers of P.  P[i] is the vector_id of A^i(p)
  int *R = PR+CA_KRYLOV_S+1;                                                                  // vector_id's of the S   Matrix Powers of R.  R[i] is the vector_id of A^i(r)

  int mMax=200;
  int m=0,n;
  int i,j,k;
  int CGFailed    = 0;
  int CGConverged = 0;

  double aj_dot_GTpaj,cj_dot_Gcj,alpha,cj_dot_Gcj_new,beta,L2_norm_of_r0,L2_norm_of_residual,delta;

  residual(level,r0_id,e_id,R_id,a,b);                                                            // r0[] = R_id[] - A(e_id)
  scale_vector(level,r_id,1.0,r0_id);                                                                // r[] = r0[]
  scale_vector(level, p_id,1.0,r0_id);                                                                // p[] = r0[]
  double norm_of_r0 = norm(level,r0_id);                                                          // the norm of the initial residual...
  if(norm_of_r0 == 0.0){CGConverged=1;}                                                          // entered CG with exact solution

  delta = dot(level,r_id,r0_id);                                                                   // delta = dot(r,r0)
  if(delta==0.0){CGConverged=1;}                                                                 // entered CG with exact solution (square of L2 norm of r_id)
  L2_norm_of_r0 = sqrt(delta);                                                                   // 



  // initialize Tp[][] ...
  for(i=0;i<2*CA_KRYLOV_S+1;i++)for(j=0;j<2*CA_KRYLOV_S+1;j++) Tp[i][j]=0;                  // zero Tp
  for(i=              0;i<  CA_KRYLOV_S  ;i++){ Tp[i+1][i]=1;}                                // monomial basis
  for(i=CA_KRYLOV_S+1;i<2*CA_KRYLOV_S  ;i++){ Tp[i+1][i]=1;}                                //

  for(i=0;i<2*CA_KRYLOV_S+1;i++){PR[i] = PRrt_id+i;}                                           // columns of PR map to the consecutive spare grids allocated for the bottom solver starting at PRrt_id


  while( (m<mMax) && (!CGFailed) && (!CGConverged) ){                                           // while(not done){
    zero(   aj,2*CA_KRYLOV_S+1);
    zero(   cj,2*CA_KRYLOV_S+1);
    zero(   ej,2*CA_KRYLOV_S+1);
    zero( Tpaj,2*CA_KRYLOV_S+1);
    zero(temp1,2*CA_KRYLOV_S+1);
    zero(temp2,2*CA_KRYLOV_S+1);
    zero(temp3,2*CA_KRYLOV_S+1);

    // Using the monomial basis, compute s+1 matrix powers on p[] and s matrix powers on r[] one power at a time
    //  (conventional approach applicable to CHOMBO and BoxLib)
    scale_vector(level,P[0],1.0, p_id);                                                             // P[0] = A^0p =  p_id
    for(n=1;n<CA_KRYLOV_S+1;n++){                                                             // naive way of calculating the monomial basis.
      apply_op(level,P[n],P[n-1],a,b);                                                          // P[n] = A(P[n-1]) = A^(n)p
    }
    scale_vector(level,R[0],1.0,r_id);                                                             // R[0] = A^0r = r_id
    for(n=1;n<CA_KRYLOV_S;n++){                                                               // naive way of calculating the monomial basis.
      apply_op(level,R[n],R[n-1],a,b);                                                          // R[n] = A(R[n-1]) = A^(n)r
    }


    // form G[][] and g[]
    level->CAKrylov_formations_of_G++;                                                         //   Record the number of times CACG formed G[][]
    matmul(level,Gbuf,PR,PR,2*CA_KRYLOV_S+1,2*CA_KRYLOV_S+1,1);                       // Compute Gbuf[][] = [P,R]^T * [P,R] (Matmul with grids but only one MPI_AllReduce)
    for(i=0,k=0;i<2*CA_KRYLOV_S+1;i++){                                                       // extract G[][] from Gbuf[]
    for(j=0    ;j<2*CA_KRYLOV_S+1;j++){G[i][j] = Gbuf[k++];}                                  // first 2*CA_KRYLOV_S+1 elements in each row go to G[][].
    }


    for(i=0;i<2*CA_KRYLOV_S+1;i++)aj[i]=0.0;aj[               0]=1.0;                         // initialized based on (???)
    for(i=0;i<2*CA_KRYLOV_S+1;i++)cj[i]=0.0;cj[CA_KRYLOV_S+1]=1.0;                          // initialized based on (???)
    for(i=0;i<2*CA_KRYLOV_S+1;i++)ej[i]=0.0;                                                  // initialized based on (???)

    for(n=0;n<CA_KRYLOV_S;n++){                                                               // for(n=0;n<CA_KRYLOV_S;n++){
      level->Krylov_iterations++;                                                               //   record number of inner-loop (j) iterations for comparison
      gemv( Tpaj,1.0,Tp,  aj,0.0, Tpaj,2*CA_KRYLOV_S+1,2*CA_KRYLOV_S+1);                  //               T'aj
      gemv(temp1,1.0, G,Tpaj,0.0,temp1,2*CA_KRYLOV_S+1,2*CA_KRYLOV_S+1);                  //    temp1[] = GT'aj
      gemv(temp2,1.0, G,  cj,0.0,temp2,2*CA_KRYLOV_S+1,2*CA_KRYLOV_S+1);                  //    temp2[] = Gcj
           aj_dot_GTpaj = vdotv(aj,temp1,2*CA_KRYLOV_S+1);                                    //   (aj,GT'aj)
             cj_dot_Gcj = vdotv(cj,temp2,2*CA_KRYLOV_S+1);                                    //   (cj,  Gcj)
      // FIX, can cj_dot_Gcj ever be zero ?
      if(aj_dot_GTpaj == 0.0){                                                                  //   pivot breakdown ???
        CGFailed=1;break;                                                                       //
      }                                                                                         //
      alpha = cj_dot_Gcj / aj_dot_GTpaj;                                                        //   alpha = (cj,Gcj) / (aj,GT'aj)
      if(isinf(alpha)){                                                                         //   alpha = big/tiny(overflow) = inf -> breakdown
        CGFailed=1;break;                                                                       // 
      }                                                                                         //
      axpy(   ej,1.0,ej,   alpha,   aj,2*CA_KRYLOV_S+1);                                    //   ej[] = ej[] + alpha*aj[]    
      axpy(   cj,1.0,cj,  -alpha, Tpaj,2*CA_KRYLOV_S+1);                                    //   cj[] = cj[] - alpha*T'*aj[]    
      gemv(temp2,1.0, G,  cj,0.0,temp2,2*CA_KRYLOV_S+1,2*CA_KRYLOV_S+1);                  //    temp2[] = Gcj
         cj_dot_Gcj_new = vdotv(cj,temp2,2*CA_KRYLOV_S+1);                                    //   (cj,  Gcj)
      // calculate the norm of the incremental residual (Saad's vector 'r') to check intra s-step convergence... == cj_dot_Gcj_new??
      L2_norm_of_residual = 0.0;if(cj_dot_Gcj_new>0)L2_norm_of_residual=sqrt(cj_dot_Gcj_new);   // finite precision can lead to the norm^2 being < 0 (Demmel says flush to 0.0)
      if(L2_norm_of_residual < desired_reduction_in_norm*L2_norm_of_r0){CGConverged=1;break;}   // terminate the inner n-loop
      if(cj_dot_Gcj_new == 0.0){                                                                //   Lanczos breakdown ???
        CGFailed=1;break;                                                                       //
      }                                                                                         //
      beta = cj_dot_Gcj_new / cj_dot_Gcj;                                                       // 
      if(isinf(beta)){CGFailed=1;break;}                                                        //   beta = inf?
      if(beta == 0.0){CGFailed=1;break;}                                                        //   beta = 0?  can't make further progress(?)
      axpy(   aj,1.0,cj,    beta,   aj,2*CA_KRYLOV_S+1);                                    //   cj[] = cj[] + beta*aj[]    

    }                                                                                           // inner n (j) loop

    // update iterates...
    for(i=0;i<2*CA_KRYLOV_S+1;i++){add_vectors(level,e_id,1.0,e_id,ej[i],PR[i]);}               // e_id[] = [P,R]ej + e_id[]
    if(!CGFailed && !CGConverged){                                                              // if we're done, then there is no point in updating these
                                   add_vectors(level,  p_id,0.0,  p_id,aj[0],PR[0]);                //    p[] = [P,R]aj
    for(i=1;i<2*CA_KRYLOV_S+1;i++){add_vectors(level,  p_id,1.0,  p_id,aj[i],PR[i]);}               //          ...
                                   add_vectors(level, r_id,0.0, r_id,cj[0],PR[0]);                //    r[] = [P,R]cj
    for(i=1;i<2*CA_KRYLOV_S+1;i++){add_vectors(level, r_id,1.0, r_id,cj[i],PR[i]);}               //          ...
    }
                              m+=CA_KRYLOV_S;                                                 //   m+=CA_KRYLOV_S;
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  }                                                                                             // } // outer m loop

}
