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
#define KRYLOV_DIAGONAL_PRECONDITION
//------------------------------------------------------------------------------------------------------------------------------
void CG(level_type * level, int x_id, int R_id, double a, double b, double desired_reduction_in_norm){
  // Algorithm 9.1 in Iterative Methods for Sparse Linear Systems(Yousef Saad)
  int  r0_id = VECTORS_RESERVED+0;
  int   r_id = VECTORS_RESERVED+1;
  int   p_id = VECTORS_RESERVED+2;
  int  Ap_id = VECTORS_RESERVED+3;
  int   z_id = VECTORS_RESERVED+4;

  int jMax=200;
  int j=0;
  int CGFailed    = 0;
  int CGConverged = 0;
  residual(level,r0_id,x_id,R_id,a,b);                                          // r0[] = R_id[] - A(x_id)
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(level->must_subtract_mean == 1){
    double mean_of_r0 = mean(level,r0_id);
    shift_vector(level,r0_id,r0_id,-mean_of_r0);
  }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  scale_vector(level,r_id,1.0,r0_id);                                           // r[] = r0[]
  #ifdef KRYLOV_DIAGONAL_PRECONDITION                                           //
  mul_vectors(level,z_id,1.0,VECTOR_DINV,r0_id);                                // z[] = Dinv[]*r0[]
  #else                                                                         //
  scale_vector(level,z_id,1.0,r0_id);                                           // z[] = I*r0[]
  #endif                                                                        //
  scale_vector(level,p_id,1.0,z_id);                                            // p[] = z[]
  double norm_of_r0 = norm(level,r_id);                                         // the norm of the initial residual...
  if(norm_of_r0 == 0.0){CGConverged=1;}                                         // entered CG with exact solution
  double r_dot_z = dot(level,r_id,z_id);                                        // r_dot_z = dot(r,z)
  while( (j<jMax) && (!CGFailed) && (!CGConverged) ){                           // while(not done){
    j++;level->Krylov_iterations++;                                             //
    apply_op(level,Ap_id,p_id,a,b);                                             //   Ap[] = A(p)
    double Ap_dot_p = dot(level,Ap_id,p_id);                                    //   Ap_dot_p = dot(Ap,p)
    if(Ap_dot_p == 0.0){CGFailed=1;break;}                                      //   pivot breakdown ???
    double alpha = r_dot_z / Ap_dot_p;                                          //   alpha = r_dot_z / Ap_dot_p
    if(isinf(alpha)){CGFailed=1;break;}                                         //   ???
    add_vectors(level,x_id,1.0,x_id, alpha,p_id );                              //   x_id[] = x_id[] + alpha*p[]
    add_vectors(level,r_id,1.0,r_id,-alpha,Ap_id);                              //   r[]    = r[]    - alpha*Ap[]   (intermediate residual?)
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if(level->must_subtract_mean == 1){
      double mean_of_r = mean(level,r_id);
      shift_vector(level,r_id,r_id,-mean_of_r);
    }
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    double norm_of_r = norm(level,r_id);                                        //   norm of intermediate residual
    if(norm_of_r == 0.0){CGConverged=1;break;}                                  //
    if(norm_of_r < desired_reduction_in_norm*norm_of_r0){CGConverged=1;break;}  //
    #ifdef KRYLOV_DIAGONAL_PRECONDITION                                         //
    mul_vectors(level,z_id,1.0,VECTOR_DINV,r_id);                               //   z[] = Dinv[]*r[]
    #else                                                                       //
    scale_vector(level,z_id,1.0,r_id);                                          //   z[] = I*r[]
    #endif                                                                      //
    double r_dot_z_new = dot(level,r_id,z_id);                                  //   r_dot_z_new = dot(r_{j+1},z_{j+1})
    if(r_dot_z_new == 0.0){CGFailed=1;break;}                                   //   Lanczos breakdown ???
    double beta = (r_dot_z_new/r_dot_z);                                        //   beta = (r_dot_z_new/r_dot_z)
    if(isinf(beta)){CGFailed=1;break;}                                          //   ???
    add_vectors(level,p_id,1.0,z_id,beta,p_id );                                //   p[] = z[] + beta*p[]
    r_dot_z = r_dot_z_new;                                                      //   r_dot_r = r_dot_r_new   (save old r_dot_r)
  }                                                                             // }
}
