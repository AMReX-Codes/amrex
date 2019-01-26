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
void BiCGStab(level_type * level, int x_id, int R_id, double a, double b, double desired_reduction_in_norm){
  // Algorithm 7.7 in Iterative Methods for Sparse Linear Systems(Yousef Saad)
  // Algorithm 1 in Analysis and Practical use of Flexible BiCGStab (Jie Chen)
  int  r0_id = VECTORS_RESERVED+0;
  int   r_id = VECTORS_RESERVED+1;
  int   p_id = VECTORS_RESERVED+2;
  int   q_id = VECTORS_RESERVED+3; // q = D^{-1}p
  int   s_id = VECTORS_RESERVED+4;
  int   t_id = VECTORS_RESERVED+5; // t = D^{-1}s
  int  Ap_id = VECTORS_RESERVED+6;
  int  As_id = VECTORS_RESERVED+7;

  int jMax=200;
  int j=0;
  int BiCGStabFailed    = 0;
  int BiCGStabConverged = 0;
  residual(level,r0_id,x_id,R_id,a,b);                                          // r0[] = R_id[] - A(x_id)
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(level->must_subtract_mean == 1){
    double mean_of_r0 = mean(level,r0_id);
    shift_vector(level,r0_id,r0_id,-mean_of_r0);
  }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  scale_vector(level,r_id,1.0,r0_id);                                           // r[] = r0[]
  scale_vector(level,p_id,1.0,r0_id);                                           // p[] = r0[]
  double r_dot_r0 = dot(level,r_id,r0_id);                                      // r_dot_r0 = dot(r,r0)
  double norm_of_r0 = norm(level,r_id);                                         // the norm of the initial residual...
  if(r_dot_r0   == 0.0){BiCGStabConverged=1;}                                   // entered BiCGStab with exact solution
  if(norm_of_r0 == 0.0){BiCGStabConverged=1;}                                   // entered BiCGStab with exact solution
  while( (j<jMax) && (!BiCGStabFailed) && (!BiCGStabConverged) ){               // while(not done){
    j++;level->Krylov_iterations++;                                             //
    #ifdef KRYLOV_DIAGONAL_PRECONDITION                                         //
    mul_vectors(level,q_id,1.0,VECTOR_DINV,p_id);                               //   q[] = Dinv[]*p[]
    #else                                                                       //
    scale_vector(level,q_id,1.0,p_id);                                          //   q[] =        p[]
    #endif                                                                      //
    apply_op(level,Ap_id,q_id,a,b);                                             //   Ap[] = AM^{-1}(p)
    double Ap_dot_r0 = dot(level,Ap_id,r0_id);                                  //   Ap_dot_r0 = dot(Ap,r0)
    if(Ap_dot_r0 == 0.0){BiCGStabFailed=1;break;}                               //   pivot breakdown ???
    double alpha = r_dot_r0 / Ap_dot_r0;                                        //   alpha = r_dot_r0 / Ap_dot_r0
    if(isinf(alpha)){BiCGStabFailed=2;break;}                                   //   pivot breakdown ???
    add_vectors(level,x_id,1.0,x_id, alpha, q_id);                              //   x_id[] = x_id[] + alpha*q[]
    add_vectors(level,s_id,1.0,r_id,-alpha,Ap_id);                              //   s[]    = r[]    - alpha*Ap[]   (intermediate residual?)
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if(level->must_subtract_mean == 1){
      double mean_of_s = mean(level,s_id);
      shift_vector(level,s_id,s_id,-mean_of_s);
    }
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    double norm_of_s = norm(level,s_id);                                        //   FIX - redundant??  norm of intermediate residual
    if(norm_of_s == 0.0){BiCGStabConverged=1;break;}                            //   FIX - redundant??  if As_dot_As==0, then As must be 0 which implies s==0
    if(norm_of_s < desired_reduction_in_norm*norm_of_r0){BiCGStabConverged=1;break;}
    #ifdef KRYLOV_DIAGONAL_PRECONDITION                                         //
    mul_vectors(level,t_id,1.0,VECTOR_DINV,s_id);                               //   t[] = Dinv[]*s[]
    #else                                                                       //
    scale_vector(level,t_id,1.0,s_id);                                          //   t[] =        s[]
    #endif                                                                      //
    apply_op(level,As_id,t_id,a,b);                                             //   As = AM^{-1}(s)
    double As_dot_As = dot(level,As_id,As_id);                                  //   As_dot_As = dot(As,As)
    double As_dot_s  = dot(level,As_id, s_id);                                  //   As_dot_s  = dot(As, s)
    if(As_dot_As == 0.0){BiCGStabConverged=1;break;}                            //   converged ?
    double omega = As_dot_s / As_dot_As;                                        //   omega = As_dot_s / As_dot_As
    if(omega == 0.0){BiCGStabFailed=3;break;}                                   //   stabilization breakdown ???
    if(isinf(omega)){BiCGStabFailed=4;break;}                                   //   stabilization breakdown ???
    add_vectors(level,x_id,1.0,x_id, omega, t_id);                              //   x_id[] = x_id[] + omega*t[]
    add_vectors(level,r_id,1.0,s_id,-omega,As_id);                              //   r[]    = s[]    - omega*As[]  (recursively computed / updated residual)
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if(level->must_subtract_mean == 1){
      double mean_of_r = mean(level,r_id);
      shift_vector(level,r_id,r_id,-mean_of_r);
    }
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    double norm_of_r = norm(level,r_id);                                        //   norm of recursively computed residual (good enough??)
    if(norm_of_r == 0.0){BiCGStabConverged=1;break;}                            //
    if(norm_of_r < desired_reduction_in_norm*norm_of_r0){BiCGStabConverged=1;break;}
    double r_dot_r0_new = dot(level,r_id,r0_id);                                //   r_dot_r0_new = dot(r,r0)
    if(r_dot_r0_new == 0.0){BiCGStabFailed=5;break;}                            //   Lanczos breakdown ???
    double beta = (r_dot_r0_new/r_dot_r0) * (alpha/omega);                      //   beta = (r_dot_r0_new/r_dot_r0) * (alpha/omega)
    if(isinf(beta)){BiCGStabFailed=6;break;}                                    //   ???
    add_vectors(level,VECTOR_TEMP,1.0,p_id,-omega,      Ap_id);                 //   VECTOR_TEMP = (p[]-omega*Ap[])
    add_vectors(level,       p_id,1.0,r_id,  beta,VECTOR_TEMP);                 //   p[] = r[] + beta*(p[]-omega*Ap[])
    r_dot_r0 = r_dot_r0_new;                                                    //   r_dot_r0 = r_dot_r0_new   (save old r_dot_r0)
  }                                                                             // }
}
