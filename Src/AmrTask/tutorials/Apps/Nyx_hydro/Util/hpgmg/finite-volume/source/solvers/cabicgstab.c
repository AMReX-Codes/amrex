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
//#define KRYLOV_DIAGONAL_PRECONDITION
//------------------------------------------------------------------------------------------------------------------------------
#ifndef    CA_KRYLOV_TELESCOPING
#define    CA_KRYLOV_TELESCOPING
#endif
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
#ifdef CA_KRYLOV_TELESCOPING
void CABiCGStab(level_type * level, int e_id, int R_id, double a, double b, double desired_reduction_in_norm){
  // based on Erin Carson/Jim Demmel/Nick Knight's s-Step BiCGStab Algorithm 3.4
  // However, the formation of [P,R] is expensive ~ 4S+1 exchanges.  Moreover, formation of G[][] requires (4S+2)(4S+1) grid operations.
  //   When the required number of iterations is small, this overhead is large and can make the s-step version slower than vanilla BiCGStab
  //   Thus, this version is a telescoping s-step method that will start out with s=1, then do s=2, then s=4
  int    rt_id = VECTORS_RESERVED+0;
  int     r_id = VECTORS_RESERVED+1;
  int     p_id = VECTORS_RESERVED+2;
  int  PRrt_id = VECTORS_RESERVED+3;


  // note: CA_KRYLOV_S should be tiny (2-8?).  As such, 4*CA_KRYLOV_S+1 is also tiny (9-33).  Just allocate on the stack...
  double  temp1[4*CA_KRYLOV_S+1];                                                               //
  double  temp2[4*CA_KRYLOV_S+1];                                                               //
  double  temp3[4*CA_KRYLOV_S+1];                                                               //
  double     Tp[4*CA_KRYLOV_S+1][4*CA_KRYLOV_S+1];                                              // T'  indexed as [row][col]
  double    Tpp[4*CA_KRYLOV_S+1][4*CA_KRYLOV_S+1];                                              // T'' indexed as [row][col]
  double     aj[4*CA_KRYLOV_S+1];                                                               //
  double     cj[4*CA_KRYLOV_S+1];                                                               //
  double     ej[4*CA_KRYLOV_S+1];                                                               //
  double   Tpaj[4*CA_KRYLOV_S+1];                                                               //
  double   Tpcj[4*CA_KRYLOV_S+1];                                                               //
  double  Tppaj[4*CA_KRYLOV_S+1];                                                               //
  double      G[4*CA_KRYLOV_S+1][4*CA_KRYLOV_S+1];                                              // extracted from first 4*CA_KRYLOV_S+1 columns of Gg[][].  indexed as [row][col]
  double      g[4*CA_KRYLOV_S+1];                                                               // extracted from last [4*CA_KRYLOV_S+1] column of Gg[][].
  double    Gg[(4*CA_KRYLOV_S+1)*(4*CA_KRYLOV_S+2)];                                            // buffer to hold the Gram-like matrix produced by matmul().  indexed as [row*(4*CA_KRYLOV_S+2) + col]
  int      PRrt[4*CA_KRYLOV_S+2];                                                               // vector_id's of the concatenation of the 2S+1 matrix powers of P, 2S matrix powers of R, and rt

  int mMax=200;
  int m=0,n;
  int i,j,k;
  int BiCGStabFailed    = 0;
  int BiCGStabConverged = 0;
  double g_dot_Tpaj,alpha,omega_numerator,omega_denominator,omega,delta,delta_next,beta;
  double L2_norm_of_rt,L2_norm_of_residual,cj_dot_Gcj,L2_norm_of_s;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  residual(level,rt_id,e_id,R_id,a,b);                                                           // rt[] = R_id[] - A(e_id)... note, if DPC, then rt = R-AD^-1De
  scale_vector(level,r_id,1.0,rt_id);                                                               // r[] = rt[]
  scale_vector(level, p_id,1.0,rt_id);                                                               // p[] = rt[]
  double norm_of_rt = norm(level,rt_id);                                                         // the norm of the initial residual...
  #ifdef VERBOSE
  if(level->my_rank==0)ffprintf(stderr,stderr,"m=%8d, norm   =%0.20f\n",m,norm_of_rt);
  #endif
  if(norm_of_rt == 0.0){BiCGStabConverged=1;}                                                   // entered BiCGStab with exact solution
  delta = dot(level,r_id,rt_id);                                                                  // delta = dot(r,rt)
  if(delta==0.0){BiCGStabConverged=1;}                                                          // entered BiCGStab with exact solution (square of L2 norm of r_id)
  L2_norm_of_rt = sqrt(delta);

  int ca_krylov_s = 1;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  while( (m<mMax) && (!BiCGStabFailed) && (!BiCGStabConverged) ){                               // while(not done){
    zero(   aj,4*ca_krylov_s+1);                                                            //
    zero(   cj,4*ca_krylov_s+1);                                                            //
    zero(   ej,4*ca_krylov_s+1);                                                            //
    zero( Tpaj,4*ca_krylov_s+1);                                                            //
    zero( Tpcj,4*ca_krylov_s+1);                                                            //
    zero(Tppaj,4*ca_krylov_s+1);                                                            //
    zero(temp1,4*ca_krylov_s+1);                                                            //
    zero(temp2,4*ca_krylov_s+1);                                                            //
    zero(temp3,4*ca_krylov_s+1);                                                            //
 
    for(i=0;i<4*ca_krylov_s+1;i++)for(j=0;j<4*ca_krylov_s+1;j++) Tp[i][j]=0;                // initialize Tp[][] and Tpp[][] ...
    for(i=0;i<4*ca_krylov_s+1;i++)for(j=0;j<4*ca_krylov_s+1;j++)Tpp[i][j]=0;                //
    for(i=              0;i<2*ca_krylov_s  ;i++){ Tp[i+1][i]=1;}                            // monomial basis... Fixed (typo in SIAM paper)
    for(i=2*ca_krylov_s+1;i<4*ca_krylov_s  ;i++){ Tp[i+1][i]=1;}                            //
    for(i=              0;i<2*ca_krylov_s-1;i++){Tpp[i+2][i]=1;}                            //
    for(i=2*ca_krylov_s+1;i<4*ca_krylov_s-1;i++){Tpp[i+2][i]=1;}                            //

    for(i=0;i<4*ca_krylov_s+1;i++){PRrt[              i] = PRrt_id+i;}                       // columns of PRrt map to the consecutive spare grid indices starting at PRrt_id
                                   PRrt[4*ca_krylov_s+1] = rt_id;                            // last column or PRrt (r tilde) maps to rt
    int *P = PRrt+              0;                                                            // vector_id's of the 2S+1 Matrix Powers of P.  P[i] is the vector_id of A^i(p)
    int *R = PRrt+2*ca_krylov_s+1;                                                            // vector_id's of the 2S   Matrix Powers of R.  R[i] is the vector_id of A^i(r)

    // Using the monomial basis, compute 2s+1 matrix powers on p[] and 2s matrix powers on r[] one power at a time 
    // (conventional approach applicable to CHOMBO and BoxLib)
    scale_vector(level,P[0],1.0, p_id);                                                             // P[0] = A^0p =  p_id
    for(n=1;n<2*ca_krylov_s+1;n++){                                                           // naive way of calculating the monomial basis.
      #ifdef KRYLOV_DIAGONAL_PRECONDITION                                                             //
      mul_vectors(level, VECTOR_TEMP,1.0, VECTOR_DINV,P[n-1]);                                                //   temp[] = Dinv[]*P[n-1]
      apply_op(level,P[n], VECTOR_TEMP,a,b);                                                          //   P[n] = AD^{-1} VECTOR_TEMP = AD^{-1}P[n-1] = ((AD^{-1})^n)p
      #else                                                                                     //
      apply_op(level,P[n],P[n-1],a,b);                                                          //   P[n] = A(P[n-1]) = (A^n)p
      #endif                                                                                    //
    }
    scale_vector(level,R[0],1.0,r_id);                                                             // R[0] = A^0r = r_id
    for(n=1;n<2*ca_krylov_s;n++){                                                             // naive way of calculating the monomial basis.
      #ifdef KRYLOV_DIAGONAL_PRECONDITION                                                             //
      mul_vectors(level, VECTOR_TEMP,1.0, VECTOR_DINV,R[n-1]);                                                //   temp[] = Dinv[]*R[n-1]
      apply_op(level,R[n], VECTOR_TEMP,a,b);                                                          //   R[n] = AD^{-1} VECTOR_TEMP = AD^{-1}R[n-1]
      #else                                                                                     //
      apply_op(level,R[n],R[n-1],a,b);                                                          //   R[n] = A(R[n-1]) = (A^n)r
      #endif                                                                                    //
    }

    // Compute Gg[][] = [P,R]^T * [P,R,rt] (Matmul with grids with ghost zones but only one MPI_AllReduce)
    level->CAKrylov_formations_of_G++;                                                         //   Record the number of times CABiCGStab formed G[][]
    matmul(level,Gg,PRrt,PRrt,4*ca_krylov_s+1,4*ca_krylov_s+2,1);
    for(i=0,k=0;i<4*ca_krylov_s+1;i++){                                                       // extract G[][] and g[] from Gg[]
    for(j=0    ;j<4*ca_krylov_s+1;j++){G[i][j] = Gg[k++];}                                    // first 4*ca_krylov_s+1 elements in each row go to G[][].
                                         g[i]    = Gg[k++];                                     // last element in row goes to g[].
    }

    for(i=0;i<4*ca_krylov_s+1;i++)aj[i]=0.0;aj[              0]=1.0;                        // initialized based on (3.26)
    for(i=0;i<4*ca_krylov_s+1;i++)cj[i]=0.0;cj[2*ca_krylov_s+1]=1.0;                        // initialized based on (3.26)
    for(i=0;i<4*ca_krylov_s+1;i++)ej[i]=0.0;                                                  // initialized based on (3.26)

    for(n=0;n<ca_krylov_s;n++){                                                               // for(n=0;n<ca_krylov_s;n++){
      level->Krylov_iterations++;                                                               // record number of inner-loop (j) iterations for comparison
      gemv( Tpaj,   1.0, Tp,   aj,   0.0, Tpaj,4*ca_krylov_s+1,4*ca_krylov_s+1);          //    T'aj
      gemv( Tpcj,   1.0, Tp,   cj,   0.0, Tpcj,4*ca_krylov_s+1,4*ca_krylov_s+1);          //    T'cj
      gemv(Tppaj,   1.0,Tpp,   aj,   0.0,Tppaj,4*ca_krylov_s+1,4*ca_krylov_s+1);          //   T''aj
                       g_dot_Tpaj = vdotv(g,Tpaj,4*ca_krylov_s+1);                            // (g,T'aj)
      if(g_dot_Tpaj == 0.0){                                                                    // pivot breakdown ???
        #ifdef VERBOSE                                                                        //
        if(level->my_rank==0){ffprintf(stderr,stderr,"g_dot_Tpaj == 0.0\n");}                                   //
        #endif                                                                                  //
        BiCGStabFailed=1;break;                                                                 //
      }                                                                                         //
      alpha = delta / g_dot_Tpaj;                                                               // delta / (g,T'aj)
      if(isinf(alpha)){                                                                         // alpha = big/tiny(overflow) = inf -> breakdown
        #ifdef VERBOSE                                                                        //
        if(level->my_rank==0){ffprintf(stderr,stderr,"alpha == inf\n");}                                        // 
        #endif                                                                                  //
        BiCGStabFailed=1;break;                                                                 // 
      }                                                                                         // 
      #if 0                                                                                     // seems to have accuracy problems in finite precision...
      gemv(temp1,-alpha,  G, Tpaj,   0.0,temp1,4*ca_krylov_s+1,4*ca_krylov_s+1);          //  temp1[] =       - alpha*GT'aj
      gemv(temp1,   1.0,  G,   cj,   1.0,temp1,4*ca_krylov_s+1,4*ca_krylov_s+1);          //  temp1[] =   Gcj - alpha*GT'aj
      gemv(temp2,-alpha,  G,Tppaj,   0.0,temp2,4*ca_krylov_s+1,4*ca_krylov_s+1);          //  temp2[] =       − alpha*GT′′aj
      gemv(temp2,   1.0,  G, Tpcj,   1.0,temp2,4*ca_krylov_s+1,4*ca_krylov_s+1);          //  temp2[] = GT′cj − alpha*GT′′aj
      axpy(temp3,   1.0,     Tpcj,-alpha,Tppaj,4*ca_krylov_s+1);                            //  temp3[] =  T′cj − alpha*T′′aj
             omega_numerator = vdotv(temp3,temp1,4*ca_krylov_s+1);                            //  (temp3,temp1) = ( T'cj-alpha*T''aj ,   Gcj-alpha*GT'aj )
           omega_denominator = vdotv(temp3,temp2,4*ca_krylov_s+1);                            //  (temp3,temp2) = ( T′cj−alpha*T′′aj , GT′cj−alpha*GT′′aj )
      #else                                                                                     // better to change the order of operations Gx-Gy -> G(x-y) ...  (note, G is symmetric)
      axpy(temp1,   1.0,     Tpcj,-alpha,Tppaj,4*ca_krylov_s+1);                            //  temp1[] =  (T'cj - alpha*T''aj)
      gemv(temp2,   1.0,  G,temp1,   0.0,temp2,4*ca_krylov_s+1,4*ca_krylov_s+1);          //  temp2[] = G(T'cj - alpha*T''aj)
      axpy(temp3,   1.0,       cj,-alpha, Tpaj,4*ca_krylov_s+1);                            //  temp3[] =     cj - alpha*T'aj
             omega_numerator = vdotv(temp3,temp2,4*ca_krylov_s+1);                            //  (temp3,temp2) = ( (  cj - alpha*T'aj ) , G(T'cj - alpha*T''aj) )
           omega_denominator = vdotv(temp1,temp2,4*ca_krylov_s+1);                            //  (temp1,temp2) = ( (T'cj - alpha*T''aj) , G(T'cj - alpha*T''aj) )
      #endif                                                                                    // 
      // NOTE: omega_numerator/omega_denominator can be 0/x or 0/0, but should never be x/0
      // If omega_numerator==0, and ||s||==0, then convergence, x=x+alpha*aj
      // If omega_numerator==0, and ||s||!=0, then stabilization breakdown

      // !!! PARTIAL UPDATE OF ej MUST HAPPEN BEFORE THE CHECK ON OMEGA TO ENSURE FORWARD PROGRESS !!!
      axpy(   ej,1.0,ej,       alpha,   aj,4*ca_krylov_s+1);                                // ej[] = ej[] + alpha*aj[]    

      // calculate the norm of Saad's vector 's' to check intra s-step convergence...
      axpy(temp1,   1.0,       cj,-alpha, Tpaj,4*ca_krylov_s+1);                            //  temp1[] =   cj - alpha*T'aj
      gemv(temp2,   1.0,  G,temp1,   0.0,temp2,4*ca_krylov_s+1,4*ca_krylov_s+1);          //  temp2[] = G(cj - alpha*T'aj)
                                 L2_norm_of_s = vdotv(temp1,temp2,4*ca_krylov_s+1);           //  (temp1,temp2) = ( (cj - alpha*T'aj) , G(cj - alpha*T'aj) )  == square of L2 norm of s in exact arithmetic
      if(L2_norm_of_s<0)L2_norm_of_s=0;else L2_norm_of_s=sqrt(L2_norm_of_s);                    // finite precision can lead to the norm^2 being < 0 (Demmel says flush to 0.0)
      #ifdef VERBOSE                                                                          //
      if(level->my_rank==0){fprintf(stderr,"m=%8d, norm(s)=%0.20f\n",m+n,L2_norm_of_s);}                //
      #endif                                                                                    //
      if(L2_norm_of_s < desired_reduction_in_norm*L2_norm_of_rt){BiCGStabConverged=1;break;}    // terminate the inner n-loop


      if(omega_denominator == 0.0){                                                             // ??? breakdown
        #ifdef VERBOSE                                                                        //
        if(level->my_rank==0){if(omega_denominator == 0.0)fprintf(stderr,"omega_denominator == 0.0\n");}//
        #endif                                                                                  //
        BiCGStabFailed=1;break;                                                                 //
      }                                                                                         //
      omega = omega_numerator / omega_denominator;                                              // 
      if(isinf(omega)){                                                                         // omega = big/tiny(oveflow) = inf
        #ifdef VERBOSE                                                                        //
        if(level->my_rank==0){if(isinf(omega))fprintf(stderr,"omega == inf\n");}                        // 
        #endif                                                                                  //
        BiCGStabFailed=1;break;                                                                 //
      }                                                                                         //
      // !!! COMPLETE THE UPDATE OF ej & cj now that omega is known to be ok                    //
      axpy(   ej,1.0,ej,       omega,   cj,4*ca_krylov_s+1);                                // ej[] = ej[] + alpha*aj[] + omega*cj[]
      axpy(   ej,1.0,ej,-omega*alpha, Tpaj,4*ca_krylov_s+1);                                // ej[] = ej[] + alpha*aj[] + omega*cj[] - omega*alpha*T'aj[]
      axpy(   cj,1.0,cj,      -omega, Tpcj,4*ca_krylov_s+1);                                // cj[] = cj[] - omega*T'cj[]
      axpy(   cj,1.0,cj,      -alpha, Tpaj,4*ca_krylov_s+1);                                // cj[] = cj[] - omega*T'cj[] - alpha*T'aj[]
      axpy(   cj,1.0,cj, omega*alpha,Tppaj,4*ca_krylov_s+1);                                // cj[] = cj[] - omega*T'cj[] - alpha*T'aj[] + omega*alpha*T''aj[]


      // calculate the norm of the incremental residual (Saad's vector 'r') to check intra s-step convergence...
      gemv(temp1,   1.0,  G,   cj,   0.0,temp1,4*ca_krylov_s+1,4*ca_krylov_s+1);          // temp1[] = Gcj
                                       cj_dot_Gcj = vdotv(cj,temp1,4*ca_krylov_s+1);          // sqrt( (cj,Gcj) ) == L2 norm of the intermediate residual in exact arithmetic
      L2_norm_of_residual = 0.0;if(cj_dot_Gcj>0)L2_norm_of_residual=sqrt(cj_dot_Gcj);           // finite precision can lead to the norm^2 being < 0 (Demmel says flush to 0.0)
      #ifdef VERBOSE 
      if(level->my_rank==0){fprintf(stderr,"m=%8d, norm(r)=%0.20f (cj_dot_Gcj=%0.20e)\n",m+n,L2_norm_of_residual,cj_dot_Gcj);}
      #endif
      if(L2_norm_of_residual < desired_reduction_in_norm*L2_norm_of_rt){BiCGStabConverged=1;break;} // terminate the inner n-loop


      delta_next = vdotv( g,cj,4*ca_krylov_s+1);                                              // (g,cj)
      #ifdef VERBOSE                                                                          //
      if(level->my_rank==0){                                                                    //
        if(isinf(delta_next)     ){fprintf(stderr,"delta == inf\n");}                                   // delta = big/tiny(overflow) = inf
        if(delta_next      == 0.0){fprintf(stderr,"delta == 0.0\n");}                                   // Lanczos breakdown
        if(omega_numerator == 0.0){fprintf(stderr,"omega_numerator == 0.0\n");}                         // stabilization breakdown
        if(omega           == 0.0){fprintf(stderr,"omega == 0.0\n");}                                   // stabilization breakdown 
      }                                                                                         //
      #endif                                                                                    //
      if(isinf(delta_next)){BiCGStabFailed   =1;break;}                                         // delta = inf?
      if(delta_next  ==0.0){BiCGStabFailed   =1;break;}                                         // Lanczos breakdown...
      if(omega       ==0.0){BiCGStabFailed   =1;break;}                                         // stabilization breakdown 
      beta = (delta_next/delta)*(alpha/omega);                                                  // (delta_next/delta)*(alpha/omega)
      #ifdef VERBOSE                                                                          //
      if(level->my_rank==0){                                                                    //
        if(isinf(beta)           ){fprintf(stderr,"beta == inf\n");}                                    // beta = inf?
        if(beta            == 0.0){fprintf(stderr,"beta == 0.0\n");}                                    // beta = 0?  can't make further progress(?)
      }                                                                                         //
      #endif                                                                                    //
      if(isinf(beta)      ){BiCGStabFailed   =1;break;}                                         // beta = inf?
      if(beta       == 0.0){BiCGStabFailed   =1;break;}                                         // beta = 0?  can't make further progress(?)
      axpy(   aj,1.0,cj,        beta,   aj,4*ca_krylov_s+1);                                // aj[] = cj[] + beta*aj[]
      axpy(   aj,1.0,aj, -omega*beta, Tpaj,4*ca_krylov_s+1);                                // aj[] = cj[] + beta*aj[] - omega*beta*T'aj
      delta = delta_next;                                                                       // delta = delta_next

    }                                                                                           // inner n (j) loop

    // update iterates...
    for(i=0;i<4*ca_krylov_s+1;i++){add_vectors(level,e_id,1.0,e_id,ej[i],PRrt[i]);}             // e_id[] = [P,R]ej + e_id[]
    if(!BiCGStabFailed && !BiCGStabConverged){                                                  // if we're done, then there is no point in updating these
                                   add_vectors(level,  p_id,0.0,  p_id,aj[0],PRrt[0]);              //    p[] = [P,R]aj
    for(i=1;i<4*ca_krylov_s+1;i++){add_vectors(level,  p_id,1.0,  p_id,aj[i],PRrt[i]);}             //          ...
                                   add_vectors(level, r_id,0.0, r_id,cj[0],PRrt[0]);              //    r[] = [P,R]cj
    for(i=1;i<4*ca_krylov_s+1;i++){add_vectors(level, r_id,1.0, r_id,cj[i],PRrt[i]);}             //          ...
    }                                                                                           //
    m+=ca_krylov_s;                                                                           //   m+=ca_krylov_s;
    ca_krylov_s*=2;if(ca_krylov_s>CA_KRYLOV_S)ca_krylov_s=CA_KRYLOV_S;
  }                                                                                             // } // outer m loop
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #ifdef KRYLOV_DIAGONAL_PRECONDITION
  mul_vectors(level,e_id,1.0, VECTOR_DINV,e_id);                                                        //   e_id[] = Dinv[]*e_id[] // i.e. e = D^{-1}e'
  #endif

}
//------------------------------------------------------------------------------------------------------------------------------
#else // CA_KRYLOV_TELESCOPING =0
void CABiCGStab(level_type * level, int e_id, int R_id, double a, double b, double desired_reduction_in_norm){
  // based on Erin Carson/Jim Demmel/Nick Knight's s-Step BiCGStab Algorithm 3.4
  int    rt_id = VECTORS_RESERVED+0;
  int     r_id = VECTORS_RESERVED+1;
  int     p_id = VECTORS_RESERVED+2;
  int  PRrt_id = VECTORS_RESERVED+3;

  // note: CA_KRYLOV_S should be tiny (2-8?).  As such, 4*CA_KRYLOV_S+1 is also tiny (9-33).  Just allocate on the stack...
  double  temp1[4*CA_KRYLOV_S+1];                                               //
  double  temp2[4*CA_KRYLOV_S+1];                                               //
  double  temp3[4*CA_KRYLOV_S+1];                                               //
  double     Tp[4*CA_KRYLOV_S+1][4*CA_KRYLOV_S+1];                              // T'  indexed as [row][col]
  double    Tpp[4*CA_KRYLOV_S+1][4*CA_KRYLOV_S+1];                              // T'' indexed as [row][col]
  double     aj[4*CA_KRYLOV_S+1];                                               //
  double     cj[4*CA_KRYLOV_S+1];                                               //
  double     ej[4*CA_KRYLOV_S+1];                                               //
  double   Tpaj[4*CA_KRYLOV_S+1];                                               //
  double   Tpcj[4*CA_KRYLOV_S+1];                                               //
  double  Tppaj[4*CA_KRYLOV_S+1];                                               //
  double      G[4*CA_KRYLOV_S+1][4*CA_KRYLOV_S+1];                              // extracted from first 4*CA_KRYLOV_S+1 columns of Gg[][].  indexed as [row][col]
  double      g[4*CA_KRYLOV_S+1];                                               // extracted from last [4*CA_KRYLOV_S+1] column of Gg[][].
  double    Gg[(4*CA_KRYLOV_S+1)*(4*CA_KRYLOV_S+2)];                            // buffer to hold the Gram-like matrix produced by matmul().  indexed as [row*(4*CA_KRYLOV_S+2) + col]
  int      PRrt[4*CA_KRYLOV_S+2];                                               // vector_id's of the concatenation of the 2S+1 matrix powers of P, 2S matrix powers of R, and rt
  int *P = PRrt+                0;                                              // vector_id's of the 2S+1 Matrix Powers of P.  P[i] is the vector_id of A^i(p)
  int *R = PRrt+2*CA_KRYLOV_S+1;                                                // vector_id's of the 2S   Matrix Powers of R.  R[i] is the vector_id of A^i(r)

  int mMax=200;
  int m=0,n;
  int i,j,k;
  int BiCGStabFailed    = 0;
  int BiCGStabConverged = 0;
  double g_dot_Tpaj,alpha,omega_numerator,omega_denominator,omega,delta,delta_next,beta;
  double L2_norm_of_rt,L2_norm_of_residual,cj_dot_Gcj,L2_norm_of_s;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  residual(level,rt_id,e_id,R_id,a,b);                                           // rt[] = R_id[] - A(e_id)... note, if DPC, then rt = R-AD^-1De
  scale_vector(level,r_id,1.0,rt_id);                                               // r[] = rt[]
  scale_vector(level, p_id,1.0,rt_id);                                               // p[] = rt[]
  double norm_of_rt = norm(level,rt_id);                                         // the norm of the initial residual...
  #ifdef VERBOSE
  if(level->my_rank==0)fprintf(stderr,"m=%8d, norm   =%0.20f\n",m,norm_of_rt);
  #endif
  if(norm_of_rt == 0.0){BiCGStabConverged=1;}                                   // entered BiCGStab with exact solution
  delta = dot(level,r_id,rt_id);                                                  // delta = dot(r,rt)
  if(delta==0.0){BiCGStabConverged=1;}                                          // entered BiCGStab with exact solution (square of L2 norm of r_id)
  L2_norm_of_rt = sqrt(delta);

  int ca_krylov_s = CA_KRYLOV_S;                                              // by making this a variable, I prevent the compiler from optimizing more than the telescoping version, thus preserving a bit-identcal result

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for(i=0;i<4*ca_krylov_s+1;i++)for(j=0;j<4*ca_krylov_s+1;j++) Tp[i][j]=0;  // initialize Tp[][] and Tpp[][] ...
  for(i=0;i<4*ca_krylov_s+1;i++)for(j=0;j<4*ca_krylov_s+1;j++)Tpp[i][j]=0;  //
  for(i=              0;i<2*ca_krylov_s  ;i++){ Tp[i+1][i]=1;}              // monomial basis... Fixed (typo in SIAM paper)
  for(i=2*ca_krylov_s+1;i<4*ca_krylov_s  ;i++){ Tp[i+1][i]=1;}              //
  for(i=              0;i<2*ca_krylov_s-1;i++){Tpp[i+2][i]=1;}              //
  for(i=2*ca_krylov_s+1;i<4*ca_krylov_s-1;i++){Tpp[i+2][i]=1;}              //

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for(i=0;i<4*ca_krylov_s+1;i++){PRrt[              i] = PRrt_id+i;}         // columns of PRrt map to the consecutive spare grid indices starting at PRrt_id
                                 PRrt[4*ca_krylov_s+1] = rt_id;              // last column or PRrt (r tilde) maps to rt

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  while( (m<mMax) && (!BiCGStabFailed) && (!BiCGStabConverged) ){               // while(not done){
    zero(   aj,4*ca_krylov_s+1);                                            //
    zero(   cj,4*ca_krylov_s+1);                                            //
    zero(   ej,4*ca_krylov_s+1);                                            //
    zero( Tpaj,4*ca_krylov_s+1);                                            //
    zero( Tpcj,4*ca_krylov_s+1);                                            //
    zero(Tppaj,4*ca_krylov_s+1);                                            //
    zero(temp1,4*ca_krylov_s+1);                                            //
    zero(temp2,4*ca_krylov_s+1);                                            //
    zero(temp3,4*ca_krylov_s+1);                                            //

    // Using the monomial basis, compute 2s+1 matrix powers on p[] and 2s matrix powers on r[] one power at a time 
    // (conventional approach applicable to CHOMBO and BoxLib)
    scale_vector(level,P[0],1.0, p_id);                                             // P[0] = A^0p =  p_id
    for(n=1;n<2*ca_krylov_s+1;n++){                                           // naive way of calculating the monomial basis.
      #ifdef KRYLOV_DIAGONAL_PRECONDITION                                             //
      mul_vectors(level, VECTOR_TEMP,1.0, VECTOR_DINV,P[n-1]);                           //   temp[] = Dinv[]*P[n-1]
      apply_op(level,P[n], VECTOR_TEMP,a,b);                                          //   P[n] = AD^{-1} VECTOR_TEMP = AD^{-1}P[n-1] = ((AD^{-1})^n)p
      #else                                                                     //
      apply_op(level,P[n],P[n-1],a,b);                                          //   P[n] = A(P[n-1]) = (A^n)p
      #endif                                                                    //
    }
    scale_vector(level,R[0],1.0,r_id);                                             // R[0] = A^0r = r_id
    for(n=1;n<2*ca_krylov_s;n++){                                             // naive way of calculating the monomial basis.
      #ifdef KRYLOV_DIAGONAL_PRECONDITION                                             //
      mul_vectors(level, VECTOR_TEMP,1.0, VECTOR_DINV,R[n-1]);                                //   temp[] = Dinv[]*R[n-1]
      apply_op(level,R[n], VECTOR_TEMP,a,b);                                          //   R[n] = AD^{-1} VECTOR_TEMP = AD^{-1}R[n-1]
      #else                                                                     //
      apply_op(level,R[n],R[n-1],a,b);                                          //   R[n] = A(R[n-1]) = (A^n)r
      #endif                                                                    //
    }

    // Compute Gg[][] = [P,R]^T * [P,R,rt] (Matmul with grids with ghost zones but only one MPI_AllReduce)
    level->CAKrylov_formations_of_G++;                                                         //   Record the number of times CABiCGStab formed G[][]
    matmul(level,Gg,PRrt,PRrt,4*ca_krylov_s+1,4*ca_krylov_s+2,1);
    for(i=0,k=0;i<4*ca_krylov_s+1;i++){                                                       // extract G[][] and g[] from Gg[]
    for(j=0    ;j<4*ca_krylov_s+1;j++){G[i][j] = Gg[k++];}                                    // first 4*ca_krylov_s+1 elements in each row go to G[][].
                                         g[i]    = Gg[k++];                                     // last element in row goes to g[].
    }

    for(i=0;i<4*ca_krylov_s+1;i++)aj[i]=0.0;aj[              0]=1.0;                        // initialized based on (3.26)
    for(i=0;i<4*ca_krylov_s+1;i++)cj[i]=0.0;cj[2*ca_krylov_s+1]=1.0;                        // initialized based on (3.26)
    for(i=0;i<4*ca_krylov_s+1;i++)ej[i]=0.0;                                                  // initialized based on (3.26)

    for(n=0;n<ca_krylov_s;n++){                                                               // for(n=0;n<ca_krylov_s;n++){
      level->Krylov_iterations++;                                                               // record number of inner-loop (j) iterations for comparison
      gemv( Tpaj,   1.0, Tp,   aj,   0.0, Tpaj,4*ca_krylov_s+1,4*ca_krylov_s+1);          //    T'aj
      gemv( Tpcj,   1.0, Tp,   cj,   0.0, Tpcj,4*ca_krylov_s+1,4*ca_krylov_s+1);          //    T'cj
      gemv(Tppaj,   1.0,Tpp,   aj,   0.0,Tppaj,4*ca_krylov_s+1,4*ca_krylov_s+1);          //   T''aj
                       g_dot_Tpaj = vdotv(g,Tpaj,4*ca_krylov_s+1);                            // (g,T'aj)
      if(g_dot_Tpaj == 0.0){                                                                    // pivot breakdown ???
        #ifdef VERBOSE                                                                        //
        if(level->my_rank==0){fprintf(stderr,"g_dot_Tpaj == 0.0\n");}                                   //
        #endif                                                                                  //
        BiCGStabFailed=1;break;                                                                 //
      }                                                                                         //
      alpha = delta / g_dot_Tpaj;                                                               // delta / (g,T'aj)
      if(isinf(alpha)){                                                                         // alpha = big/tiny(overflow) = inf -> breakdown
        #ifdef VERBOSE                                                                        //
        if(level->my_rank==0){fprintf(stderr,"alpha == inf\n");}                                        // 
        #endif                                                                                  //
        BiCGStabFailed=1;break;                                                                 // 
      }                                                                                         // 
      #if 0                                                                                     // seems to have accuracy problems in finite precision...
      gemv(temp1,-alpha,  G, Tpaj,   0.0,temp1,4*ca_krylov_s+1,4*ca_krylov_s+1);          //  temp1[] =       - alpha*GT'aj
      gemv(temp1,   1.0,  G,   cj,   1.0,temp1,4*ca_krylov_s+1,4*ca_krylov_s+1);          //  temp1[] =   Gcj - alpha*GT'aj
      gemv(temp2,-alpha,  G,Tppaj,   0.0,temp2,4*ca_krylov_s+1,4*ca_krylov_s+1);          //  temp2[] =       − alpha*GT′′aj
      gemv(temp2,   1.0,  G, Tpcj,   1.0,temp2,4*ca_krylov_s+1,4*ca_krylov_s+1);          //  temp2[] = GT′cj − alpha*GT′′aj
      axpy(temp3,   1.0,     Tpcj,-alpha,Tppaj,4*ca_krylov_s+1);                            //  temp3[] =  T′cj − alpha*T′′aj
             omega_numerator = vdotv(temp3,temp1,4*ca_krylov_s+1);                            //  (temp3,temp1) = ( T'cj-alpha*T''aj ,   Gcj-alpha*GT'aj )
           omega_denominator = vdotv(temp3,temp2,4*ca_krylov_s+1);                            //  (temp3,temp2) = ( T′cj−alpha*T′′aj , GT′cj−alpha*GT′′aj )
      #else                                                                                     // better to change the order of operations Gx-Gy -> G(x-y) ...  (note, G is symmetric)
      axpy(temp1,   1.0,     Tpcj,-alpha,Tppaj,4*ca_krylov_s+1);                            //  temp1[] =  (T'cj - alpha*T''aj)
      gemv(temp2,   1.0,  G,temp1,   0.0,temp2,4*ca_krylov_s+1,4*ca_krylov_s+1);          //  temp2[] = G(T'cj - alpha*T''aj)
      axpy(temp3,   1.0,       cj,-alpha, Tpaj,4*ca_krylov_s+1);                            //  temp3[] =     cj - alpha*T'aj
             omega_numerator = vdotv(temp3,temp2,4*ca_krylov_s+1);                            //  (temp3,temp2) = ( (  cj - alpha*T'aj ) , G(T'cj - alpha*T''aj) )
           omega_denominator = vdotv(temp1,temp2,4*ca_krylov_s+1);                            //  (temp1,temp2) = ( (T'cj - alpha*T''aj) , G(T'cj - alpha*T''aj) )
      #endif                                                                                    // 
      // NOTE: omega_numerator/omega_denominator can be 0/x or 0/0, but should never be x/0
      // If omega_numerator==0, and ||s||==0, then convergence, x=x+alpha*aj
      // If omega_numerator==0, and ||s||!=0, then stabilization breakdown

      // !!! PARTIAL UPDATE OF ej MUST HAPPEN BEFORE THE CHECK ON OMEGA TO ENSURE FORWARD PROGRESS !!!
      axpy(   ej,1.0,ej,       alpha,   aj,4*ca_krylov_s+1);                                // ej[] = ej[] + alpha*aj[]    

      // calculate the norm of Saad's vector 's' to check intra s-step convergence...
      axpy(temp1,   1.0,       cj,-alpha, Tpaj,4*ca_krylov_s+1);                            //  temp1[] =   cj - alpha*T'aj
      gemv(temp2,   1.0,  G,temp1,   0.0,temp2,4*ca_krylov_s+1,4*ca_krylov_s+1);          //  temp2[] = G(cj - alpha*T'aj)
                                 L2_norm_of_s = vdotv(temp1,temp2,4*ca_krylov_s+1);           //  (temp1,temp2) = ( (cj - alpha*T'aj) , G(cj - alpha*T'aj) )  == square of L2 norm of s in exact arithmetic
      if(L2_norm_of_s<0)L2_norm_of_s=0;else L2_norm_of_s=sqrt(L2_norm_of_s);                    // finite precision can lead to the norm^2 being < 0 (Demmel says flush to 0.0)
      #ifdef VERBOSE                                                                          //
      if(level->my_rank==0){fprintf(stderr,"m=%8d, norm(s)=%0.20f\n",m+n,L2_norm_of_s);}                //
      #endif                                                                                    //
      if(L2_norm_of_s < desired_reduction_in_norm*L2_norm_of_rt){BiCGStabConverged=1;break;}    // terminate the inner n-loop


      if(omega_denominator == 0.0){                                                             // ??? breakdown
        #ifdef VERBOSE                                                                        //
        if(level->my_rank==0){if(omega_denominator == 0.0)fprintf(stderr,"omega_denominator == 0.0\n");}//
        #endif                                                                                  //
        BiCGStabFailed=1;break;                                                                 //
      }                                                                                         //
      omega = omega_numerator / omega_denominator;                                              // 
      if(isinf(omega)){                                                                         // omega = big/tiny(oveflow) = inf
        #ifdef VERBOSE                                                                        //
        if(level->my_rank==0){if(isinf(omega))fprintf(stderr,"omega == inf\n");}                        // 
        #endif                                                                                  //
        BiCGStabFailed=1;break;                                                                 //
      }                                                                                         //
      // !!! COMPLETE THE UPDATE OF ej & cj now that omega is known to be ok                    //
      axpy(   ej,1.0,ej,       omega,   cj,4*ca_krylov_s+1);                                // ej[] = ej[] + alpha*aj[] + omega*cj[]
      axpy(   ej,1.0,ej,-omega*alpha, Tpaj,4*ca_krylov_s+1);                                // ej[] = ej[] + alpha*aj[] + omega*cj[] - omega*alpha*T'aj[]
      axpy(   cj,1.0,cj,      -omega, Tpcj,4*ca_krylov_s+1);                                // cj[] = cj[] - omega*T'cj[]
      axpy(   cj,1.0,cj,      -alpha, Tpaj,4*ca_krylov_s+1);                                // cj[] = cj[] - omega*T'cj[] - alpha*T'aj[]
      axpy(   cj,1.0,cj, omega*alpha,Tppaj,4*ca_krylov_s+1);                                // cj[] = cj[] - omega*T'cj[] - alpha*T'aj[] + omega*alpha*T''aj[]


      // calculate the norm of the incremental residual (Saad's vector 'r') to check intra s-step convergence...
      gemv(temp1,   1.0,  G,   cj,   0.0,temp1,4*ca_krylov_s+1,4*ca_krylov_s+1);          // temp1[] = Gcj
                                       cj_dot_Gcj = vdotv(cj,temp1,4*ca_krylov_s+1);          // sqrt( (cj,Gcj) ) == L2 norm of the intermediate residual in exact arithmetic
      L2_norm_of_residual = 0.0;if(cj_dot_Gcj>0)L2_norm_of_residual=sqrt(cj_dot_Gcj);           // finite precision can lead to the norm^2 being < 0 (Demmel says flush to 0.0)
      #ifdef VERBOSE 
      if(level->my_rank==0){fprintf(stderr,"m=%8d, norm(r)=%0.20f (cj_dot_Gcj=%0.20e)\n",m+n,L2_norm_of_residual,cj_dot_Gcj);}
      #endif
      if(L2_norm_of_residual < desired_reduction_in_norm*L2_norm_of_rt){BiCGStabConverged=1;break;} // terminate the inner n-loop


      delta_next = vdotv( g,cj,4*ca_krylov_s+1);                                              // (g,cj)
      #ifdef VERBOSE                                                                          //
      if(level->my_rank==0){                                                                    //
        if(isinf(delta_next)     ){fprintf(stderr,"delta == inf\n");}                                   // delta = big/tiny(overflow) = inf
        if(delta_next      == 0.0){fprintf(stderr,"delta == 0.0\n");}                                   // Lanczos breakdown
        if(omega_numerator == 0.0){fprintf(stderr,"omega_numerator == 0.0\n");}                         // stabilization breakdown
        if(omega           == 0.0){fprintf(stderr,"omega == 0.0\n");}                                   // stabilization breakdown 
      }                                                                                         //
      #endif                                                                                    //
      if(isinf(delta_next)){BiCGStabFailed   =1;break;}                                         // delta = inf?
      if(delta_next  ==0.0){BiCGStabFailed   =1;break;}                                         // Lanczos breakdown...
      if(omega       ==0.0){BiCGStabFailed   =1;break;}                                         // stabilization breakdown 
      beta = (delta_next/delta)*(alpha/omega);                                                  // (delta_next/delta)*(alpha/omega)
      #ifdef VERBOSE                                                                          //
      if(level->my_rank==0){                                                                    //
        if(isinf(beta)           ){fprintf(stderr,"beta == inf\n");}                                    // beta = inf?
        if(beta            == 0.0){fprintf(stderr,"beta == 0.0\n");}                                    // beta = 0?  can't make further progress(?)
      }                                                                                         //
      #endif                                                                                    //
      if(isinf(beta)      ){BiCGStabFailed   =1;break;}                                         // beta = inf?
      if(beta       == 0.0){BiCGStabFailed   =1;break;}                                         // beta = 0?  can't make further progress(?)
      axpy(   aj,1.0,cj,        beta,   aj,4*ca_krylov_s+1);                                // aj[] = cj[] + beta*aj[]
      axpy(   aj,1.0,aj, -omega*beta, Tpaj,4*ca_krylov_s+1);                                // aj[] = cj[] + beta*aj[] - omega*beta*T'aj
      delta = delta_next;                                                                       // delta = delta_next

    }                                                                                           // inner n (j) loop

    // update iterates...
    for(i=0;i<4*ca_krylov_s+1;i++){add_vectors(level,e_id,1.0,e_id,ej[i],PRrt[i]);}             // e_id[] = [P,R]ej + e_id[]
    if(!BiCGStabFailed && !BiCGStabConverged){                                                  // if we're done, then there is no point in updating these
                                   add_vectors(level,  p_id,0.0,  p_id,aj[0],PRrt[0]);              //    p[] = [P,R]aj
    for(i=1;i<4*ca_krylov_s+1;i++){add_vectors(level,  p_id,1.0,  p_id,aj[i],PRrt[i]);}             //          ...
                                   add_vectors(level, r_id,0.0, r_id,cj[0],PRrt[0]);              //    r[] = [P,R]cj
    for(i=1;i<4*ca_krylov_s+1;i++){add_vectors(level, r_id,1.0, r_id,cj[i],PRrt[i]);}             //          ...
    }                                                                                           //
    m+=ca_krylov_s;                                                                           //   m+=ca_krylov_s;
  }                                                                                             // } // outer m loop
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #ifdef KRYLOV_DIAGONAL_PRECONDITION
  mul_vectors(level,e_id,1.0, VECTOR_DINV,e_id);                                                        //   e_id[] = Dinv[]*e_id[] // i.e. e = D^{-1}e'
  #endif
}
#endif // CA_KRYLOV_TELESCOPING
//------------------------------------------------------------------------------------------------------------------------------
