//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#ifndef M_PI
#define M_PI 3.14159265358979323846 // in case math.h doesn't define it
#endif
double evaluateBeta(double x, double y, double z, double h, int add_Bxx, int add_Byy, int add_Bzz){
  double b = 0.25;
  double a = 2.0*M_PI; // one period on [0,1]^3

  double B    = 1.0 + b*sin(a*x)*sin(a*y)*sin(a*z);
//double Bx   =     a*b*cos(a*x)*sin(a*y)*sin(a*z);
//double By   =     a*b*sin(a*x)*cos(a*y)*sin(a*z);
//double Bz   =     a*b*sin(a*x)*sin(a*y)*cos(a*z);
  double Bxx  =  -a*a*b*sin(a*x)*sin(a*y)*sin(a*z);
  double Byy  =  -a*a*b*sin(a*x)*sin(a*y)*sin(a*z);
  double Bzz   = -a*a*b*sin(a*x)*sin(a*y)*sin(a*z);

  // 4th order correction to approximate the conversion of cell-centered values to cell-averaged...
  if(add_Bxx)B+=(h*h/24.0)*Bxx;
  if(add_Byy)B+=(h*h/24.0)*Byy;
  if(add_Bzz)B+=(h*h/24.0)*Bzz;
  return(B);
}


//------------------------------------------------------------------------------------------------------------------------------
double evaluateF(double x, double y, double z, double h, int add_Fxx, int add_Fyy, int add_Fzz){
  #if 0 // harder problem... not sure I manually differentiated this right...
  // 8 'poles', one per octant
  double    cx = 0.75;
  double    cy = 0.75;
  double    cz = 0.75,sign = 1.0;
  if(x<0.5){cx = 0.25;sign*=-1.0;}
  if(y<0.5){cy = 0.25;sign*=-1.0;}
  if(z<0.5){cz = 0.25;sign*=-1.0;}

  double r0  = 0.1;
  double a   = M_PI/2/r0;
  double r   =  pow( (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz) ,  0.5); // euclidean distance
  double rx  =  pow( (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz) , -0.5)*(x-cx); // dr/dx
  double ry  =  pow( (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz) , -0.5)*(y-cy);
  double rz  =  pow( (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz) , -0.5)*(z-cz);
  double rxx = -pow( (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz) , -1.5)*(x-cx)*(x-cx) + pow( (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz) , -0.5); // d2r/dx2
  double ryy = -pow( (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz) , -1.5)*(y-cy)*(y-cy) + pow( (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz) , -0.5);
  double rzz = -pow( (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz) , -1.5)*(z-cz)*(z-cz) + pow( (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz) , -0.5);

  double p   = 6.0;
  double F   = sign*(        pow(cos(a*r),p  )    );
  double Fx  = sign*(   -a*p*pow(cos(a*r),p-1)*sin(a*r)*rx );
  double Fy  = sign*(   -a*p*pow(cos(a*r),p-1)*sin(a*r)*ry );
  double Fz  = sign*(   -a*p*pow(cos(a*r),p-1)*sin(a*r)*rz );
  double Fxx = sign*( -a*a*p*pow(cos(a*r),p  )*rx*rx  +  a*a*p*(p-1)*pow(cos(a*r),p-2)*pow(sin(a*r),2)*rx*rx  -  a*p*pow(cos(a*r),p-1)*sin(a*r)*rxx );
  double Fyy = sign*( -a*a*p*pow(cos(a*r),p  )*ry*ry  +  a*a*p*(p-1)*pow(cos(a*r),p-2)*pow(sin(a*r),2)*ry*ry  -  a*p*pow(cos(a*r),p-1)*sin(a*r)*ryy );
  double Fzz = sign*( -a*a*p*pow(cos(a*r),p  )*rz*rz  +  a*a*p*(p-1)*pow(cos(a*r),p-2)*pow(sin(a*r),2)*rz*rz  -  a*p*pow(cos(a*r),p-1)*sin(a*r)*rzz );

  if(r>=r0){
    F   = 0.0;
    Fx  = 0.0;
    Fy  = 0.0;
    Fz  = 0.0;
    Fxx = 0.0;
    Fyy = 0.0;
    Fzz = 0.0;
  }
  #else
  double a = 2.0*M_PI;
  double p = 7.0;
  double F   =        pow(sin(a*x),p  )*pow(sin(a*y),p  )*pow(sin(a*z),p  );
//double Fx  =    a*p*pow(sin(a*x),p-1)*pow(sin(a*y),p  )*pow(sin(a*z),p  )*cos(a*x);
//double Fy  =    a*p*pow(sin(a*x),p  )*pow(sin(a*y),p-1)*pow(sin(a*z),p  )*cos(a*y);
//double Fz  =    a*p*pow(sin(a*x),p  )*pow(sin(a*y),p  )*pow(sin(a*z),p-1)*cos(a*z);
  double Fxx = -a*a*p*pow(sin(a*x),p  )*pow(sin(a*y),p  )*pow(sin(a*z),p  )  +  a*a*p*(p-1)*pow(sin(a*x),p-2)*pow(sin(a*y),p  )*pow(sin(a*z),p  )*pow(cos(a*x),2);
  double Fyy = -a*a*p*pow(sin(a*x),p  )*pow(sin(a*y),p  )*pow(sin(a*z),p  )  +  a*a*p*(p-1)*pow(sin(a*x),p  )*pow(sin(a*y),p-2)*pow(sin(a*z),p  )*pow(cos(a*y),2);
  double Fzz = -a*a*p*pow(sin(a*x),p  )*pow(sin(a*y),p  )*pow(sin(a*z),p  )  +  a*a*p*(p-1)*pow(sin(a*x),p  )*pow(sin(a*y),p  )*pow(sin(a*z),p-2)*pow(cos(a*z),2);
  #endif

  // 4th order correction to approximate the conversion of cell-centered values to cell-averaged...
  if(add_Fxx)F+=(h*h/24.0)*Fxx;
  if(add_Fyy)F+=(h*h/24.0)*Fyy;
  if(add_Fzz)F+=(h*h/24.0)*Fzz;

  return(F);
}


//------------------------------------------------------------------------------------------------------------------------------
void initialize_problem(level_type * level, double hLevel, double a, double b){
  level->h = hLevel;

  int box;
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int   dim_i = level->my_boxes[box].dim;
    const int   dim_j = level->my_boxes[box].dim;
    const int   dim_k = level->my_boxes[box].dim;
    #ifdef _OPENMP
    #pragma omp parallel for private(k,j,i) collapse(3)
    #endif
    for(k=0;k<=dim_k;k++){ // include high face
    for(j=0;j<=dim_j;j++){ // include high face
    for(i=0;i<=dim_i;i++){ // include high face
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      int ijk = (i+ghosts) + (j+ghosts)*jStride + (k+ghosts)*kStride;
      double x = hLevel*( (double)(i+level->my_boxes[box].low.i) + 0.5 ); // +0.5 to get to the center of cell
      double y = hLevel*( (double)(j+level->my_boxes[box].low.j) + 0.5 );
      double z = hLevel*( (double)(k+level->my_boxes[box].low.k) + 0.5 );
      double A,Bi,Bj,Bk;
      //double A,B,Bx,By,Bz,Bi,Bj,Bk;
      //double U,Ux,Uy,Uz,Uxx,Uyy,Uzz;
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      A  = 1.0;
      Bi = 1.0;
      Bj = 1.0;
      Bk = 1.0;
      #ifdef STENCIL_VARIABLE_COEFFICIENT // variable coefficient problem...
      Bi=evaluateBeta(x-hLevel*0.5,y           ,z           ,hLevel,0,1,1); // face-centered value of Beta for beta_i
      Bj=evaluateBeta(x           ,y-hLevel*0.5,z           ,hLevel,1,0,1); // face-centered value of Beta for beta_j
      Bk=evaluateBeta(x           ,y           ,z-hLevel*0.5,hLevel,1,1,0); // face-centered value of Beta for beta_k
      #endif
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      double F=evaluateF(x,y,z,hLevel,1,1,1);
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      level->my_boxes[box].vectors[VECTOR_ALPHA ][ijk] = A;
      level->my_boxes[box].vectors[VECTOR_BETA_I][ijk] = Bi;
      level->my_boxes[box].vectors[VECTOR_BETA_J][ijk] = Bj;
      level->my_boxes[box].vectors[VECTOR_BETA_K][ijk] = Bk;
      level->my_boxes[box].vectors[VECTOR_F     ][ijk] = F;
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    }}}
  }

}
//------------------------------------------------------------------------------------------------------------------------------
