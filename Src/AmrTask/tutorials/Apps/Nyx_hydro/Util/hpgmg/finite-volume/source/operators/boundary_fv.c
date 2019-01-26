//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
void apply_BCs_v1(level_type * level, int x_id, int shape){
  // For cell-centered, we need to fill in the ghost zones to apply any BC's
  // This code does a simple linear interpolation for homogeneous dirichlet (0 on boundary)
  // Nominally, this is first performed across faces, then to edges, then to corners.  
  // In this implementation, these three steps are fused
  //
  //   . . . . . . . . . .        . . . . . . . . . .
  //   .       .       .          .       .       .
  //   .   ?   .   ?   .          .+x(0,0).-x(0,0).
  //   .       .       .          .       .       .
  //   . . . . +---0---+--        . . . . +-------+--
  //   .       |       |          .       |       |
  //   .   ?   0 x(0,0)|          .-x(0,0)| x(0,0)|
  //   .       |       |          .       |       |
  //   . . . . +-------+--        . . . . +-------+--
  //   .       |       |          .       |       |
  //
  //
  if(shape>=STENCIL_MAX_SHAPES)shape=STENCIL_SHAPE_BOX;  // shape must be < STENCIL_MAX_SHAPES in order to safely index into boundary_condition.blocks[]
  if(level->boundary_condition.type == BC_PERIODIC)return; // no BC's to apply !

  const int   faces[27] = {0,0,0,0,1,0,0,0,0,  0,1,0,1,0,1,0,1,0,  0,0,0,0,1,0,0,0,0};
  const int   edges[27] = {0,1,0,1,0,1,0,1,0,  1,0,1,0,0,0,1,0,1,  0,1,0,1,0,1,0,1,0};
  const int corners[27] = {1,0,1,0,0,0,1,0,1,  0,0,0,0,0,0,0,0,0,  1,0,1,0,0,0,1,0,1};

  int buffer;
  double _timeStart = getTime();
  PRAGMA_THREAD_ACROSS_BLOCKS(level,buffer,level->boundary_condition.num_blocks[shape])
  for(buffer=0;buffer<level->boundary_condition.num_blocks[shape];buffer++){
    double scale = 1.0;
    if(  faces[level->boundary_condition.blocks[shape][buffer].subtype])scale=-1.0;
    if(  edges[level->boundary_condition.blocks[shape][buffer].subtype])scale= 1.0;
    if(corners[level->boundary_condition.blocks[shape][buffer].subtype])scale=-1.0;

    int i,j,k;
    const int       box = level->boundary_condition.blocks[shape][buffer].read.box; 
    const int     dim_i = level->boundary_condition.blocks[shape][buffer].dim.i;
    const int     dim_j = level->boundary_condition.blocks[shape][buffer].dim.j;
    const int     dim_k = level->boundary_condition.blocks[shape][buffer].dim.k;
    const int       ilo = level->boundary_condition.blocks[shape][buffer].read.i;
    const int       jlo = level->boundary_condition.blocks[shape][buffer].read.j;
    const int       klo = level->boundary_condition.blocks[shape][buffer].read.k;
    const int normal = 26-level->boundary_condition.blocks[shape][buffer].subtype; // invert the normal vector
 
    // hard code for box to box BC's 
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    double * __restrict__  x = level->my_boxes[box].vectors[x_id] + level->box_ghosts*(1+jStride+kStride);

    // convert normal vector into pointer offsets...
    const int di = (((normal % 3)  )-1);
    const int dj = (((normal % 9)/3)-1);
    const int dk = (((normal / 9)  )-1);
    const int stride = di + dj*jStride + dk*kStride;

    if(dim_i==1){
      for(k=0;k<dim_k;k++){
      for(j=0;j<dim_j;j++){
        int ijk = (  ilo) + (j+jlo)*jStride + (k+klo)*kStride;
        x[ijk] = scale*x[ijk+stride]; // homogeneous linear = 1pt stencil
      }}
    }else if(dim_j==1){
      for(k=0;k<dim_k;k++){
      for(i=0;i<dim_i;i++){
        int ijk = (i+ilo) + (  jlo)*jStride + (k+klo)*kStride;
        x[ijk] = scale*x[ijk+stride]; // homogeneous linear = 1pt stencil
      }}
    }else if(dim_k==1){
      for(j=0;j<dim_j;j++){
      for(i=0;i<dim_i;i++){
        int ijk = (i+ilo) + (j+jlo)*jStride + (  klo)*kStride;
        x[ijk] = scale*x[ijk+stride]; // homogeneous linear = 1pt stencil
      }}
    }else{
      for(k=0;k<dim_k;k++){
      for(j=0;j<dim_j;j++){
      for(i=0;i<dim_i;i++){
        int ijk = (i+ilo) + (j+jlo)*jStride + (k+klo)*kStride;
        x[ijk] = scale*x[ijk+stride]; // homogeneous linear = 1pt stencil
      }}}
    }

  }
  level->timers.boundary_conditions += (double)(getTime()-_timeStart);
}

//------------------------------------------------------------------------------------------------------------------------------
// For cell-centered/averaged, one must fill in a ghost zone in order to affect a boundary condition
// The argument shape indicates on which regions of the domain (not the individual boxes) must the boundary condition be enforced.
//   If shape exceeds the range of defined shapes, the boundary condition will be applied to all faces, edges, and corners
// This code performs a simple quadratic volume averages extrapolation for homogeneous dirichlet (0 on boundary)
//   Nominally, this is first performed across faces, then to edges, then to corners.  
//   In this implementation, these three steps are fused
// This code will apply the BC only to the first ghost zone.  Subsequent (2nd, 3rd, ...) ghost zones will be zero'd
// This code will drop order if one attempts to apply quadratic BC's to boxes of less than 2^3
void apply_BCs_v2(level_type * level, int x_id, int shape){
  const int box_dim    = level->box_dim;
  const int box_ghosts = level->box_ghosts;
  if(shape>=STENCIL_MAX_SHAPES)shape=STENCIL_SHAPE_BOX;  // shape must be < STENCIL_MAX_SHAPES in order to safely index into boundary_condition.blocks[]
  if(level->boundary_condition.type == BC_PERIODIC)return; // no BC's to apply !
  if(level->box_dim<2){apply_BCs_v1(level,x_id,shape);return;}

  const int   faces[27] = {0,0,0,0,1,0,0,0,0,  0,1,0,1,0,1,0,1,0,  0,0,0,0,1,0,0,0,0};
  const int   edges[27] = {0,1,0,1,0,1,0,1,0,  1,0,1,0,0,0,1,0,1,  0,1,0,1,0,1,0,1,0};
  const int corners[27] = {1,0,1,0,0,0,1,0,1,  0,0,0,0,0,0,0,0,0,  1,0,1,0,0,0,1,0,1};

  int buffer;
  double _timeStart = getTime();
  PRAGMA_THREAD_ACROSS_BLOCKS(level,buffer,level->boundary_condition.num_blocks[shape])
  for(buffer=0;buffer<level->boundary_condition.num_blocks[shape];buffer++){
    int i,j,k;
    const int       box = level->boundary_condition.blocks[shape][buffer].read.box; 
    const int     dim_i = level->boundary_condition.blocks[shape][buffer].dim.i;
    const int     dim_j = level->boundary_condition.blocks[shape][buffer].dim.j;
    const int     dim_k = level->boundary_condition.blocks[shape][buffer].dim.k;
    const int       ilo = level->boundary_condition.blocks[shape][buffer].read.i;
    const int       jlo = level->boundary_condition.blocks[shape][buffer].read.j;
    const int       klo = level->boundary_condition.blocks[shape][buffer].read.k;
    const int   subtype = level->boundary_condition.blocks[shape][buffer].subtype;
  //const int    normal = 26-subtype;
 
    // hard code for box to box BC's 
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const double * __restrict__ x  = level->my_boxes[box].vectors[x_id] + level->box_ghosts*(1+jStride+kStride);
          double * __restrict__ xn = level->my_boxes[box].vectors[x_id] + level->box_ghosts*(1+jStride+kStride); // physically the same, but use different pointers for read/write

    // zero out entire ghost region when not all points will be updated...
    if(box_ghosts>1){
    for(k=0;k<dim_k;k++){
    for(j=0;j<dim_j;j++){
    for(i=0;i<dim_i;i++){
      int ijk = (i+ilo) + (j+jlo)*jStride + (k+klo)*kStride;
      xn[ijk] = 0.0;
    }}}}

    // apply the appropriate BC subtype (face, edge, corner)...
    if(faces[subtype]){
      //
      //    :.......|.......:.......:
      //    :   ?   | -5/2  :  1/2  :
      //    :.......|.......:.......:
      //
      int r=-1,rStride=-1,dim_r=-1,rlo=-1;
      int s=-1,sStride=-1,dim_s=-1,slo=-1;
      int t=-1,tStride=-1,dt=-1;
    
      // the two 4-point stencils can point in 6 different directions...
      switch(subtype){
        case  4:rlo=ilo;dim_r=dim_i;rStride=      1;slo=jlo;dim_s=dim_j;sStride=jStride;t=     -1;tStride=kStride;dt= tStride;break; // ij face, low k
        case 10:rlo=ilo;dim_r=dim_i;rStride=      1;slo=klo;dim_s=dim_k;sStride=kStride;t=     -1;tStride=jStride;dt= tStride;break; // ik face, low j
        case 12:rlo=jlo;dim_r=dim_j;rStride=jStride;slo=klo;dim_s=dim_k;sStride=kStride;t=     -1;tStride=      1;dt= tStride;break; // jk face, low i
        case 14:rlo=jlo;dim_r=dim_j;rStride=jStride;slo=klo;dim_s=dim_k;sStride=kStride;t=box_dim;tStride=      1;dt=-tStride;break; // jk face, high i
        case 16:rlo=ilo;dim_r=dim_i;rStride=      1;slo=klo;dim_s=dim_k;sStride=kStride;t=box_dim;tStride=jStride;dt=-tStride;break; // ik face, high j
        case 22:rlo=ilo;dim_r=dim_i;rStride=      1;slo=jlo;dim_s=dim_j;sStride=jStride;t=box_dim;tStride=kStride;dt=-tStride;break; // ij face, high k
      }
      // FIX... optimize for rStride==1 (unit-stride)
      for(s=0;s<dim_s;s++){
      for(r=0;r<dim_r;r++){
        int ijk = (r+rlo)*rStride + (s+slo)*sStride + (t)*tStride;
        xn[ijk] = -2.5*x[ijk+dt] + 0.5*x[ijk+2*dt];
      }}
    }else
    if(edges[subtype]){
      //
      //          r   +---+---+ dt
      //          ^  /   /   /|/
      //          | +---+---+ |
      //          |/   /   /|/
      //      +---+---+---+ |
      //     /   /|   |   |/
      //    +---+ |---+---+---> ds
      //    |   |/   
      //    +---+    
      //
      int r=-1,rStride=-1,dim_r=-1,rlo=-1;
      int s=-1,sStride=-1,ds=-1;
      int t=-1,tStride=-1,dt=-1;
      // the four 16-point stencils (symmetry allows you to view it as 12 4-point) can point in 12 different directions...
      switch(subtype){
        case  1:rlo=ilo;dim_r=dim_i;rStride=      1;s=     -1;sStride=jStride;t=     -1;tStride=kStride;ds= sStride;dt= tStride;break; // i-edge,  low j,  low k
        case  3:rlo=jlo;dim_r=dim_j;rStride=jStride;s=     -1;sStride=      1;t=     -1;tStride=kStride;ds= sStride;dt= tStride;break; // j-edge,  low i,  low k
        case  5:rlo=jlo;dim_r=dim_j;rStride=jStride;s=box_dim;sStride=      1;t=     -1;tStride=kStride;ds=-sStride;dt= tStride;break; // j-edge, high i,  low k
        case  7:rlo=ilo;dim_r=dim_i;rStride=      1;s=box_dim;sStride=jStride;t=     -1;tStride=kStride;ds=-sStride;dt= tStride;break; // i-edge, high j,  low k
        case  9:rlo=klo;dim_r=dim_k;rStride=kStride;s=     -1;sStride=      1;t=     -1;tStride=jStride;ds= sStride;dt= tStride;break; // k-edge,  low i,  low j
        case 11:rlo=klo;dim_r=dim_k;rStride=kStride;s=box_dim;sStride=      1;t=     -1;tStride=jStride;ds=-sStride;dt= tStride;break; // k-edge, high i,  low j
        case 15:rlo=klo;dim_r=dim_k;rStride=kStride;s=     -1;sStride=      1;t=box_dim;tStride=jStride;ds= sStride;dt=-tStride;break; // k-edge,  low i, high j
        case 17:rlo=klo;dim_r=dim_k;rStride=kStride;s=box_dim;sStride=      1;t=box_dim;tStride=jStride;ds=-sStride;dt=-tStride;break; // k-edge, high i, high j
        case 19:rlo=ilo;dim_r=dim_i;rStride=      1;s=     -1;sStride=jStride;t=box_dim;tStride=kStride;ds= sStride;dt=-tStride;break; // i-edge,  low j, high k
        case 21:rlo=jlo;dim_r=dim_j;rStride=jStride;s=     -1;sStride=      1;t=box_dim;tStride=kStride;ds= sStride;dt=-tStride;break; // j-edge,  low i, high k
        case 23:rlo=jlo;dim_r=dim_j;rStride=jStride;s=box_dim;sStride=      1;t=box_dim;tStride=kStride;ds=-sStride;dt=-tStride;break; // j-edge, high i, high k
        case 25:rlo=ilo;dim_r=dim_i;rStride=      1;s=box_dim;sStride=jStride;t=box_dim;tStride=kStride;ds=-sStride;dt=-tStride;break; // i-edge, high j, high k
      }
      // FIX... optimize for rStride==1 (unit-stride)
      for(r=0;r<dim_r;r++){
        int ijk = (r+rlo)*rStride + (s)*sStride + (t)*tStride;
        xn[ijk] =   6.25*x[ijk+  ds+  dt] 
                  - 1.25*x[ijk+2*ds+  dt]
                  - 1.25*x[ijk+  ds+2*dt]
                  + 0.25*x[ijk+2*ds+2*dt];
      }
    }else
    if(corners[subtype]){
      //
      //                  +---+---+
      //                 /   /   /|
      //                +---+---+ |
      //               /   /   /|/|
      //              +---+---+ | | 
      //              |   |   |/|/
      //              +---+---+ |
      //              |   |   |/
      //          +---+---+---+
      //         /   /|
      //        +---+ |
      //        |   |/
      //        +---+
      //
      int i=-1,di=-1;
      int j=-1,dj=-1;
      int k=-1,dk=-1;
      // the eight 64-point stencils (symmetry allows you to view it as 56 4-point) can point in 8 different directions...
      switch(subtype){
        case  0:i=     -1;j=     -1;k=     -1;di= 1;dj= jStride;dk= kStride;break; //  low i,  low j,  low k
        case  2:i=box_dim;j=     -1;k=     -1;di=-1;dj= jStride;dk= kStride;break; // high i,  low j,  low k
        case  6:i=     -1;j=box_dim;k=     -1;di= 1;dj=-jStride;dk= kStride;break; //  low i, high j,  low k
        case  8:i=box_dim;j=box_dim;k=     -1;di=-1;dj=-jStride;dk= kStride;break; // high i, high j,  low k
        case 18:i=     -1;j=     -1;k=box_dim;di= 1;dj= jStride;dk=-kStride;break; //  low i,  low j, high k
        case 20:i=box_dim;j=     -1;k=box_dim;di=-1;dj= jStride;dk=-kStride;break; // high i,  low j, high k
        case 24:i=     -1;j=box_dim;k=box_dim;di= 1;dj=-jStride;dk=-kStride;break; //  low i, high j, high k
        case 26:i=box_dim;j=box_dim;k=box_dim;di=-1;dj=-jStride;dk=-kStride;break; // high i, high j, high k
      }
      int ijk = (i) + (j)*jStride + (k)*kStride;
      xn[ijk] =  -15.625*x[ijk+  di+  dj+  dk] 
                 + 3.125*x[ijk+2*di+  dj+  dk] 
                 + 3.125*x[ijk+  di+2*dj+  dk] 
                 + 3.125*x[ijk+  di+  dj+2*dk] 
                 - 0.625*x[ijk+2*di+2*dj+  dk] 
                 - 0.625*x[ijk+  di+2*dj+2*dk] 
                 - 0.625*x[ijk+2*di+  dj+2*dk] 
                 + 0.125*x[ijk+2*di+2*dj+2*dk];
    }
  }
  level->timers.boundary_conditions += (double)(getTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
// For cell-centered/averaged, one must fill in a ghost zone in order to affect a boundary condition
// The argument shape indicates on which regions of the domain (not the individual boxes) must the boundary condition be enforced.
//   If shape exceeds the range of defined shapes, the boundary condition will be applied to all faces, edges, and corners
// This code performs a simple quartic volume averages extrapolation for homogeneous dirichlet (0 on boundary)
//   Nominally, this is first performed across faces, then to edges, then to corners.  
//   In this implementation, these three steps are fused
// It is considered an error to call this routine if the domain has less that two ghost zones
// This code will drop order if one attempts to apply quartic BC's to boxes of less than 4^3
void apply_BCs_v4(level_type * level, int x_id, int shape){
  const int box_dim    = level->box_dim;
  const int box_ghosts = level->box_ghosts;
  if(shape>=STENCIL_MAX_SHAPES)shape=STENCIL_SHAPE_BOX;  // shape must be < STENCIL_MAX_SHAPES in order to safely index into boundary_condition.blocks[]
  if(level->boundary_condition.type == BC_PERIODIC)return; // no BC's to apply !
  if(box_ghosts<2){fprintf(stderr,"called quartic BC's with only 1 ghost zone!!!\n");exit(0);}
//if(box_dim   <4){fprintf(stderr,"called quartic BC's with boxes < 4^3 \n");exit(0);}
  if(box_dim   <4){apply_BCs_v2(level,x_id,shape);return;} // FIX... is it safe to drop order on the boundary on coarse grids ??

  const int   faces[27] = {0,0,0,0,1,0,0,0,0,  0,1,0,1,0,1,0,1,0,  0,0,0,0,1,0,0,0,0};
  const int   edges[27] = {0,1,0,1,0,1,0,1,0,  1,0,1,0,0,0,1,0,1,  0,1,0,1,0,1,0,1,0};
  const int corners[27] = {1,0,1,0,0,0,1,0,1,  0,0,0,0,0,0,0,0,0,  1,0,1,0,0,0,1,0,1};

  int buffer;
  double _timeStart = getTime();
  PRAGMA_THREAD_ACROSS_BLOCKS(level,buffer,level->boundary_condition.num_blocks[shape])
  for(buffer=0;buffer<level->boundary_condition.num_blocks[shape];buffer++){
    int i,j,k;
    const int       box = level->boundary_condition.blocks[shape][buffer].read.box; 
    const int     dim_i = level->boundary_condition.blocks[shape][buffer].dim.i;
    const int     dim_j = level->boundary_condition.blocks[shape][buffer].dim.j;
    const int     dim_k = level->boundary_condition.blocks[shape][buffer].dim.k;
    const int       ilo = level->boundary_condition.blocks[shape][buffer].read.i;
    const int       jlo = level->boundary_condition.blocks[shape][buffer].read.j;
    const int       klo = level->boundary_condition.blocks[shape][buffer].read.k;
    const int   subtype = level->boundary_condition.blocks[shape][buffer].subtype;
  //const int    normal = 26-subtype;
 
    // hard code for box to box BC's 
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const double * __restrict__ x  = level->my_boxes[box].vectors[x_id] + level->box_ghosts*(1+jStride+kStride);
          double * __restrict__ xn = level->my_boxes[box].vectors[x_id] + level->box_ghosts*(1+jStride+kStride); // physically the same, but use different pointers for read/write

    double OneTwelfth = 1.0/12.0;

    // zero out entire ghost region when not all points will be updated...
    if(box_ghosts>2){
    for(k=0;k<dim_k;k++){
    for(j=0;j<dim_j;j++){
    for(i=0;i<dim_i;i++){
      int ijk = (i+ilo) + (j+jlo)*jStride + (k+klo)*kStride;
      xn[ijk] = 0.0;
    }}}}

    // apply the appropriate BC subtype (face, edge, corner)...
    if(faces[subtype]){
      //
      //    :....:....|....:....:....:....:.
      //    : ?? : ?? | x1 : x2 : x3 : x4 :
      //    :....:....|....:....:....:....:.
      //
      int r=-1,rStride=-1,dim_r=-1,rlo=-1;
      int s=-1,sStride=-1,dim_s=-1,slo=-1;
      int t=-1,tStride=-1,dt=-1;
    
      // the two 4-point stencils can point in 6 different directions...
      switch(subtype){
        case  4:rlo=ilo;dim_r=dim_i;rStride=      1;slo=jlo;dim_s=dim_j;sStride=jStride;t=     -1;tStride=kStride;dt= tStride;break; // ij face, low k
        case 10:rlo=ilo;dim_r=dim_i;rStride=      1;slo=klo;dim_s=dim_k;sStride=kStride;t=     -1;tStride=jStride;dt= tStride;break; // ik face, low j
        case 12:rlo=jlo;dim_r=dim_j;rStride=jStride;slo=klo;dim_s=dim_k;sStride=kStride;t=     -1;tStride=      1;dt= tStride;break; // jk face, low i
        case 14:rlo=jlo;dim_r=dim_j;rStride=jStride;slo=klo;dim_s=dim_k;sStride=kStride;t=box_dim;tStride=      1;dt=-tStride;break; // jk face, high i
        case 16:rlo=ilo;dim_r=dim_i;rStride=      1;slo=klo;dim_s=dim_k;sStride=kStride;t=box_dim;tStride=jStride;dt=-tStride;break; // ik face, high j
        case 22:rlo=ilo;dim_r=dim_i;rStride=      1;slo=jlo;dim_s=dim_j;sStride=jStride;t=box_dim;tStride=kStride;dt=-tStride;break; // ij face, high k
      }
      // FIX... optimize for rStride==1 (unit-stride)
      // FIX... optimize for dt==+/-1
      double * __restrict__  ghost0 = (double * __restrict__)(x   ); // convince the compiler that read (box) & write (ghost zone) are disjoint
      double * __restrict__  ghost1 = (double * __restrict__)(x-dt); // convince the compiler that read (box) & write (ghost zone) are disjoint
      for(s=0;s<dim_s;s++){
      for(r=0;r<dim_r;r++){
        int ijk = (r+rlo)*rStride + (s+slo)*sStride + (t)*tStride;
        double x1=x[ijk+dt], x2=x[ijk+2*dt], x3=x[ijk+3*dt], x4=x[ijk+4*dt];
        #if 0 // compiler cannot disambiguate xn[ijk] and xn[ijk-dt]
        xn[ijk   ] = OneTwelfth*(  -77.0*x1 +  43.0*x2 -  17.0*x3 +  3.0*x4 );
        xn[ijk-dt] = OneTwelfth*( -505.0*x1 + 335.0*x2 - 145.0*x3 + 27.0*x4 );
        #else
        ghost0[ijk] = OneTwelfth*(  -77.0*x1 +  43.0*x2 -  17.0*x3 +  3.0*x4 );
        ghost1[ijk] = OneTwelfth*( -505.0*x1 + 335.0*x2 - 145.0*x3 + 27.0*x4 );
        #endif
      }}
    }else
    if(edges[subtype]){
      //
      //                        +---+---+---+---+ dt
      //                       /   /   /   /   /|/
      //                      +---+---+---+---+ |
      //                r    /   /   /   /   /|/
      //                ^   +---+---+---+---+ |
      //                |  /   /   /   /   /|/
      //                | +---+---+---+---+ |
      //                |/   /   /   /   /|/
      //        +---+---+---+---+---+---+ |
      //       /   /   /|   |   |   |   |/
      //      +---+---+ |---+---+---+---+---> ds
      //     /   /   /|/ 
      //    +---+---+ |  
      //    |   |   |/   
      //    +---+---+    
      //
      //              ^ dt
      //              |
      //    :....:....|....:....:....:....:.
      //    : f4 : n4 | 14 : 24 : 34 : 44 :
      //    :....:....|....:....:....:....:.
      //    : f3 : n3 | 13 : 23 : 33 : 43 :
      //    :....:....|....:....:....:....:.
      //    : f2 : n2 | 12 : 22 : 32 : 42 :
      //    :....:....|....:....:....:....:.
      //    : f1 : n1 | 11 : 21 : 31 : 41 :
      //    ----------+---------------------> ds
      //    : ?? : ?? |
      //    :....:....|
      //    : ?? : ?? |
      //    :....:....|
      //
      int r=-1,rStride=-1,dim_r=-1,rlo=-1;
      int s=-1,sStride=-1,ds=-1;
      int t=-1,tStride=-1,dt=-1;
      // the four 16-point stencils (symmetry allows you to view it as 12 4-point) can point in 12 different directions...
      switch(subtype){
        case  1:rlo=ilo;dim_r=dim_i;rStride=      1;s=     -1;sStride=jStride;t=     -1;tStride=kStride;ds= sStride;dt= tStride;break; // i-edge,  low j,  low k
        case  3:rlo=jlo;dim_r=dim_j;rStride=jStride;s=     -1;sStride=      1;t=     -1;tStride=kStride;ds= sStride;dt= tStride;break; // j-edge,  low i,  low k
        case  5:rlo=jlo;dim_r=dim_j;rStride=jStride;s=box_dim;sStride=      1;t=     -1;tStride=kStride;ds=-sStride;dt= tStride;break; // j-edge, high i,  low k
        case  7:rlo=ilo;dim_r=dim_i;rStride=      1;s=box_dim;sStride=jStride;t=     -1;tStride=kStride;ds=-sStride;dt= tStride;break; // i-edge, high j,  low k
        case  9:rlo=klo;dim_r=dim_k;rStride=kStride;s=     -1;sStride=      1;t=     -1;tStride=jStride;ds= sStride;dt= tStride;break; // k-edge,  low i,  low j
        case 11:rlo=klo;dim_r=dim_k;rStride=kStride;s=box_dim;sStride=      1;t=     -1;tStride=jStride;ds=-sStride;dt= tStride;break; // k-edge, high i,  low j
        case 15:rlo=klo;dim_r=dim_k;rStride=kStride;s=     -1;sStride=      1;t=box_dim;tStride=jStride;ds= sStride;dt=-tStride;break; // k-edge,  low i, high j
        case 17:rlo=klo;dim_r=dim_k;rStride=kStride;s=box_dim;sStride=      1;t=box_dim;tStride=jStride;ds=-sStride;dt=-tStride;break; // k-edge, high i, high j
        case 19:rlo=ilo;dim_r=dim_i;rStride=      1;s=     -1;sStride=jStride;t=box_dim;tStride=kStride;ds= sStride;dt=-tStride;break; // i-edge,  low j, high k
        case 21:rlo=jlo;dim_r=dim_j;rStride=jStride;s=     -1;sStride=      1;t=box_dim;tStride=kStride;ds= sStride;dt=-tStride;break; // j-edge,  low i, high k
        case 23:rlo=jlo;dim_r=dim_j;rStride=jStride;s=box_dim;sStride=      1;t=box_dim;tStride=kStride;ds=-sStride;dt=-tStride;break; // j-edge, high i, high k
        case 25:rlo=ilo;dim_r=dim_i;rStride=      1;s=box_dim;sStride=jStride;t=box_dim;tStride=kStride;ds=-sStride;dt=-tStride;break; // i-edge, high j, high k
      }
      // FIX... optimize for rStride==1 (unit-stride)
      // FIX... optimize for ds==+/-1
      double * __restrict__  ghost00 = (double * __restrict__)(x      ); // convince the compiler that read (box) & write (ghost zone) are disjoint
      double * __restrict__  ghost01 = (double * __restrict__)(x   -dt); // convince the compiler that read (box) & write (ghost zone) are disjoint
      double * __restrict__  ghost10 = (double * __restrict__)(x-ds   ); // convince the compiler that read (box) & write (ghost zone) are disjoint
      double * __restrict__  ghost11 = (double * __restrict__)(x-ds-dt); // convince the compiler that read (box) & write (ghost zone) are disjoint
      for(r=0;r<dim_r;r++){
        int ijk = (r+rlo)*rStride + (s)*sStride + (t)*tStride;
        double x11 = x[ijk+  ds+  dt], x21 = x[ijk+2*ds+  dt], x31 = x[ijk+3*ds+  dt], x41 = x[ijk+4*ds+  dt];
        double x12 = x[ijk+  ds+2*dt], x22 = x[ijk+2*ds+2*dt], x32 = x[ijk+3*ds+2*dt], x42 = x[ijk+4*ds+2*dt];
        double x13 = x[ijk+  ds+3*dt], x23 = x[ijk+2*ds+3*dt], x33 = x[ijk+3*ds+3*dt], x43 = x[ijk+4*ds+3*dt];
        double x14 = x[ijk+  ds+4*dt], x24 = x[ijk+2*ds+4*dt], x34 = x[ijk+3*ds+4*dt], x44 = x[ijk+4*ds+4*dt];
            double n1 = OneTwelfth*(  -77.0*x11 +  43.0*x21 -  17.0*x31 +  3.0*x41 );
            double n2 = OneTwelfth*(  -77.0*x12 +  43.0*x22 -  17.0*x32 +  3.0*x42 );
            double n3 = OneTwelfth*(  -77.0*x13 +  43.0*x23 -  17.0*x33 +  3.0*x43 );
            double n4 = OneTwelfth*(  -77.0*x14 +  43.0*x24 -  17.0*x34 +  3.0*x44 );
            double f1 = OneTwelfth*( -505.0*x11 + 335.0*x21 - 145.0*x31 + 27.0*x41 );
            double f2 = OneTwelfth*( -505.0*x12 + 335.0*x22 - 145.0*x32 + 27.0*x42 );
            double f3 = OneTwelfth*( -505.0*x13 + 335.0*x23 - 145.0*x33 + 27.0*x43 );
            double f4 = OneTwelfth*( -505.0*x14 + 335.0*x24 - 145.0*x34 + 27.0*x44 );
        #if 0 // compiler cannot disambiguate pointers
        xn[ijk      ] = OneTwelfth*(  -77.0*n1  +  43.0*n2  -  17.0*n3  +  3.0*n4  );
        xn[ijk   -dt] = OneTwelfth*( -505.0*n1  + 335.0*n2  - 145.0*n3  + 27.0*n4  );
        xn[ijk-ds   ] = OneTwelfth*(  -77.0*f1  +  43.0*f2  -  17.0*f3  +  3.0*f4  );
        xn[ijk-ds-dt] = OneTwelfth*( -505.0*f1  + 335.0*f2  - 145.0*f3  + 27.0*f4  );
        #else
        ghost00[ijk] = OneTwelfth*(  -77.0*n1  +  43.0*n2  -  17.0*n3  +  3.0*n4  );
        ghost01[ijk] = OneTwelfth*( -505.0*n1  + 335.0*n2  - 145.0*n3  + 27.0*n4  );
        ghost10[ijk] = OneTwelfth*(  -77.0*f1  +  43.0*f2  -  17.0*f3  +  3.0*f4  );
        ghost11[ijk] = OneTwelfth*( -505.0*f1  + 335.0*f2  - 145.0*f3  + 27.0*f4  );
        #endif
      }
    }else
    if(corners[subtype]){
      //
      //                        +---+---+---+---+
      //                       /   /   /   /   /|
      //                      +---+---+---+---+ |
      //                     /   /   /   /   /|/|
      //                    +---+---+---+---+ | |
      //                   /   /   /   /   /|/|/|
      //                  +---+---+---+---+ | | |
      //                 /   /   /   /   /|/|/|/|
      //                +---+---+---+---+ | | | |
      //                |   |   |   |   |/|/|/|/
      //                +---+---+---+---+ | | |
      //                |   |   |   |   |/|/|/
      //                +---+---+---+---+ | |
      //                |   |   |   |   |/|/
      //                +---+---+---+---+ |
      //                |   |   |   |   |/
      //        +---+---+---+---+---+---+
      //       /   /   /|
      //      +---+---+ |
      //     /   /   /|/|
      //    +---+---+ | |
      //    |   |   |/|/
      //    +---+---+ |
      //    |   |   |/
      //    +---+---+
      //
      int i=-1,di=-1;
      int j=-1,dj=-1;
      int k=-1,dk=-1;
      // the eight 64-point stencils (symmetry allows you to view it as 56 4-point) can point in 8 different directions...
      // FIX... optimize for di==+/-1
      switch(subtype){
        case  0:i=     -1;j=     -1;k=     -1;di= 1;dj= jStride;dk= kStride;break; //  low i,  low j,  low k
        case  2:i=box_dim;j=     -1;k=     -1;di=-1;dj= jStride;dk= kStride;break; // high i,  low j,  low k
        case  6:i=     -1;j=box_dim;k=     -1;di= 1;dj=-jStride;dk= kStride;break; //  low i, high j,  low k
        case  8:i=box_dim;j=box_dim;k=     -1;di=-1;dj=-jStride;dk= kStride;break; // high i, high j,  low k
        case 18:i=     -1;j=     -1;k=box_dim;di= 1;dj= jStride;dk=-kStride;break; //  low i,  low j, high k
        case 20:i=box_dim;j=     -1;k=box_dim;di=-1;dj= jStride;dk=-kStride;break; // high i,  low j, high k
        case 24:i=     -1;j=box_dim;k=box_dim;di= 1;dj=-jStride;dk=-kStride;break; //  low i, high j, high k
        case 26:i=box_dim;j=box_dim;k=box_dim;di=-1;dj=-jStride;dk=-kStride;break; // high i, high j, high k
      }
      int ijk = (i) + (j)*jStride + (k)*kStride;
      double x144 = x[ijk+  di+4*dj+4*dk];double x244 = x[ijk+2*di+4*dj+4*dk];double x344 = x[ijk+3*di+4*dj+4*dk];double x444 = x[ijk+4*di+4*dj+4*dk];
      double x134 = x[ijk+  di+3*dj+4*dk];double x234 = x[ijk+2*di+3*dj+4*dk];double x334 = x[ijk+3*di+3*dj+4*dk];double x434 = x[ijk+4*di+3*dj+4*dk];
      double x124 = x[ijk+  di+2*dj+4*dk];double x224 = x[ijk+2*di+2*dj+4*dk];double x324 = x[ijk+3*di+2*dj+4*dk];double x424 = x[ijk+4*di+2*dj+4*dk];
      double x114 = x[ijk+  di+  dj+4*dk];double x214 = x[ijk+2*di+  dj+4*dk];double x314 = x[ijk+3*di+  dj+4*dk];double x414 = x[ijk+4*di+  dj+4*dk];

      double x143 = x[ijk+  di+4*dj+3*dk];double x243 = x[ijk+2*di+4*dj+3*dk];double x343 = x[ijk+3*di+4*dj+3*dk];double x443 = x[ijk+4*di+4*dj+3*dk];
      double x133 = x[ijk+  di+3*dj+3*dk];double x233 = x[ijk+2*di+3*dj+3*dk];double x333 = x[ijk+3*di+3*dj+3*dk];double x433 = x[ijk+4*di+3*dj+3*dk];
      double x123 = x[ijk+  di+2*dj+3*dk];double x223 = x[ijk+2*di+2*dj+3*dk];double x323 = x[ijk+3*di+2*dj+3*dk];double x423 = x[ijk+4*di+2*dj+3*dk];
      double x113 = x[ijk+  di+  dj+3*dk];double x213 = x[ijk+2*di+  dj+3*dk];double x313 = x[ijk+3*di+  dj+3*dk];double x413 = x[ijk+4*di+  dj+3*dk];

      double x142 = x[ijk+  di+4*dj+2*dk];double x242 = x[ijk+2*di+4*dj+2*dk];double x342 = x[ijk+3*di+4*dj+2*dk];double x442 = x[ijk+4*di+4*dj+2*dk];
      double x132 = x[ijk+  di+3*dj+2*dk];double x232 = x[ijk+2*di+3*dj+2*dk];double x332 = x[ijk+3*di+3*dj+2*dk];double x432 = x[ijk+4*di+3*dj+2*dk];
      double x122 = x[ijk+  di+2*dj+2*dk];double x222 = x[ijk+2*di+2*dj+2*dk];double x322 = x[ijk+3*di+2*dj+2*dk];double x422 = x[ijk+4*di+2*dj+2*dk];
      double x112 = x[ijk+  di+  dj+2*dk];double x212 = x[ijk+2*di+  dj+2*dk];double x312 = x[ijk+3*di+  dj+2*dk];double x412 = x[ijk+4*di+  dj+2*dk];

      double x141 = x[ijk+  di+4*dj+  dk];double x241 = x[ijk+2*di+4*dj+  dk];double x341 = x[ijk+3*di+4*dj+  dk];double x441 = x[ijk+4*di+4*dj+  dk];
      double x131 = x[ijk+  di+3*dj+  dk];double x231 = x[ijk+2*di+3*dj+  dk];double x331 = x[ijk+3*di+3*dj+  dk];double x431 = x[ijk+4*di+3*dj+  dk];
      double x121 = x[ijk+  di+2*dj+  dk];double x221 = x[ijk+2*di+2*dj+  dk];double x321 = x[ijk+3*di+2*dj+  dk];double x421 = x[ijk+4*di+2*dj+  dk];
      double x111 = x[ijk+  di+  dj+  dk];double x211 = x[ijk+2*di+  dj+  dk];double x311 = x[ijk+3*di+  dj+  dk];double x411 = x[ijk+4*di+  dj+  dk];

      // 32 stencils in i...
      double n11 = OneTwelfth*(  -77.0*x111 +  43.0*x211 -  17.0*x311 +  3.0*x411 );
      double n21 = OneTwelfth*(  -77.0*x121 +  43.0*x221 -  17.0*x321 +  3.0*x421 );
      double n31 = OneTwelfth*(  -77.0*x131 +  43.0*x231 -  17.0*x331 +  3.0*x431 );
      double n41 = OneTwelfth*(  -77.0*x141 +  43.0*x241 -  17.0*x341 +  3.0*x441 );
      double n12 = OneTwelfth*(  -77.0*x112 +  43.0*x212 -  17.0*x312 +  3.0*x412 );
      double n22 = OneTwelfth*(  -77.0*x122 +  43.0*x222 -  17.0*x322 +  3.0*x422 );
      double n32 = OneTwelfth*(  -77.0*x132 +  43.0*x232 -  17.0*x332 +  3.0*x432 );
      double n42 = OneTwelfth*(  -77.0*x142 +  43.0*x242 -  17.0*x342 +  3.0*x442 );
      double n13 = OneTwelfth*(  -77.0*x113 +  43.0*x213 -  17.0*x313 +  3.0*x413 );
      double n23 = OneTwelfth*(  -77.0*x123 +  43.0*x223 -  17.0*x323 +  3.0*x423 );
      double n33 = OneTwelfth*(  -77.0*x133 +  43.0*x233 -  17.0*x333 +  3.0*x433 );
      double n43 = OneTwelfth*(  -77.0*x143 +  43.0*x243 -  17.0*x343 +  3.0*x443 );
      double n14 = OneTwelfth*(  -77.0*x114 +  43.0*x214 -  17.0*x314 +  3.0*x414 );
      double n24 = OneTwelfth*(  -77.0*x124 +  43.0*x224 -  17.0*x324 +  3.0*x424 );
      double n34 = OneTwelfth*(  -77.0*x134 +  43.0*x234 -  17.0*x334 +  3.0*x434 );
      double n44 = OneTwelfth*(  -77.0*x144 +  43.0*x244 -  17.0*x344 +  3.0*x444 );

      double f11 = OneTwelfth*( -505.0*x111 + 335.0*x211 - 145.0*x311 +  27.0*x411 );
      double f21 = OneTwelfth*( -505.0*x121 + 335.0*x221 - 145.0*x321 +  27.0*x421 );
      double f31 = OneTwelfth*( -505.0*x131 + 335.0*x231 - 145.0*x331 +  27.0*x431 );
      double f41 = OneTwelfth*( -505.0*x141 + 335.0*x241 - 145.0*x341 +  27.0*x441 );
      double f12 = OneTwelfth*( -505.0*x112 + 335.0*x212 - 145.0*x312 +  27.0*x412 );
      double f22 = OneTwelfth*( -505.0*x122 + 335.0*x222 - 145.0*x322 +  27.0*x422 );
      double f32 = OneTwelfth*( -505.0*x132 + 335.0*x232 - 145.0*x332 +  27.0*x432 );
      double f42 = OneTwelfth*( -505.0*x142 + 335.0*x242 - 145.0*x342 +  27.0*x442 );
      double f13 = OneTwelfth*( -505.0*x113 + 335.0*x213 - 145.0*x313 +  27.0*x413 );
      double f23 = OneTwelfth*( -505.0*x123 + 335.0*x223 - 145.0*x323 +  27.0*x423 );
      double f33 = OneTwelfth*( -505.0*x133 + 335.0*x233 - 145.0*x333 +  27.0*x433 );
      double f43 = OneTwelfth*( -505.0*x143 + 335.0*x243 - 145.0*x343 +  27.0*x443 );
      double f14 = OneTwelfth*( -505.0*x114 + 335.0*x214 - 145.0*x314 +  27.0*x414 );
      double f24 = OneTwelfth*( -505.0*x124 + 335.0*x224 - 145.0*x324 +  27.0*x424 );
      double f34 = OneTwelfth*( -505.0*x134 + 335.0*x234 - 145.0*x334 +  27.0*x434 );
      double f44 = OneTwelfth*( -505.0*x144 + 335.0*x244 - 145.0*x344 +  27.0*x444 );

      // 16 stencils in j...
      double nn1 = OneTwelfth*(  -77.0*n11 +  43.0*n21 -  17.0*n31 +  3.0*n41 );
      double nn2 = OneTwelfth*(  -77.0*n12 +  43.0*n22 -  17.0*n32 +  3.0*n42 );
      double nn3 = OneTwelfth*(  -77.0*n13 +  43.0*n23 -  17.0*n33 +  3.0*n43 );
      double nn4 = OneTwelfth*(  -77.0*n14 +  43.0*n24 -  17.0*n34 +  3.0*n44 );
      double nf1 = OneTwelfth*( -505.0*n11 + 335.0*n21 - 145.0*n31 + 27.0*n41 );
      double nf2 = OneTwelfth*( -505.0*n12 + 335.0*n22 - 145.0*n32 + 27.0*n42 );
      double nf3 = OneTwelfth*( -505.0*n13 + 335.0*n23 - 145.0*n33 + 27.0*n43 );
      double nf4 = OneTwelfth*( -505.0*n14 + 335.0*n24 - 145.0*n34 + 27.0*n44 );

      double fn1 = OneTwelfth*(  -77.0*f11 +  43.0*f21 -  17.0*f31 +  3.0*f41 );
      double fn2 = OneTwelfth*(  -77.0*f12 +  43.0*f22 -  17.0*f32 +  3.0*f42 );
      double fn3 = OneTwelfth*(  -77.0*f13 +  43.0*f23 -  17.0*f33 +  3.0*f43 );
      double fn4 = OneTwelfth*(  -77.0*f14 +  43.0*f24 -  17.0*f34 +  3.0*f44 );
      double ff1 = OneTwelfth*( -505.0*f11 + 335.0*f21 - 145.0*f31 + 27.0*f41 );
      double ff2 = OneTwelfth*( -505.0*f12 + 335.0*f22 - 145.0*f32 + 27.0*f42 );
      double ff3 = OneTwelfth*( -505.0*f13 + 335.0*f23 - 145.0*f33 + 27.0*f43 );
      double ff4 = OneTwelfth*( -505.0*f14 + 335.0*f24 - 145.0*f34 + 27.0*f44 );

      //  8 stencils in k...
      double nnn = OneTwelfth*(  -77.0*nn1 +  43.0*nn2 -  17.0*nn3 +  3.0*nn4 );
      double nnf = OneTwelfth*( -505.0*nn1 + 335.0*nn2 - 145.0*nn3 + 27.0*nn4 );
      double nfn = OneTwelfth*(  -77.0*nf1 +  43.0*nf2 -  17.0*nf3 +  3.0*nf4 );
      double nff = OneTwelfth*( -505.0*nf1 + 335.0*nf2 - 145.0*nf3 + 27.0*nf4 );
      double fnn = OneTwelfth*(  -77.0*fn1 +  43.0*fn2 -  17.0*fn3 +  3.0*fn4 );
      double fnf = OneTwelfth*( -505.0*fn1 + 335.0*fn2 - 145.0*fn3 + 27.0*fn4 );
      double ffn = OneTwelfth*(  -77.0*ff1 +  43.0*ff2 -  17.0*ff3 +  3.0*ff4 );
      double fff = OneTwelfth*( -505.0*ff1 + 335.0*ff2 - 145.0*ff3 + 27.0*ff4 );

      // commit to the 8 ghost zones in this corner...
      xn[ijk         ] = nnn;
      xn[ijk      -dk] = nnf;
      xn[ijk   -dj   ] = nfn;
      xn[ijk   -dj-dk] = nff;
      xn[ijk-di      ] = fnn;
      xn[ijk-di   -dk] = fnf;
      xn[ijk-di-dj   ] = ffn;
      xn[ijk-di-dj-dk] = fff;
    }
  }
  level->timers.boundary_conditions += (double)(getTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
void extrapolate_betas(level_type * level){
  if(level->boundary_condition.type == BC_PERIODIC)return; // no BC's to apply !
  int shape=0;

  int buffer;
  double _timeStart = getTime();
  PRAGMA_THREAD_ACROSS_BLOCKS(level,buffer,level->boundary_condition.num_blocks[shape])
  for(buffer=0;buffer<level->boundary_condition.num_blocks[shape];buffer++){
    int i,j,k;
    const int       box = level->boundary_condition.blocks[shape][buffer].read.box; 
    const int     dim_i = level->boundary_condition.blocks[shape][buffer].dim.i;
    const int     dim_j = level->boundary_condition.blocks[shape][buffer].dim.j;
    const int     dim_k = level->boundary_condition.blocks[shape][buffer].dim.k;
    const int       ilo = level->boundary_condition.blocks[shape][buffer].read.i;
    const int       jlo = level->boundary_condition.blocks[shape][buffer].read.j;
    const int       klo = level->boundary_condition.blocks[shape][buffer].read.k;

    // total hack/reuse of the existing boundary list...
    //   however, whereas boundary subtype represents the normal to the domain at that point, 
    //   one needs the box-relative (not domain-relative) normal when extending the face averaged beta's into the ghost zones
    //   Thus, I reuse the list to tell me which areas are beyond the domain boundary, but must calculate their normals here
          int   subtype = 13;
    if(ilo <               0)subtype-=1;
    if(jlo <               0)subtype-=3;
    if(klo <               0)subtype-=9;
    if(ilo >= level->box_dim)subtype+=1;
    if(jlo >= level->box_dim)subtype+=3;
    if(klo >= level->box_dim)subtype+=9;
    const int    normal = 26-subtype; // invert the normal vector
 
    // hard code for box to box BC's 
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    double * __restrict__  beta_i = level->my_boxes[box].vectors[VECTOR_BETA_I] + level->box_ghosts*(1+jStride+kStride);
    double * __restrict__  beta_j = level->my_boxes[box].vectors[VECTOR_BETA_J] + level->box_ghosts*(1+jStride+kStride);
    double * __restrict__  beta_k = level->my_boxes[box].vectors[VECTOR_BETA_K] + level->box_ghosts*(1+jStride+kStride);

    // convert normal vector into pointer offsets...
    const int di = (((normal % 3)  )-1);
    const int dj = (((normal % 9)/3)-1);
    const int dk = (((normal / 9)  )-1);

    // beta_i should be extrapolated in the j- and k-directions, but not i
    // beta_j should be extrapolated in the i- and k-directions, but not j
    // beta_k should be extrapolated in the i- and j-directions, but not k
    // e.g.
    //                  .................................
    //                 .       .       .       .       .
    //                .       .  ???  .  ???  .       .
    //               .       .       .       .       .
    //              ........+-------+-------+........
    //             .       /       /       /       .
    //            .  ???  /<betaK>/<betaK>/  ???  .
    //           .       /       /       /       .
    //          ........+-------+-------+........
    //         .       /       /       /       .
    //        .  ???  /<betaK>/<betaK>/  ???  .
    //       .       /       /       /       .
    //      ........+-------+-------+........   k   j
    //     .       .       .       .       .    ^  ^   
    //    .       .  ???  .  ???  .       .     | /
    //   .       .       .       .       .      |/
    //  .................................       +-----> i
    //
    const int biStride =      dj*jStride + dk*kStride;
    const int bjStride = di              + dk*kStride;
    const int bkStride = di + dj*jStride             ;

    // note, 
    //   the face values normal to i should have been filled via RESTRICT_I (skip them)
    //   the face values normal to j should have been filled via RESTRICT_J (skip them)
    //   the face values normal to k should have been filled via RESTRICT_K (skip them)
    if(level->box_dim>=5){
      // quartic extrapolation... 
      for(k=0;k<dim_k;k++){
      for(j=0;j<dim_j;j++){
      for(i=0;i<dim_i;i++){
        int ijk = (i+ilo) + (j+jlo)*jStride + (k+klo)*kStride;
        if( (subtype!=14) && (subtype!=12) ){beta_i[ijk] = 5.0*beta_i[ijk+biStride] - 10.0*beta_i[ijk+2*biStride] + 10.0*beta_i[ijk+3*biStride] - 5.0*beta_i[ijk+4*biStride] + beta_i[ijk+5*biStride];}
        if( (subtype!=16) && (subtype!=10) ){beta_j[ijk] = 5.0*beta_j[ijk+bjStride] - 10.0*beta_j[ijk+2*bjStride] + 10.0*beta_j[ijk+3*bjStride] - 5.0*beta_j[ijk+4*bjStride] + beta_j[ijk+5*bjStride];}
        if( (subtype!=22) && (subtype!= 4) ){beta_k[ijk] = 5.0*beta_k[ijk+bkStride] - 10.0*beta_k[ijk+2*bkStride] + 10.0*beta_k[ijk+3*bkStride] - 5.0*beta_k[ijk+4*bkStride] + beta_k[ijk+5*bkStride];}
      }}}
    }else 
    if(level->box_dim>=4){
      // cubic extrapolation... 
      for(k=0;k<dim_k;k++){
      for(j=0;j<dim_j;j++){
      for(i=0;i<dim_i;i++){
        int ijk = (i+ilo) + (j+jlo)*jStride + (k+klo)*kStride;
        if( (subtype!=14) && (subtype!=12) ){beta_i[ijk] = 4.0*beta_i[ijk+biStride] - 6.0*beta_i[ijk+2*biStride] + 4.0*beta_i[ijk+3*biStride] - beta_i[ijk+4*biStride];}
        if( (subtype!=16) && (subtype!=10) ){beta_j[ijk] = 4.0*beta_j[ijk+bjStride] - 6.0*beta_j[ijk+2*bjStride] + 4.0*beta_j[ijk+3*bjStride] - beta_j[ijk+4*bjStride];}
        if( (subtype!=22) && (subtype!= 4) ){beta_k[ijk] = 4.0*beta_k[ijk+bkStride] - 6.0*beta_k[ijk+2*bkStride] + 4.0*beta_k[ijk+3*bkStride] - beta_k[ijk+4*bkStride];}
      }}}
    }else 
    if(level->box_dim>=2){
      // linear extrapolation...
      for(k=0;k<dim_k;k++){
      for(j=0;j<dim_j;j++){
      for(i=0;i<dim_i;i++){
        int ijk = (i+ilo) + (j+jlo)*jStride + (k+klo)*kStride;
        if( (subtype!=14) && (subtype!=12) ){beta_i[ijk] = 2.0*beta_i[ijk+biStride] - beta_i[ijk+2*biStride];}
        if( (subtype!=16) && (subtype!=10) ){beta_j[ijk] = 2.0*beta_j[ijk+bjStride] - beta_j[ijk+2*bjStride];}
        if( (subtype!=22) && (subtype!= 4) ){beta_k[ijk] = 2.0*beta_k[ijk+bkStride] - beta_k[ijk+2*bkStride];}
      }}}
    }

  }
  level->timers.boundary_conditions += (double)(getTime()-_timeStart);
}

//------------------------------------------------------------------------------------------------------------------------------
