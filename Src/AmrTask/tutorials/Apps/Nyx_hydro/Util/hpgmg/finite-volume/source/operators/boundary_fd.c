//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
void apply_BCs_p1(level_type * level, int x_id, int shape){
  // For cell-centered, we need to fill in the ghost zones to apply any BC's
  // This code does a simple piecewise linear interpolation for homogeneous dirichlet (0 on boundary)
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
    double * __restrict__  x = level->my_boxes[box].vectors[x_id] + level->my_boxes[box].ghosts*(1+jStride+kStride);

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
void apply_BCs_p2(level_type * level, int x_id, int shape){
  // For cell-centered, we need to fill in the ghost zones to apply any BC's
  // This code does a simple piecewise quadratic interpolation for homogeneous dirichlet (0 on boundary)
  // Nominally, this is first performed across faces, then to edges, then to corners.  
  //
  if(shape>=STENCIL_MAX_SHAPES)shape=STENCIL_SHAPE_BOX;  // shape must be < STENCIL_MAX_SHAPES in order to safely index into boundary_condition.blocks[]
  if(level->boundary_condition.type == BC_PERIODIC)return; // no BC's to apply !
  if(level->box_dim<2){apply_BCs_p1(level,x_id,shape);return;}

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
    const int normal = 26-level->boundary_condition.blocks[shape][buffer].subtype; // invert the normal vector
 
    // hard code for box to box BC's 
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    double * __restrict__  x = level->my_boxes[box].vectors[x_id] + level->my_boxes[box].ghosts*(1+jStride+kStride);

    // convert normal vector into pointer offsets...
    const int di = (((normal % 3)  )-1)*1;
    const int dj = (((normal % 9)/3)-1)*jStride;
    const int dk = (((normal / 9)  )-1)*kStride;

    if(faces[normal]){
      //
      //      /------/------/------/
      //     /  ??  /  -2  /  1/3 /
      //    /------/------/------/
      //
      const int stride = di+dj+dk;
      const int stride2 = stride*2;
      for(k=0;k<dim_k;k++){
      for(j=0;j<dim_j;j++){
      for(i=0;i<dim_i;i++){
        int ijk = (i+ilo) + (j+jlo)*jStride + (k+klo)*kStride;
        x[ijk] = -2.0*x[ijk+stride] + 0.333333333333333333*x[ijk+stride2]; // 2pt stencil
      }}}
    }else if(edges[normal]){
      //
      //                 /------/------/
      //                / -2/3 /  1/9 /
      //               /------/------/
      //              /   4  / -2/3 /
      //      /------/------/------/
      //     /  ??  /
      //    /------/
      //
      int dr=-1;
      int ds=-1;
      if(di==0){dr=dj;ds=dk;}
      if(dj==0){dr=di;ds=dk;}
      if(dk==0){dr=di;ds=dj;}
      for(k=0;k<dim_k;k++){
      for(j=0;j<dim_j;j++){
      for(i=0;i<dim_i;i++){
        // 4pt stencil...
        int ijk = (i+ilo) + (j+jlo)*jStride + (k+klo)*kStride;
        x[ijk] =   4.000000000000000000*x[ijk+  dr+  ds] 
                 - 0.666666666666666667*x[ijk+2*dr+  ds]
                 - 0.666666666666666667*x[ijk+  dr+2*ds]
                 + 0.111111111111111111*x[ijk+2*dr+2*ds];
      }}}
    }else if(corners[normal]){
      //
      //              /------/------/
      //             / -2/9 / 1/27 /.
      //            /------/------/ .
      //           /  4/3 / -2/9 /  .
      //          /------/------/   .
      //          .                 .
      //          .   /------/------/
      //          .  /  4/3 / -2/9 /.
      //          . /------/------/ .
      //          ./  -8  /  4/3 /  .
      //          /------/------/   .
      //          .                 .
      //          .   /------/------/
      //          .  /   0  /   0  /
      //          . /------/------/
      //          ./   0  /   0  /
      //   /------/------/------/
      //  /  ??  /
      // /------/
      //
      // 4pt stencil...
      int ijk = (ilo) + (jlo)*jStride + (klo)*kStride;
      x[ijk] =  -8.000000000000000000*x[ijk+  di+  dj+  dk] 
                +1.333333333333333333*x[ijk+2*di+  dj+  dk] 
                +1.333333333333333333*x[ijk+  di+2*dj+  dk] 
                +1.333333333333333333*x[ijk+  di+  dj+2*dk] 
                -0.222222222222222222*x[ijk+2*di+2*dj+  dk] 
                -0.222222222222222222*x[ijk+  di+2*dj+2*dk] 
                -0.222222222222222222*x[ijk+2*di+  dj+2*dk] 
                +0.037037037037037037*x[ijk+2*di+2*dj+2*dk];
    }

  }
  level->timers.boundary_conditions += (double)(getTime()-_timeStart);
}
