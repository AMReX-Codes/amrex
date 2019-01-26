//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <math.h>
//------------------------------------------------------------------------------------------------------------------------------
static inline void interpolation_v4_block(level_type *level_f, int id_f, double prescale_f, level_type *level_c, int id_c, blockCopy_type *block){
  // interpolate 3D array from read_i,j,k of read[] to write_i,j,k in write[] using volume averaged quartic prolongation
  int write_dim_i   = block->dim.i<<1; // calculate the dimensions of the resultant fine block
  int write_dim_j   = block->dim.j<<1;
  int write_dim_k   = block->dim.k<<1;

  int  read_i       = block->read.i;
  int  read_j       = block->read.j;
  int  read_k       = block->read.k;
  int  read_jStride = block->read.jStride;
  int  read_kStride = block->read.kStride;

  int write_i       = block->write.i;
  int write_j       = block->write.j;
  int write_k       = block->write.k;
  int write_jStride = block->write.jStride;
  int write_kStride = block->write.kStride;

  const double * __restrict__  read = block->read.ptr;
        double * __restrict__ write = block->write.ptr;

  if(block->read.box >=0){
     read_jStride = level_c->my_boxes[block->read.box ].jStride;
     read_kStride = level_c->my_boxes[block->read.box ].kStride;
     read = level_c->my_boxes[ block->read.box].vectors[id_c] + level_c->box_ghosts*(1+ read_jStride+ read_kStride);
  }
  if(block->write.box>=0){
    write_jStride = level_f->my_boxes[block->write.box].jStride;
    write_kStride = level_f->my_boxes[block->write.box].kStride;
    write = level_f->my_boxes[block->write.box].vectors[id_f] + level_f->box_ghosts*(1+write_jStride+write_kStride);
  }
 

  #ifdef USE_NAIVE_INTERP
  // naive 125pt per fine grid cell
  int i,j,k;
  double c2 = -3.0/128.0;
  double c1 = 22.0/128.0;
  int dj  =   read_jStride;
  int dk  =   read_kStride;
  int dj2 = 2*read_jStride;
  int dk2 = 2*read_kStride;
  for(k=0;k<write_dim_k;k++){double sk1=c1,sk2=c2;if(k&0x1){sk1=-c1;sk2=-c2;}
  for(j=0;j<write_dim_j;j++){double sj1=c1,sj2=c2;if(j&0x1){sj1=-c1;sj2=-c2;}
  for(i=0;i<write_dim_i;i++){double si1=c1,si2=c2;if(i&0x1){si1=-c1;si2=-c2;}
    int write_ijk = ((i   )+write_i) + (((j   )+write_j)*write_jStride) + (((k   )+write_k)*write_kStride);
    int  read_ijk = ((i>>1)+ read_i) + (((j>>1)+ read_j)* read_jStride) + (((k>>1)+ read_k)* read_kStride);
    //
    // |   -3/128  |  +22/128  |    1.0    |  -22/128  |   +3/128  | coarse grid
    // |-----+-----|-----+-----|-----+-----|-----+-----|-----+-----|
    // |     |     |     |     |?????|     |     |     |     |     | fine grid
    //
    write[write_ijk] = prescale_f*write[write_ijk] +
                       + sk2*( + sj2*( si2*read[read_ijk-2-dj2-dk2] + si1*read[read_ijk-1-dj2-dk2] + read[read_ijk-dj2-dk2] - si1*read[read_ijk+1-dj2-dk2] - si2*read[read_ijk+2-dj2-dk2] )
                               + sj1*( si2*read[read_ijk-2-dj -dk2] + si1*read[read_ijk-1-dj -dk2] + read[read_ijk-dj -dk2] - si1*read[read_ijk+1-dj -dk2] - si2*read[read_ijk+2-dj -dk2] )
                               +     ( si2*read[read_ijk-2    -dk2] + si1*read[read_ijk-1    -dk2] + read[read_ijk    -dk2] - si1*read[read_ijk+1    -dk2] - si2*read[read_ijk+2    -dk2] )
                               - sj1*( si2*read[read_ijk-2+dj -dk2] + si1*read[read_ijk-1+dj -dk2] + read[read_ijk+dj -dk2] - si1*read[read_ijk+1+dj -dk2] - si2*read[read_ijk+2+dj -dk2] )
                               - sj2*( si2*read[read_ijk-2+dj2-dk2] + si1*read[read_ijk-1+dj2-dk2] + read[read_ijk+dj2-dk2] - si1*read[read_ijk+1+dj2-dk2] - si2*read[read_ijk+2+dj2-dk2] ) )
                       + sk1*( + sj2*( si2*read[read_ijk-2-dj2-dk ] + si1*read[read_ijk-1-dj2-dk ] + read[read_ijk-dj2-dk ] - si1*read[read_ijk+1-dj2-dk ] - si2*read[read_ijk+2-dj2-dk ] )
                               + sj1*( si2*read[read_ijk-2-dj -dk ] + si1*read[read_ijk-1-dj -dk ] + read[read_ijk-dj -dk ] - si1*read[read_ijk+1-dj -dk ] - si2*read[read_ijk+2-dj -dk ] )
                               +     ( si2*read[read_ijk-2    -dk ] + si1*read[read_ijk-1    -dk ] + read[read_ijk    -dk ] - si1*read[read_ijk+1    -dk ] - si2*read[read_ijk+2    -dk ] )
                               - sj1*( si2*read[read_ijk-2+dj -dk ] + si1*read[read_ijk-1+dj -dk ] + read[read_ijk+dj -dk ] - si1*read[read_ijk+1+dj -dk ] - si2*read[read_ijk+2+dj -dk ] )
                               - sj2*( si2*read[read_ijk-2+dj2-dk ] + si1*read[read_ijk-1+dj2-dk ] + read[read_ijk+dj2-dk ] - si1*read[read_ijk+1+dj2-dk ] - si2*read[read_ijk+2+dj2-dk ] ) )
                       +     ( + sj2*( si2*read[read_ijk-2-dj2    ] + si1*read[read_ijk-1-dj2    ] + read[read_ijk-dj2    ] - si1*read[read_ijk+1-dj2    ] - si2*read[read_ijk+2-dj2    ] )
                               + sj1*( si2*read[read_ijk-2-dj     ] + si1*read[read_ijk-1-dj     ] + read[read_ijk-dj     ] - si1*read[read_ijk+1-dj     ] - si2*read[read_ijk+2-dj     ] )
                               +     ( si2*read[read_ijk-2        ] + si1*read[read_ijk-1        ] + read[read_ijk        ] - si1*read[read_ijk+1        ] - si2*read[read_ijk+2        ] )
                               - sj1*( si2*read[read_ijk-2+dj     ] + si1*read[read_ijk-1+dj     ] + read[read_ijk+dj     ] - si1*read[read_ijk+1+dj     ] - si2*read[read_ijk+2+dj     ] )
                               - sj2*( si2*read[read_ijk-2+dj2    ] + si1*read[read_ijk-1+dj2    ] + read[read_ijk+dj2    ] - si1*read[read_ijk+1+dj2    ] - si2*read[read_ijk+2+dj2    ] ) )
                       - sk1*( + sj2*( si2*read[read_ijk-2-dj2+dk ] + si1*read[read_ijk-1-dj2+dk ] + read[read_ijk-dj2+dk ] - si1*read[read_ijk+1-dj2+dk ] - si2*read[read_ijk+2-dj2+dk ] )
                               + sj1*( si2*read[read_ijk-2-dj +dk ] + si1*read[read_ijk-1-dj +dk ] + read[read_ijk-dj +dk ] - si1*read[read_ijk+1-dj +dk ] - si2*read[read_ijk+2-dj +dk ] )
                               +     ( si2*read[read_ijk-2    +dk ] + si1*read[read_ijk-1    +dk ] + read[read_ijk    +dk ] - si1*read[read_ijk+1    +dk ] - si2*read[read_ijk+2    +dk ] )
                               - sj1*( si2*read[read_ijk-2+dj +dk ] + si1*read[read_ijk-1+dj +dk ] + read[read_ijk+dj +dk ] - si1*read[read_ijk+1+dj +dk ] - si2*read[read_ijk+2+dj +dk ] )
                               - sj2*( si2*read[read_ijk-2+dj2+dk ] + si1*read[read_ijk-1+dj2+dk ] + read[read_ijk+dj2+dk ] - si1*read[read_ijk+1+dj2+dk ] - si2*read[read_ijk+2+dj2+dk ] ) )
                       - sk2*( + sj2*( si2*read[read_ijk-2-dj2+dk2] + si1*read[read_ijk-1-dj2+dk2] + read[read_ijk-dj2+dk2] - si1*read[read_ijk+1-dj2+dk2] - si2*read[read_ijk+2-dj2+dk2] )
                               + sj1*( si2*read[read_ijk-2-dj +dk2] + si1*read[read_ijk-1-dj +dk2] + read[read_ijk-dj +dk2] - si1*read[read_ijk+1-dj +dk2] - si2*read[read_ijk+2-dj +dk2] )
                               +     ( si2*read[read_ijk-2    +dk2] + si1*read[read_ijk-1    +dk2] + read[read_ijk    +dk2] - si1*read[read_ijk+1    +dk2] - si2*read[read_ijk+2    +dk2] )
                               - sj1*( si2*read[read_ijk-2+dj +dk2] + si1*read[read_ijk-1+dj +dk2] + read[read_ijk+dj +dk2] - si1*read[read_ijk+1+dj +dk2] - si2*read[read_ijk+2+dj +dk2] )
                               - sj2*( si2*read[read_ijk-2+dj2+dk2] + si1*read[read_ijk-1+dj2+dk2] + read[read_ijk+dj2+dk2] - si1*read[read_ijk+1+dj2+dk2] - si2*read[read_ijk+2+dj2+dk2] ) );
  }}}
  #else
  // exploit tensor product symmetry and perform 8 fine grid interpolations at a time...
  //   50 x 5pt for i
  //   20 x 5pt for j 
  //    8 x 5pt for k
  // ----------------
  //   78 x 5pt for 8 cells (vs 8x125pt = 200x5pt in naive)
  int i,j,k;
  int ii,jj,kk;
  double c2 = -3.0/128.0;
  double c1 = 22.0/128.0;
  int dj  =   read_jStride;
  int dk  =   read_kStride;
  int dj2 = 2*read_jStride;
  int dk2 = 2*read_kStride;
  for(k=0,kk=0;k<write_dim_k;k+=2,kk++){
  for(j=0,jj=0;j<write_dim_j;j+=2,jj++){
  // compiler cannot infer/speculate write[ijk+write_jStride] is disjoint from write[ijk], so create a unique restrict pointers for each nonliteral offset...
  double * __restrict__ write00 = write + write_i + (write_j+j+0)*write_jStride + (write_k+k+0)*write_kStride;
  double * __restrict__ write10 = write + write_i + (write_j+j+1)*write_jStride + (write_k+k+0)*write_kStride;
  double * __restrict__ write01 = write + write_i + (write_j+j+0)*write_jStride + (write_k+k+1)*write_kStride;
  double * __restrict__ write11 = write + write_i + (write_j+j+1)*write_jStride + (write_k+k+1)*write_kStride;
  for(i=0,ii=0;i<write_dim_i;i+=2,ii++){
    int write_ijk = ( i+write_i) + ( j+write_j)*write_jStride + ( k+write_k)*write_kStride;
    int  read_ijk = (ii+ read_i) + (jj+ read_j)* read_jStride + (kk+ read_k)* read_kStride;
    //
    // |   -3/128  |  +22/128  |    1.0    |  -22/128  |   +3/128  | coarse grid
    // |-----+-----|-----+-----|-----+-----|-----+-----|-----+-----|
    // |     |     |     |     |?????|     |     |     |     |     | fine grid
    //

    // grab all coarse grid points...
    const double c000=read[read_ijk-2-dj2-dk2], c100=read[read_ijk-1-dj2-dk2], c200=read[read_ijk-dj2-dk2], c300=read[read_ijk+1-dj2-dk2], c400=read[read_ijk+2-dj2-dk2];
    const double c010=read[read_ijk-2-dj -dk2], c110=read[read_ijk-1-dj -dk2], c210=read[read_ijk-dj -dk2], c310=read[read_ijk+1-dj -dk2], c410=read[read_ijk+2-dj -dk2];
    const double c020=read[read_ijk-2    -dk2], c120=read[read_ijk-1    -dk2], c220=read[read_ijk    -dk2], c320=read[read_ijk+1    -dk2], c420=read[read_ijk+2    -dk2];
    const double c030=read[read_ijk-2+dj -dk2], c130=read[read_ijk-1+dj -dk2], c230=read[read_ijk+dj -dk2], c330=read[read_ijk+1+dj -dk2], c430=read[read_ijk+2+dj -dk2];
    const double c040=read[read_ijk-2+dj2-dk2], c140=read[read_ijk-1+dj2-dk2], c240=read[read_ijk+dj2-dk2], c340=read[read_ijk+1+dj2-dk2], c440=read[read_ijk+2+dj2-dk2];

    const double c001=read[read_ijk-2-dj2-dk ], c101=read[read_ijk-1-dj2-dk ], c201=read[read_ijk-dj2-dk ], c301=read[read_ijk+1-dj2-dk ], c401=read[read_ijk+2-dj2-dk ];
    const double c011=read[read_ijk-2-dj -dk ], c111=read[read_ijk-1-dj -dk ], c211=read[read_ijk-dj -dk ], c311=read[read_ijk+1-dj -dk ], c411=read[read_ijk+2-dj -dk ];
    const double c021=read[read_ijk-2    -dk ], c121=read[read_ijk-1    -dk ], c221=read[read_ijk    -dk ], c321=read[read_ijk+1    -dk ], c421=read[read_ijk+2    -dk ];
    const double c031=read[read_ijk-2+dj -dk ], c131=read[read_ijk-1+dj -dk ], c231=read[read_ijk+dj -dk ], c331=read[read_ijk+1+dj -dk ], c431=read[read_ijk+2+dj -dk ];
    const double c041=read[read_ijk-2+dj2-dk ], c141=read[read_ijk-1+dj2-dk ], c241=read[read_ijk+dj2-dk ], c341=read[read_ijk+1+dj2-dk ], c441=read[read_ijk+2+dj2-dk ];

    const double c002=read[read_ijk-2-dj2    ], c102=read[read_ijk-1-dj2    ], c202=read[read_ijk-dj2    ], c302=read[read_ijk+1-dj2    ], c402=read[read_ijk+2-dj2    ];
    const double c012=read[read_ijk-2-dj     ], c112=read[read_ijk-1-dj     ], c212=read[read_ijk-dj     ], c312=read[read_ijk+1-dj     ], c412=read[read_ijk+2-dj     ];
    const double c022=read[read_ijk-2        ], c122=read[read_ijk-1        ], c222=read[read_ijk        ], c322=read[read_ijk+1        ], c422=read[read_ijk+2        ];
    const double c032=read[read_ijk-2+dj     ], c132=read[read_ijk-1+dj     ], c232=read[read_ijk+dj     ], c332=read[read_ijk+1+dj     ], c432=read[read_ijk+2+dj     ];
    const double c042=read[read_ijk-2+dj2    ], c142=read[read_ijk-1+dj2    ], c242=read[read_ijk+dj2    ], c342=read[read_ijk+1+dj2    ], c442=read[read_ijk+2+dj2    ];

    const double c003=read[read_ijk-2-dj2+dk ], c103=read[read_ijk-1-dj2+dk ], c203=read[read_ijk-dj2+dk ], c303=read[read_ijk+1-dj2+dk ], c403=read[read_ijk+2-dj2+dk ];
    const double c013=read[read_ijk-2-dj +dk ], c113=read[read_ijk-1-dj +dk ], c213=read[read_ijk-dj +dk ], c313=read[read_ijk+1-dj +dk ], c413=read[read_ijk+2-dj +dk ];
    const double c023=read[read_ijk-2    +dk ], c123=read[read_ijk-1    +dk ], c223=read[read_ijk    +dk ], c323=read[read_ijk+1    +dk ], c423=read[read_ijk+2    +dk ];
    const double c033=read[read_ijk-2+dj +dk ], c133=read[read_ijk-1+dj +dk ], c233=read[read_ijk+dj +dk ], c333=read[read_ijk+1+dj +dk ], c433=read[read_ijk+2+dj +dk ];
    const double c043=read[read_ijk-2+dj2+dk ], c143=read[read_ijk-1+dj2+dk ], c243=read[read_ijk+dj2+dk ], c343=read[read_ijk+1+dj2+dk ], c443=read[read_ijk+2+dj2+dk ];

    const double c004=read[read_ijk-2-dj2+dk2], c104=read[read_ijk-1-dj2+dk2], c204=read[read_ijk-dj2+dk2], c304=read[read_ijk+1-dj2+dk2], c404=read[read_ijk+2-dj2+dk2];
    const double c014=read[read_ijk-2-dj +dk2], c114=read[read_ijk-1-dj +dk2], c214=read[read_ijk-dj +dk2], c314=read[read_ijk+1-dj +dk2], c414=read[read_ijk+2-dj +dk2];
    const double c024=read[read_ijk-2    +dk2], c124=read[read_ijk-1    +dk2], c224=read[read_ijk    +dk2], c324=read[read_ijk+1    +dk2], c424=read[read_ijk+2    +dk2];
    const double c034=read[read_ijk-2+dj +dk2], c134=read[read_ijk-1+dj +dk2], c234=read[read_ijk+dj +dk2], c334=read[read_ijk+1+dj +dk2], c434=read[read_ijk+2+dj +dk2];
    const double c044=read[read_ijk-2+dj2+dk2], c144=read[read_ijk-1+dj2+dk2], c244=read[read_ijk+dj2+dk2], c344=read[read_ijk+1+dj2+dk2], c444=read[read_ijk+2+dj2+dk2];

    // interpolate in i to create fine i / coarse jk points...
    const double f0c00 = ( c200 + c1*(c100-c300) + c2*(c000-c400) ); // same as original 5pt stencil...  f0c00 = ( c2*c000 + c1*c100 + c200 - c1*c300 - c2*c400 )
    const double f1c00 = ( c200 - c1*(c100-c300) - c2*(c000-c400) );
    const double f0c10 = ( c210 + c1*(c110-c310) + c2*(c010-c410) );
    const double f1c10 = ( c210 - c1*(c110-c310) - c2*(c010-c410) );
    const double f0c20 = ( c220 + c1*(c120-c320) + c2*(c020-c420) );
    const double f1c20 = ( c220 - c1*(c120-c320) - c2*(c020-c420) );
    const double f0c30 = ( c230 + c1*(c130-c330) + c2*(c030-c430) );
    const double f1c30 = ( c230 - c1*(c130-c330) - c2*(c030-c430) );
    const double f0c40 = ( c240 + c1*(c140-c340) + c2*(c040-c440) );
    const double f1c40 = ( c240 - c1*(c140-c340) - c2*(c040-c440) );

    const double f0c01 = ( c201 + c1*(c101-c301) + c2*(c001-c401) );
    const double f1c01 = ( c201 - c1*(c101-c301) - c2*(c001-c401) );
    const double f0c11 = ( c211 + c1*(c111-c311) + c2*(c011-c411) );
    const double f1c11 = ( c211 - c1*(c111-c311) - c2*(c011-c411) );
    const double f0c21 = ( c221 + c1*(c121-c321) + c2*(c021-c421) );
    const double f1c21 = ( c221 - c1*(c121-c321) - c2*(c021-c421) );
    const double f0c31 = ( c231 + c1*(c131-c331) + c2*(c031-c431) );
    const double f1c31 = ( c231 - c1*(c131-c331) - c2*(c031-c431) );
    const double f0c41 = ( c241 + c1*(c141-c341) + c2*(c041-c441) );
    const double f1c41 = ( c241 - c1*(c141-c341) - c2*(c041-c441) );

    const double f0c02 = ( c202 + c1*(c102-c302) + c2*(c002-c402) );
    const double f1c02 = ( c202 - c1*(c102-c302) - c2*(c002-c402) );
    const double f0c12 = ( c212 + c1*(c112-c312) + c2*(c012-c412) );
    const double f1c12 = ( c212 - c1*(c112-c312) - c2*(c012-c412) );
    const double f0c22 = ( c222 + c1*(c122-c322) + c2*(c022-c422) );
    const double f1c22 = ( c222 - c1*(c122-c322) - c2*(c022-c422) );
    const double f0c32 = ( c232 + c1*(c132-c332) + c2*(c032-c432) );
    const double f1c32 = ( c232 - c1*(c132-c332) - c2*(c032-c432) );
    const double f0c42 = ( c242 + c1*(c142-c342) + c2*(c042-c442) );
    const double f1c42 = ( c242 - c1*(c142-c342) - c2*(c042-c442) );

    const double f0c03 = ( c203 + c1*(c103-c303) + c2*(c003-c403) );
    const double f1c03 = ( c203 - c1*(c103-c303) - c2*(c003-c403) );
    const double f0c13 = ( c213 + c1*(c113-c313) + c2*(c013-c413) );
    const double f1c13 = ( c213 - c1*(c113-c313) - c2*(c013-c413) );
    const double f0c23 = ( c223 + c1*(c123-c323) + c2*(c023-c423) );
    const double f1c23 = ( c223 - c1*(c123-c323) - c2*(c023-c423) );
    const double f0c33 = ( c233 + c1*(c133-c333) + c2*(c033-c433) );
    const double f1c33 = ( c233 - c1*(c133-c333) - c2*(c033-c433) );
    const double f0c43 = ( c243 + c1*(c143-c343) + c2*(c043-c443) );
    const double f1c43 = ( c243 - c1*(c143-c343) - c2*(c043-c443) );

    const double f0c04 = ( c204 + c1*(c104-c304) + c2*(c004-c404) );
    const double f1c04 = ( c204 - c1*(c104-c304) - c2*(c004-c404) );
    const double f0c14 = ( c214 + c1*(c114-c314) + c2*(c014-c414) );
    const double f1c14 = ( c214 - c1*(c114-c314) - c2*(c014-c414) );
    const double f0c24 = ( c224 + c1*(c124-c324) + c2*(c024-c424) );
    const double f1c24 = ( c224 - c1*(c124-c324) - c2*(c024-c424) );
    const double f0c34 = ( c234 + c1*(c134-c334) + c2*(c034-c434) );
    const double f1c34 = ( c234 - c1*(c134-c334) - c2*(c034-c434) );
    const double f0c44 = ( c244 + c1*(c144-c344) + c2*(c044-c444) );
    const double f1c44 = ( c244 - c1*(c144-c344) - c2*(c044-c444) );

    // interpolate in j to create fine ij / coarse k points...
    const double f00c0 = (f0c20 + c1*(f0c10-f0c30) + c2*(f0c00-f0c40) );
    const double f10c0 = (f1c20 + c1*(f1c10-f1c30) + c2*(f1c00-f1c40) );
    const double f01c0 = (f0c20 - c1*(f0c10-f0c30) - c2*(f0c00-f0c40) );
    const double f11c0 = (f1c20 - c1*(f1c10-f1c30) - c2*(f1c00-f1c40) );

    const double f00c1 = (f0c21 + c1*(f0c11-f0c31) + c2*(f0c01-f0c41) );
    const double f10c1 = (f1c21 + c1*(f1c11-f1c31) + c2*(f1c01-f1c41) );
    const double f01c1 = (f0c21 - c1*(f0c11-f0c31) - c2*(f0c01-f0c41) );
    const double f11c1 = (f1c21 - c1*(f1c11-f1c31) - c2*(f1c01-f1c41) );

    const double f00c2 = (f0c22 + c1*(f0c12-f0c32) + c2*(f0c02-f0c42) );
    const double f10c2 = (f1c22 + c1*(f1c12-f1c32) + c2*(f1c02-f1c42) );
    const double f01c2 = (f0c22 - c1*(f0c12-f0c32) - c2*(f0c02-f0c42) );
    const double f11c2 = (f1c22 - c1*(f1c12-f1c32) - c2*(f1c02-f1c42) );

    const double f00c3 = (f0c23 + c1*(f0c13-f0c33) + c2*(f0c03-f0c43) );
    const double f10c3 = (f1c23 + c1*(f1c13-f1c33) + c2*(f1c03-f1c43) );
    const double f01c3 = (f0c23 - c1*(f0c13-f0c33) - c2*(f0c03-f0c43) );
    const double f11c3 = (f1c23 - c1*(f1c13-f1c33) - c2*(f1c03-f1c43) );

    const double f00c4 = (f0c24 + c1*(f0c14-f0c34) + c2*(f0c04-f0c44) );
    const double f10c4 = (f1c24 + c1*(f1c14-f1c34) + c2*(f1c04-f1c44) );
    const double f01c4 = (f0c24 - c1*(f0c14-f0c34) - c2*(f0c04-f0c44) );
    const double f11c4 = (f1c24 - c1*(f1c14-f1c34) - c2*(f1c04-f1c44) );

    // interpolate in k to create fine ijk points...
    const double f000 = (f00c2 + c1*(f00c1-f00c3) + c2*(f00c0-f00c4) );
    const double f100 = (f10c2 + c1*(f10c1-f10c3) + c2*(f10c0-f10c4) );
    const double f010 = (f01c2 + c1*(f01c1-f01c3) + c2*(f01c0-f01c4) );
    const double f110 = (f11c2 + c1*(f11c1-f11c3) + c2*(f11c0-f11c4) );
    const double f001 = (f00c2 - c1*(f00c1-f00c3) - c2*(f00c0-f00c4) );
    const double f101 = (f10c2 - c1*(f10c1-f10c3) - c2*(f10c0-f10c4) );
    const double f011 = (f01c2 - c1*(f01c1-f01c3) - c2*(f01c0-f01c4) );
    const double f111 = (f11c2 - c1*(f11c1-f11c3) - c2*(f11c0-f11c4) );

    // commit to memory...
    #if 0 // compiler cannot infer/speculate write[ijk+write_jStride] is disjoint from write[ijk], and thus cannot vectorize...
    write[write_ijk                              ] = prescale_f*write[write_ijk                              ] + f000;
    write[write_ijk+1                            ] = prescale_f*write[write_ijk+1                            ] + f100;
    write[write_ijk  +write_jStride              ] = prescale_f*write[write_ijk  +write_jStride              ] + f010;
    write[write_ijk+1+write_jStride              ] = prescale_f*write[write_ijk+1+write_jStride              ] + f110;
    write[write_ijk                +write_kStride] = prescale_f*write[write_ijk                +write_kStride] + f001;
    write[write_ijk+1              +write_kStride] = prescale_f*write[write_ijk+1              +write_kStride] + f101;
    write[write_ijk  +write_jStride+write_kStride] = prescale_f*write[write_ijk  +write_jStride+write_kStride] + f011;
    write[write_ijk+1+write_jStride+write_kStride] = prescale_f*write[write_ijk+1+write_jStride+write_kStride] + f111;
    #else // use a unique restrict pointer for each pencil...
    write00[i  ] = prescale_f*write00[i  ] + f000;
    write00[i+1] = prescale_f*write00[i+1] + f100;
    write10[i  ] = prescale_f*write10[i  ] + f010;
    write10[i+1] = prescale_f*write10[i+1] + f110;
    write01[i  ] = prescale_f*write01[i  ] + f001;
    write01[i+1] = prescale_f*write01[i+1] + f101;
    write11[i  ] = prescale_f*write11[i  ] + f011;
    write11[i+1] = prescale_f*write11[i+1] + f111;
    #endif

  }}}
  #endif

}


//------------------------------------------------------------------------------------------------------------------------------
// perform a (inter-level) volumetric quartic interpolation on vector id_c of the coarse level and increments prescale_f*vector id_f on the fine level by the result
// i.e. id_f = prescale_f*id_f + P*id_c
// prescale_f is nominally 1.0 or 0.0
// quartic interpolation requires a full ghost zone exchange and boundary condition
// This is a rather bulk synchronous implementation which packs all MPI buffers before initiating any sends
// Similarly, it waits for all remote data before copying any into local boxes.
// It does however attempt to overlap local interpolation with MPI
void interpolation_v4(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c){
    exchange_boundary(level_c,id_c,STENCIL_SHAPE_BOX);
         apply_BCs_v4(level_c,id_c,STENCIL_SHAPE_BOX);

  double _timeCommunicationStart = getTime();
  double _timeStart,_timeEnd;
  int buffer=0;
  int n;
  int my_tag = (level_f->tag<<4) | 0x7;


  #ifdef USE_MPI
  // by convention, level_f allocates a combined array of requests for both level_f recvs and level_c sends...
  int nMessages = level_c->interpolation.num_sends + level_f->interpolation.num_recvs;
  MPI_Request *recv_requests = level_f->interpolation.requests;
  MPI_Request *send_requests = level_f->interpolation.requests + level_f->interpolation.num_recvs;


  // loop through packed list of MPI receives and prepost Irecv's...
  if(level_f->interpolation.num_recvs>0){
    _timeStart = getTime();
    #ifdef USE_MPI_THREAD_MULTIPLE
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for(n=0;n<level_f->interpolation.num_recvs;n++){
      MPI_Irecv(level_f->interpolation.recv_buffers[n],
                level_f->interpolation.recv_sizes[n],
                MPI_DOUBLE,
                level_f->interpolation.recv_ranks[n],
                my_tag,
                MPI_COMM_WORLD,
                &recv_requests[n]
      );
    }
    _timeEnd = getTime();
    level_f->timers.interpolation_recv += (_timeEnd-_timeStart);
  }


  // pack MPI send buffers...
  if(level_c->interpolation.num_blocks[0]>0){
    _timeStart = getTime();
    PRAGMA_THREAD_ACROSS_BLOCKS(level_f,buffer,level_c->interpolation.num_blocks[0])
    for(buffer=0;buffer<level_c->interpolation.num_blocks[0];buffer++){
      // !!! prescale==0 because you don't want to increment the MPI buffer
      interpolation_v4_block(level_f,id_f,0.0,level_c,id_c,&level_c->interpolation.blocks[0][buffer]);
    }
    _timeEnd = getTime();
    level_f->timers.interpolation_pack += (_timeEnd-_timeStart);
  }


  // loop through MPI send buffers and post Isend's...
  if(level_c->interpolation.num_sends>0){
    _timeStart = getTime();
    #ifdef USE_MPI_THREAD_MULTIPLE
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for(n=0;n<level_c->interpolation.num_sends;n++){
      MPI_Isend(level_c->interpolation.send_buffers[n],
                level_c->interpolation.send_sizes[n],
                MPI_DOUBLE,
                level_c->interpolation.send_ranks[n],
                my_tag,
                MPI_COMM_WORLD,
                &send_requests[n]
      );
    }
    _timeEnd = getTime();
    level_f->timers.interpolation_send += (_timeEnd-_timeStart);
  }
  #endif


  // perform local interpolation... try and hide within Isend latency... 
  if(level_c->interpolation.num_blocks[1]>0){
    _timeStart = getTime();
    PRAGMA_THREAD_ACROSS_BLOCKS(level_f,buffer,level_c->interpolation.num_blocks[1])
    for(buffer=0;buffer<level_c->interpolation.num_blocks[1];buffer++){
      interpolation_v4_block(level_f,id_f,prescale_f,level_c,id_c,&level_c->interpolation.blocks[1][buffer]);
    }
    _timeEnd = getTime();
    level_f->timers.interpolation_local += (_timeEnd-_timeStart);
  }


  // wait for MPI to finish...
  #ifdef USE_MPI 
  if(nMessages>0){
    _timeStart = getTime();
    MPI_Waitall(nMessages,level_f->interpolation.requests,level_f->interpolation.status);
    _timeEnd = getTime();
    level_f->timers.interpolation_wait += (_timeEnd-_timeStart);
  }


  // unpack MPI receive buffers 
  if(level_f->interpolation.num_blocks[2]>0){
    _timeStart = getTime();
    PRAGMA_THREAD_ACROSS_BLOCKS(level_f,buffer,level_f->interpolation.num_blocks[2])
    for(buffer=0;buffer<level_f->interpolation.num_blocks[2];buffer++){
      IncrementBlock(level_f,id_f,prescale_f,&level_f->interpolation.blocks[2][buffer]);
    }
    _timeEnd = getTime();
    level_f->timers.interpolation_unpack += (_timeEnd-_timeStart);
  }
  #endif 
 
 
  level_f->timers.interpolation_total += (double)(getTime()-_timeCommunicationStart);
}
