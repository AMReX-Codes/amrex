//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#if 0
#if (_OPENMP>=201107) // OpenMP 3.1 supports max reductions...
  #define PRAGMA_THREAD_ACROSS_BOXES(    level,box)           MyPragma(omp parallel for private(box)   if(level->concurrent_boxes>1) num_threads(level->concurrent_boxes)                                                       )
  #define PRAGMA_THREAD_ACROSS_BOXES_SUM(level,box,level_sum) MyPragma(omp parallel for private(box)   if(level->concurrent_boxes>1) num_threads(level->concurrent_boxes)             reduction(  +:level_sum) schedule(static) )
  #define PRAGMA_THREAD_ACROSS_BOXES_MAX(level,box,level_max) MyPragma(omp parallel for private(box)   if(level->concurrent_boxes>1) num_threads(level->concurrent_boxes)             reduction(max:level_max) schedule(static) )
  #define PRAGMA_THREAD_WITHIN_A_BOX(    level,i,j,k)         MyPragma(omp parallel for private(i,j,k) if(level->threads_per_box >1) num_threads(level->threads_per_box ) collapse(2)                                           ) 
  #define PRAGMA_THREAD_WITHIN_A_BOX_SUM(level,i,j,k,box_sum) MyPragma(omp parallel for private(i,j,k) if(level->threads_per_box >1) num_threads(level->threads_per_box ) collapse(2) reduction(  +:  box_sum) schedule(static) ) 
  #define PRAGMA_THREAD_WITHIN_A_BOX_MAX(level,i,j,k,box_max) MyPragma(omp parallel for private(i,j,k) if(level->threads_per_box >1) num_threads(level->threads_per_box ) collapse(2) reduction(max:  box_max) schedule(static) ) 
#elif _OPENMP // older OpenMP versions don't support the max reduction clause
  #define PRAGMA_THREAD_ACROSS_BOXES(    level,box)           MyPragma(omp parallel for private(box)   if(level->concurrent_boxes>1) num_threads(level->concurrent_boxes) )
  #define PRAGMA_THREAD_ACROSS_BOXES_SUM(level,box,level_sum) MyPragma(omp parallel for private(box)   if(level->concurrent_boxes>1) num_threads(level->concurrent_boxes)             reduction(  +:level_sum) schedule(static) )
  #define PRAGMA_THREAD_ACROSS_BOXES_MAX(level,box,level_max) #warning Threading max reductions requires OpenMP 3.1 (July 2011).  Please upgrade your compiler.
  #define PRAGMA_THREAD_WITHIN_A_BOX(    level,i,j,k)         MyPragma(omp parallel for private(i,j,k) if(level->threads_per_box >1) num_threads(level->threads_per_box ) collapse(2) ) 
  #define PRAGMA_THREAD_WITHIN_A_BOX_SUM(level,i,j,k,box_sum) MyPragma(omp parallel for private(i,j,k) if(level->threads_per_box >1) num_threads(level->threads_per_box ) collapse(2) reduction(  +:  box_sum) schedule(static) ) 
  #define PRAGMA_THREAD_WITHIN_A_BOX_MAX(level,i,j,k,box_max) #warning Threading max reductions requires OpenMP 3.1 (July 2011).  Please upgrade your compiler.
#else // flat MPI should not define any threading...
  #define PRAGMA_THREAD_ACROSS_BOXES(    level,box)          
  #define PRAGMA_THREAD_ACROSS_BOXES_SUM(level,box,level_sum)
  #define PRAGMA_THREAD_ACROSS_BOXES_MAX(level,box,level_max)
  #define PRAGMA_THREAD_WITHIN_A_BOX(    level,i,j,k)        
  #define PRAGMA_THREAD_WITHIN_A_BOX_SUM(level,i,j,k,box_sum)
  #define PRAGMA_THREAD_WITHIN_A_BOX_MAX(level,i,j,k,box_max)
#endif
#else
#if (_OPENMP>=201107) // OpenMP 3.1 supports max reductions...
  #define PRAGMA_THREAD_ACROSS_BOXES(    level,box)           
  #define PRAGMA_THREAD_ACROSS_BOXES_SUM(level,box,level_sum) 
  #define PRAGMA_THREAD_ACROSS_BOXES_MAX(level,box,level_max) 
  #define PRAGMA_THREAD_WITHIN_A_BOX(    level,i,j,k)         MyPragma(omp parallel for private(i,j,k) collapse(2)                                           ) 
  #define PRAGMA_THREAD_WITHIN_A_BOX_SUM(level,i,j,k,box_sum) MyPragma(omp parallel for private(i,j,k) collapse(2) reduction(  +:  box_sum) schedule(static) ) 
  #define PRAGMA_THREAD_WITHIN_A_BOX_MAX(level,i,j,k,box_max) MyPragma(omp parallel for private(i,j,k) collapse(2) reduction(max:  box_max) schedule(static) ) 
#elif _OPENMP // older OpenMP versions don't support the max reduction clause
  #define PRAGMA_THREAD_ACROSS_BOXES(    level,box)           
  #define PRAGMA_THREAD_ACROSS_BOXES_SUM(level,box,level_sum) 
  #define PRAGMA_THREAD_ACROSS_BOXES_MAX(level,box,level_max) #warning Threading max reductions requires OpenMP 3.1 (July 2011).  Please upgrade your compiler.
  #define PRAGMA_THREAD_WITHIN_A_BOX(    level,i,j,k)         MyPragma(omp parallel for private(i,j,k) collapse(2)                                           ) 
  #define PRAGMA_THREAD_WITHIN_A_BOX_SUM(level,i,j,k,box_sum) MyPragma(omp parallel for private(i,j,k) collapse(2) reduction(  +:  box_sum) schedule(static) ) 
  #define PRAGMA_THREAD_WITHIN_A_BOX_MAX(level,i,j,k,box_max) #warning Threading max reductions requires OpenMP 3.1 (July 2011).  Please upgrade your compiler.
#else // flat MPI should not define any threading...
  #define PRAGMA_THREAD_ACROSS_BOXES(    level,box)          
  #define PRAGMA_THREAD_ACROSS_BOXES_SUM(level,box,level_sum)
  #define PRAGMA_THREAD_ACROSS_BOXES_MAX(level,box,level_max)
  #define PRAGMA_THREAD_WITHIN_A_BOX(    level,i,j,k)        
  #define PRAGMA_THREAD_WITHIN_A_BOX_SUM(level,i,j,k,box_sum)
  #define PRAGMA_THREAD_WITHIN_A_BOX_MAX(level,i,j,k,box_max)
#endif
#endif
//------------------------------------------------------------------------------------------------------------------------------
