#include <mpi.h>
#include <BoxLib/Profiler.H>

#ifdef __cplusplus
extern "C"
#endif
int MPI_Allgather( 
void * sendbuf,
int sendcount,
MPI_Datatype sendtype,
void * recvbuf,
int recvcount,
MPI_Datatype recvtype,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Allgather()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Allgather( sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Allgatherv( 
void * sendbuf,
int sendcount,
MPI_Datatype sendtype,
void * recvbuf,
int * recvcounts,
int * displs,
MPI_Datatype recvtype,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Allgatherv()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Allgatherv( sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Allreduce( 
void * sendbuf,
void * recvbuf,
int count,
MPI_Datatype datatype,
MPI_Op op,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Allreduce()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Allreduce( sendbuf, recvbuf, count, datatype, op, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Alltoall( 
void * sendbuf,
int sendcount,
MPI_Datatype sendtype,
void * recvbuf,
int recvcnt,
MPI_Datatype recvtype,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Alltoall()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Alltoall( sendbuf, sendcount, sendtype, recvbuf, recvcnt, recvtype, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Alltoallv( 
void * sendbuf,
int * sendcnts,
int * sdispls,
MPI_Datatype sendtype,
void * recvbuf,
int * recvcnts,
int * rdispls,
MPI_Datatype recvtype,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Alltoallv()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Alltoallv( sendbuf, sendcnts, sdispls, sendtype, recvbuf, recvcnts, rdispls, recvtype, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Barrier( 
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Barrier()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Barrier( comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Bcast( 
void * buffer,
int count,
MPI_Datatype datatype,
int root,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Bcast()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Bcast( buffer, count, datatype, root, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Gather( 
void * sendbuf,
int sendcnt,
MPI_Datatype sendtype,
void * recvbuf,
int recvcount,
MPI_Datatype recvtype,
int root,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Gather()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Gather( sendbuf, sendcnt, sendtype, recvbuf, recvcount, recvtype, root, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Gatherv( 
void * sendbuf,
int sendcnt,
MPI_Datatype sendtype,
void * recvbuf,
int * recvcnts,
int * displs,
MPI_Datatype recvtype,
int root,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Gatherv()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Gatherv( sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, root, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Op_create( 
MPI_User_function * function,
int commute,
MPI_Op * op )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Op_create()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Op_create( function, commute, op );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Op_free( 
MPI_Op * op )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Op_free()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Op_free( op );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Reduce_scatter( 
void * sendbuf,
void * recvbuf,
int * recvcnts,
MPI_Datatype datatype,
MPI_Op op,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Reduce_scatter()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Reduce_scatter( sendbuf, recvbuf, recvcnts, datatype, op, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Reduce( 
void * sendbuf,
void * recvbuf,
int count,
MPI_Datatype datatype,
MPI_Op op,
int root,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Reduce()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Reduce( sendbuf, recvbuf, count, datatype, op, root, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Scan( 
void * sendbuf,
void * recvbuf,
int count,
MPI_Datatype datatype,
MPI_Op op,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Scan()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Scan( sendbuf, recvbuf, count, datatype, op, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Scatter( 
void * sendbuf,
int sendcnt,
MPI_Datatype sendtype,
void * recvbuf,
int recvcnt,
MPI_Datatype recvtype,
int root,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Scatter()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Scatter( sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Scatterv( 
void * sendbuf,
int * sendcnts,
int * displs,
MPI_Datatype sendtype,
void * recvbuf,
int recvcnt,
MPI_Datatype recvtype,
int root,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Scatterv()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Scatterv( sendbuf, sendcnts, displs, sendtype, recvbuf, recvcnt, recvtype, root, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Attr_delete( 
MPI_Comm comm,
int keyval )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Attr_delete()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Attr_delete( comm, keyval );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Attr_get( 
MPI_Comm comm,
int keyval,
void * attr_value,
int * flag )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Attr_get()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Attr_get( comm, keyval, attr_value, flag );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Attr_put( 
MPI_Comm comm,
int keyval,
void * attr_value )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Attr_put()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Attr_put( comm, keyval, attr_value );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Comm_compare( 
MPI_Comm comm1,
MPI_Comm comm2,
int * result )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Comm_compare()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Comm_compare( comm1, comm2, result );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Comm_create( 
MPI_Comm comm,
MPI_Group group,
MPI_Comm * comm_out )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Comm_create()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Comm_create( comm, group, comm_out );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Comm_dup( 
MPI_Comm comm,
MPI_Comm * comm_out )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Comm_dup()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Comm_dup( comm, comm_out );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Comm_free( 
MPI_Comm * comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Comm_free()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Comm_free( comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Comm_group( 
MPI_Comm comm,
MPI_Group * group )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Comm_group()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Comm_group( comm, group );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Comm_rank( 
MPI_Comm comm,
int * rank )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Comm_rank()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Comm_rank( comm, rank );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Comm_remote_group( 
MPI_Comm comm,
MPI_Group * group )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Comm_remote_group()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Comm_remote_group( comm, group );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Comm_remote_size( 
MPI_Comm comm,
int * size )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Comm_remote_size()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Comm_remote_size( comm, size );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Comm_size( 
MPI_Comm comm,
int * size )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Comm_size()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Comm_size( comm, size );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Comm_split( 
MPI_Comm comm,
int color,
int key,
MPI_Comm * comm_out )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Comm_split()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Comm_split( comm, color, key, comm_out );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Comm_test_inter( 
MPI_Comm comm,
int * flag )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Comm_test_inter()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Comm_test_inter( comm, flag );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Group_compare( 
MPI_Group group1,
MPI_Group group2,
int * result )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Group_compare()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Group_compare( group1, group2, result );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Group_difference( 
MPI_Group group1,
MPI_Group group2,
MPI_Group * group_out )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Group_difference()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Group_difference( group1, group2, group_out );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Group_excl( 
MPI_Group group,
int n,
int * ranks,
MPI_Group * newgroup )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Group_excl()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Group_excl( group, n, ranks, newgroup );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Group_free( 
MPI_Group * group )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Group_free()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Group_free( group );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Group_incl( 
MPI_Group group,
int n,
int * ranks,
MPI_Group * group_out )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Group_incl()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Group_incl( group, n, ranks, group_out );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Group_intersection( 
MPI_Group group1,
MPI_Group group2,
MPI_Group * group_out )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Group_intersection()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Group_intersection( group1, group2, group_out );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Group_rank( 
MPI_Group group,
int * rank )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Group_rank()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Group_rank( group, rank );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Group_range_excl( 
MPI_Group group,
int n,
int ranges[][3],
MPI_Group * newgroup )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Group_range_excl()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Group_range_excl( group, n, ranges, newgroup );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Group_range_incl( 
MPI_Group group,
int n,
int ranges[][3],
MPI_Group * newgroup )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Group_range_incl()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Group_range_incl( group, n, ranges, newgroup );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Group_size( 
MPI_Group group,
int * size )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Group_size()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Group_size( group, size );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Group_translate_ranks( 
MPI_Group group_a,
int n,
int * ranks_a,
MPI_Group group_b,
int * ranks_b )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Group_translate_ranks()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Group_translate_ranks( group_a, n, ranks_a, group_b, ranks_b );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Group_union( 
MPI_Group group1,
MPI_Group group2,
MPI_Group * group_out )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Group_union()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Group_union( group1, group2, group_out );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Intercomm_create( 
MPI_Comm local_comm,
int local_leader,
MPI_Comm peer_comm,
int remote_leader,
int tag,
MPI_Comm * comm_out )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Intercomm_create()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Intercomm_create( local_comm, local_leader, peer_comm, remote_leader, tag, comm_out );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Intercomm_merge( 
MPI_Comm comm,
int high,
MPI_Comm * comm_out )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Intercomm_merge()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Intercomm_merge( comm, high, comm_out );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Keyval_create( 
MPI_Copy_function * copy_fn,
MPI_Delete_function * delete_fn,
int * keyval,
void * extra_state )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Keyval_create()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Keyval_create( copy_fn, delete_fn, keyval, extra_state );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Keyval_free( 
int * keyval )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Keyval_free()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Keyval_free( keyval );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Abort( 
MPI_Comm comm,
int errorcode )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Abort()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Abort( comm, errorcode );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Error_class( 
int errorcode,
int * errorclass )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Error_class()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Error_class( errorcode, errorclass );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Errhandler_create( 
MPI_Handler_function * function,
MPI_Errhandler * errhandler )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Errhandler_create()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Errhandler_create( function, errhandler );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Errhandler_free( 
MPI_Errhandler * errhandler )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Errhandler_free()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Errhandler_free( errhandler );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Errhandler_get( 
MPI_Comm comm,
MPI_Errhandler * errhandler )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Errhandler_get()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Errhandler_get( comm, errhandler );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Error_string( 
int errorcode,
char * string,
int * resultlen )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Error_string()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Error_string( errorcode, string, resultlen );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Errhandler_set( 
MPI_Comm comm,
MPI_Errhandler errhandler )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Errhandler_set()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Errhandler_set( comm, errhandler );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Finalize(  )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Finalize()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Finalize(  );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Get_processor_name( 
char * name,
int * resultlen )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Get_processor_name()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Get_processor_name( name, resultlen );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Init( 
int * argc,
char *** argv )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Init()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Init( argc, argv );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Initialized( 
int * flag )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Initialized()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Initialized( flag );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
double MPI_Wtick(  )
{
  double returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Wtick()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Wtick(  );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Address( 
void * location,
MPI_Aint * address )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Address()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Address( location, address );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Bsend( 
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Bsend()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Bsend( buf, count, datatype, dest, tag, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Bsend_init( 
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm,
MPI_Request * request )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Bsend_init()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Bsend_init( buf, count, datatype, dest, tag, comm, request );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Buffer_attach( 
void * buffer,
int size )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Buffer_attach()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Buffer_attach( buffer, size );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Buffer_detach( 
void * buffer,
int * size )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Buffer_detach()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Buffer_detach( buffer, size );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Cancel( 
MPI_Request * request )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Cancel()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Cancel( request );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Request_free( 
MPI_Request * request )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Request_free()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Request_free( request );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Recv_init( 
void * buf,
int count,
MPI_Datatype datatype,
int source,
int tag,
MPI_Comm comm,
MPI_Request * request )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Recv_init()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Recv_init( buf, count, datatype, source, tag, comm, request );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Send_init( 
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm,
MPI_Request * request )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Send_init()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Send_init( buf, count, datatype, dest, tag, comm, request );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Get_elements( 
MPI_Status * status,
MPI_Datatype datatype,
int * elements )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Get_elements()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Get_elements( status, datatype, elements );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Get_count( 
MPI_Status * status,
MPI_Datatype datatype,
int * count )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Get_count()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Get_count( status, datatype, count );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Ibsend( 
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm,
MPI_Request * request )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Ibsend()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Ibsend( buf, count, datatype, dest, tag, comm, request );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Iprobe( 
int source,
int tag,
MPI_Comm comm,
int * flag,
MPI_Status * status )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Iprobe()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Iprobe( source, tag, comm, flag, status );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Irecv( 
void * buf,
int count,
MPI_Datatype datatype,
int source,
int tag,
MPI_Comm comm,
MPI_Request * request )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Irecv()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Irecv( buf, count, datatype, source, tag, comm, request );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Irsend( 
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm,
MPI_Request * request )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Irsend()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Irsend( buf, count, datatype, dest, tag, comm, request );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Isend( 
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm,
MPI_Request * request )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Isend()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Isend( buf, count, datatype, dest, tag, comm, request );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Issend( 
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm,
MPI_Request * request )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Issend()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Issend( buf, count, datatype, dest, tag, comm, request );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Pack( 
void * inbuf,
int incount,
MPI_Datatype type,
void * outbuf,
int outcount,
int * position,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Pack()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Pack( inbuf, incount, type, outbuf, outcount, position, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Pack_size( 
int incount,
MPI_Datatype datatype,
MPI_Comm comm,
int * size )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Pack_size()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Pack_size( incount, datatype, comm, size );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Probe( 
int source,
int tag,
MPI_Comm comm,
MPI_Status * status )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Probe()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Probe( source, tag, comm, status );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Recv( 
void * buf,
int count,
MPI_Datatype datatype,
int source,
int tag,
MPI_Comm comm,
MPI_Status * status )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Recv()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Recv( buf, count, datatype, source, tag, comm, status );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Rsend( 
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Rsend()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Rsend( buf, count, datatype, dest, tag, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Rsend_init( 
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm,
MPI_Request * request )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Rsend_init()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Rsend_init( buf, count, datatype, dest, tag, comm, request );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Send( 
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Send()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Send( buf, count, datatype, dest, tag, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Sendrecv( 
void * sendbuf,
int sendcount,
MPI_Datatype sendtype,
int dest,
int sendtag,
void * recvbuf,
int recvcount,
MPI_Datatype recvtype,
int source,
int recvtag,
MPI_Comm comm,
MPI_Status * status )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Sendrecv()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Sendrecv( sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Sendrecv_replace( 
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int sendtag,
int source,
int recvtag,
MPI_Comm comm,
MPI_Status * status )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Sendrecv_replace()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Sendrecv_replace( buf, count, datatype, dest, sendtag, source, recvtag, comm, status );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Ssend( 
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Ssend()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Ssend( buf, count, datatype, dest, tag, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Ssend_init( 
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm,
MPI_Request * request )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Ssend_init()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Ssend_init( buf, count, datatype, dest, tag, comm, request );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Start( 
MPI_Request * request )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Start()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Start( request );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Startall( 
int count,
MPI_Request * array_of_requests )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Startall()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Startall( count, array_of_requests );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Test( 
MPI_Request * request,
int * flag,
MPI_Status * status )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Test()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Test( request, flag, status );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Testall( 
int count,
MPI_Request * array_of_requests,
int * flag,
MPI_Status * array_of_statuses )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Testall()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Testall( count, array_of_requests, flag, array_of_statuses );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Testany( 
int count,
MPI_Request * array_of_requests,
int * index,
int * flag,
MPI_Status * status )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Testany()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Testany( count, array_of_requests, index, flag, status );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Test_cancelled( 
MPI_Status * status,
int * flag )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Test_cancelled()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Test_cancelled( status, flag );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Testsome( 
int incount,
MPI_Request * array_of_requests,
int * outcount,
int * array_of_indices,
MPI_Status * array_of_statuses )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Testsome()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Testsome( incount, array_of_requests, outcount, array_of_indices, array_of_statuses );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Type_commit( 
MPI_Datatype * datatype )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Type_commit()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Type_commit( datatype );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Type_contiguous( 
int count,
MPI_Datatype old_type,
MPI_Datatype * newtype )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Type_contiguous()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Type_contiguous( count, old_type, newtype );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Type_extent( 
MPI_Datatype datatype,
MPI_Aint * extent )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Type_extent()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Type_extent( datatype, extent );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Type_free( 
MPI_Datatype * datatype )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Type_free()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Type_free( datatype );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Type_hindexed( 
int count,
int * blocklens,
MPI_Aint * indices,
MPI_Datatype old_type,
MPI_Datatype * newtype )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Type_hindexed()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Type_hindexed( count, blocklens, indices, old_type, newtype );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Type_hvector( 
int count,
int blocklen,
MPI_Aint stride,
MPI_Datatype old_type,
MPI_Datatype * newtype )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Type_hvector()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Type_hvector( count, blocklen, stride, old_type, newtype );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Type_indexed( 
int count,
int * blocklens,
int * indices,
MPI_Datatype old_type,
MPI_Datatype * newtype )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Type_indexed()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Type_indexed( count, blocklens, indices, old_type, newtype );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Type_lb( 
MPI_Datatype datatype,
MPI_Aint * displacement )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Type_lb()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Type_lb( datatype, displacement );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Type_size( 
MPI_Datatype datatype,
int * size )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Type_size()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Type_size( datatype, size );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Type_struct( 
int count,
int * blocklens,
MPI_Aint * indices,
MPI_Datatype * old_types,
MPI_Datatype * newtype )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Type_struct()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Type_struct( count, blocklens, indices, old_types, newtype );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Type_ub( 
MPI_Datatype datatype,
MPI_Aint * displacement )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Type_ub()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Type_ub( datatype, displacement );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Type_vector( 
int count,
int blocklen,
int stride,
MPI_Datatype old_type,
MPI_Datatype * newtype )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Type_vector()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Type_vector( count, blocklen, stride, old_type, newtype );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Unpack( 
void * inbuf,
int insize,
int * position,
void * outbuf,
int outcount,
MPI_Datatype type,
MPI_Comm comm )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Unpack()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Unpack( inbuf, insize, position, outbuf, outcount, type, comm );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Wait( 
MPI_Request * request,
MPI_Status * status )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Wait()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Wait( request, status );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Waitall( 
int count,
MPI_Request * array_of_requests,
MPI_Status * array_of_statuses )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Waitall()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Waitall( count, array_of_requests, array_of_statuses );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Waitany( 
int count,
MPI_Request * array_of_requests,
int * index,
MPI_Status * status )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Waitany()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Waitany( count, array_of_requests, index, status );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Waitsome( 
int incount,
MPI_Request * array_of_requests,
int * outcount,
int * array_of_indices,
MPI_Status * array_of_statuses )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Waitsome()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Waitsome( incount, array_of_requests, outcount, array_of_indices, array_of_statuses );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Cart_coords( 
MPI_Comm comm,
int rank,
int maxdims,
int * coords )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Cart_coords()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Cart_coords( comm, rank, maxdims, coords );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Cart_create( 
MPI_Comm comm_old,
int ndims,
int * dims,
int * periods,
int reorder,
MPI_Comm * comm_cart )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Cart_create()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Cart_create( comm_old, ndims, dims, periods, reorder, comm_cart );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Cart_get( 
MPI_Comm comm,
int maxdims,
int * dims,
int * periods,
int * coords )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Cart_get()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Cart_get( comm, maxdims, dims, periods, coords );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Cart_map( 
MPI_Comm comm_old,
int ndims,
int * dims,
int * periods,
int * newrank )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Cart_map()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Cart_map( comm_old, ndims, dims, periods, newrank );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Cart_rank( 
MPI_Comm comm,
int * coords,
int * rank )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Cart_rank()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Cart_rank( comm, coords, rank );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Cart_shift( 
MPI_Comm comm,
int direction,
int displ,
int * source,
int * dest )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Cart_shift()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Cart_shift( comm, direction, displ, source, dest );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Cart_sub( 
MPI_Comm comm,
int * remain_dims,
MPI_Comm * comm_new )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Cart_sub()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Cart_sub( comm, remain_dims, comm_new );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Cartdim_get( 
MPI_Comm comm,
int * ndims )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Cartdim_get()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Cartdim_get( comm, ndims );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Dims_create( 
int nnodes,
int ndims,
int * dims )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Dims_create()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Dims_create( nnodes, ndims, dims );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Graph_create( 
MPI_Comm comm_old,
int nnodes,
int * index,
int * edges,
int reorder,
MPI_Comm * comm_graph )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Graph_create()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Graph_create( comm_old, nnodes, index, edges, reorder, comm_graph );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Graph_get( 
MPI_Comm comm,
int maxindex,
int maxedges,
int * index,
int * edges )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Graph_get()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Graph_get( comm, maxindex, maxedges, index, edges );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Graph_map( 
MPI_Comm comm_old,
int nnodes,
int * index,
int * edges,
int * newrank )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Graph_map()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Graph_map( comm_old, nnodes, index, edges, newrank );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Graph_neighbors( 
MPI_Comm comm,
int rank,
int maxneighbors,
int * neighbors )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Graph_neighbors()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Graph_neighbors( comm, rank, maxneighbors, neighbors );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Graph_neighbors_count( 
MPI_Comm comm,
int rank,
int * nneighbors )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Graph_neighbors_count()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Graph_neighbors_count( comm, rank, nneighbors );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Graphdims_get( 
MPI_Comm comm,
int * nnodes,
int * nedges )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Graphdims_get()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Graphdims_get( comm, nnodes, nedges );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}

#ifdef __cplusplus
extern "C"
#endif
int MPI_Topo_test( 
MPI_Comm comm,
int * top_type )
{
  int returnVal;

  BL_PROFILE_TIMER( bltimer, "MPI_Topo_test()" );	
  BL_PROFILE_START( bltimer );
  
  returnVal = PMPI_Topo_test( comm, top_type );

  BL_PROFILE_STOP( bltimer );

  return returnVal;
}
