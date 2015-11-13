
       subroutine ext_src(lo,hi, &
                          old_state,old_state_l1,old_state_l2,old_state_h1,old_state_h2,&
                          new_state,new_state_l1,new_state_l2,new_state_h1,new_state_h2,&
                          src,src_l1,src_l2,src_h1,src_h2,dx,time,dt)

       use meth_params_module, only : NVAR

       implicit none
       integer         , intent(in   ) :: lo(2),hi(2)
       integer         , intent(in   ) :: old_state_l1,old_state_l2,old_state_h1,old_state_h2
       integer         , intent(in   ) :: new_state_l1,new_state_l2,new_state_h1,new_state_h2
       integer         , intent(in   ) :: src_l1,src_l2,src_h1,src_h2
       double precision, intent(in   ) :: old_state(old_state_l1:old_state_h1,old_state_l2:old_state_h2,NVAR)
       double precision, intent(in   ) :: new_state(new_state_l1:new_state_h1,new_state_l2:new_state_h2,NVAR)
       double precision, intent(  out) :: src(    src_l1:  src_h1,  src_l2:src_h2  ,NVAR)
       double precision, intent(in   ) :: dx(2),time,dt

       integer          :: i,j

       src = 0.d0

       end subroutine ext_src
