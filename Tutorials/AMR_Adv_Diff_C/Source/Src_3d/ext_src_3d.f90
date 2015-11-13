
       subroutine ext_src(lo,hi, &
                          old_state,os_l1,os_l2,os_l3,os_h1,os_h2,os_h3,&
                          new_state,ns_l1,ns_l2,ns_l3,ns_h1,ns_h2,ns_h3,&
                          src      ,sr_l1,sr_l2,sr_l3,sr_h1,sr_h2,sr_h3,dx,time,dt)

       use meth_params_module, only : NVAR

       implicit none
       integer         , intent(in   ) :: lo(3),hi(3)
       integer         , intent(in   ) :: os_l1,os_l2,os_l3,os_h1,os_h2,os_h3
       integer         , intent(in   ) :: ns_l1,ns_l2,ns_l3,ns_h1,ns_h2,ns_h3
       integer         , intent(in   ) :: sr_l1,sr_l2,sr_l3,sr_h1,sr_h2,sr_h3
       double precision, intent(in   ) :: old_state(os_l1:os_h1,os_l2:os_h2,os_l3:os_h3,NVAR)
       double precision, intent(in   ) :: new_state(ns_l1:ns_h1,ns_l2:ns_h2,ns_l3:ns_h3,NVAR)
       double precision, intent(  out) :: src      (sr_l1:sr_h1,sr_l2:sr_h2,sr_l3:sr_h3,NVAR)
       double precision, intent(in   ) :: dx(3),time,dt

       integer          :: i,j,k

       src = 0.d0

       end subroutine ext_src
