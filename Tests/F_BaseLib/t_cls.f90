subroutine t_cluster

  use cluster_module
  use multifab_module

  type(lmultifab) tags
  type(lmultifab) ctags
  type(boxarray)  ba, rba
  type(layout)    la, cla
  type(box)       bx
  integer         ratio

  bx = make_box((/0,0,0/),(/7,7,7/))

  ratio = 2

  call build(ba,bx)
  call build(la, ba, bx)
  call destroy(ba)
  call build(tags,la,1,0)
  call setval(tags,.true.,all=.true.)
!  call setval(tags,.true.)

!  call print(tags, 'tags')

!  call tagboxes_coarsen(tags,ctags,ratio)

  call cluster(rba, tags, 2)

!  call print(ctags, 'ctags')
!  print*, ''

  call print(rba, "grids")

  call destroy(rba)

!  cla = ctags%la

!  call destroy(ctags)
!  call destroy(cla)


  call destroy(tags)
  call destroy(la)

end subroutine t_cluster
