add_sources ( AMReX_EB2_F.H           AMReX_EB2_IF_Box.H        AMReX_EB2_IF_Lathe.H )
add_sources ( AMReX_EB2_IF_Union.H    AMReX_EB2_GeometryShop.H  AMReX_EB2_IF_Complement.H )
add_sources ( AMReX_EB2_IF_Plane.H    AMReX_EB2_IndexSpaceI.H   AMReX_EB2_Graph.H )
add_sources ( AMReX_EB2_IF_Cylinder.H AMReX_EB2_IF_Scale.H      AMReX_EB2_Level.H )
add_sources ( AMReX_EB2.H             AMReX_EB2_IF_Ellipsoid.H  AMReX_EB2_IF_Sphere.H )
add_sources ( AMReX_EB2_MultiGFab.H   AMReX_EB2_IF_AllRegular.H AMReX_EB2_IF_Intersection.H )
add_sources ( AMReX_EB2_IF_Translation.H AMReX_EB2_IF_Rotation.H)

add_sources( AMReX_EB2.cpp  AMReX_EB2_Level.cpp  AMReX_EB2_MultiGFab.cpp )

add_sources ( AMReX_eb2_${DIM}d.F90 )
