.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran


.. _sec:EB:EBOverview:

Overview of Embedded Boundary Description
=========================================

For computations with complex geometries, AMReX provides data structures and
algorithms to employ an embedded boundary (EB) approach to PDE discretizations.
In this approach, the underlying computational mesh is uniform and
block-structured, but the boundary of the irregular-shaped computational domain
conceptually cuts through this mesh. Each cell in the mesh becomes labeled as
regular, cut or covered, and the finite-volume based discretization methods
traditionally used in AMReX applications can be modified to incorporate these
cell shapes. See :numref:`fig::ebexample` for an illustration.

.. raw:: latex

   \centering

.. _fig::ebexample:

.. figure:: ./EB/EB_example.png
   :width: 50.0%

   : In the embedded boundary approach to discretizing PDEs, the (uniform)
   rectangular mesh is cut by the irregular shape of the computational domain.
   The cells in the mesh are label as regular, cut or covered.

Because this is a relatively simple grid generation technique, computational
meshes for rather complex geometries can be generated quickly and robustly.
However, the technique can produce arbitrarily small cut cells in the domain.
In practice such small cells can have significant impact on the robustness and
stability of traditional finite volume methods. In this chapter we overview a
class of approaches to deal with this “small cell” problem in a robust and
efficient way, and discuss the tools and data that AMReX provides in order to
implement them.

Note that in a completely general implementation of the EB approach, there
would be no restrictions on the shape or complexity of the EB surface.  With
this generality comes the possibility that the process of "cutting" the cells
results in a single :math:`(i,j,k)` cell being broken into multiple cell
fragments.  The current release of AMReX does not support multi-valued cells,
thus there is a practical restriction on the complexity of domains (and
numerical algorithms) supported.  AMReX support for EB with AMR will be
available by early 2018; EB support for multi-valued cells will follow.

This chapter discusses the EB tools, data structures and algorithms currently
supported by AMReX to enable the construction of discretizations of
conservation law systems. The discussion will focus on general requirements
associated with building fluxes and taking divergences of them to advance such
systems. We also give examples of how to initialize the geometry data
structures and access them to build the numerical difference operators.

Finite Volume Discretizations
-----------------------------

Consider a system of PDEs to advance a conserved quantity :math:`U` with fluxes
:math:`F`:

.. math:: \frac{\partial U}{\partial t} + \nabla \cdot F = 0.
  :label: eqn::hypsys

A conservative, finite volume discretization starts with the divergence theorm

.. math:: \int_V \nabla \cdot F dV = \int_{\partial V} F \cdot n dA.

In an embedded boundary cell, the “conservative divergence” is discretized (as
:math:`D^c(F)`) as follows

.. math::
  :label: eqn::ebdiv

   D^c(F) = \frac{1}{\kappa h} \left( \sum^D_{d = 1}
     (F_{d, \mathrm{hi}} \, A_{d, \mathrm{hi}} - F_{d, \mathrm{lo}}\, A_{d, \mathrm{lo}})
     + F^{EB} A^{EB} \right).

Geometry is discretely represented by volumes (:math:`V = \kappa h^d`) and
apertures (:math:`A= \alpha h^{d-1}`), where :math:`h` is the (uniform) mesh
spacing at that AMR level, :math:`\kappa` is the volume fraction and
:math:`\alpha` are the area fractions.  Without multivalued cells the volume
fractions, area fractions and cell and face centroids (see
:numref:`fig::volume`) are the only geometric information needed to compute
second-order fluxes centered at the face centroids, and to infer the
connectivity of the cells.  Cells are connected if adjacent on the Cartesian
mesh, and only via coordinate-aligned faces on the mesh. If an aperture,
:math:`\alpha = 0`, between two cells, they are not directly connected to each
other.

.. raw:: latex

   \centering

.. |a| image:: ./EB/areas_and_volumes.png
       :width: 100%

.. |b| image:: ./EB/eb_fluxes.png
       :width: 100%

.. _fig::volume:

.. table:: Illustration of embedded boundary cutting a two-dimensional cell.
   :align: center

   +-----------------------------------------------------+------------------------------------------------------+
   |                        |a|                          |                        |b|                           |
   +-----------------------------------------------------+------------------------------------------------------+
   | | A typical two-dimensional uniform cell that is    | | Fluxes in a cut cell.                              |
   | | cut by the embedded boundary. The grey area       | |                                                    |
   | | represents the region excluded from the           | |                                                    |
   | | calculation. The portion of the cell faces        | |                                                    |
   | | faces (labelled with A) through which fluxes flow | |                                                    |
   | | are the "uncovered" regions of the full cell      | |                                                    |
   | | faces. The volume (labelled V) is the uncovered   | |                                                    |
   | | region of the interior.                           | |                                                    |
   +-----------------------------------------------------+------------------------------------------------------+



Small Cells And Stability
-------------------------

In the context of time-explicit advance methods for, say hyperbolic
conservation laws, a naive discretization in time of :eq:`eqn::hypsys` using
:eq:`eqn::ebdiv`,

.. math:: U^{n+1} = U^{n} - \delta t D^c(F)

would have a time step constraint :math:`\delta t \sim h \kappa^{1/D}/V_m`,
which goes to zero as the size of the smallest volume fraction :math:`\kappa`
in the calculation. Since EB volume fractions can be arbitrarily small, this is
an unacceptable constraint. One way to remedy this is to create
“non-conservative” approximation to the divergence :math:`D^{nc}`, which at a
cell :math:`{\bf i}`, can be formed as an average of the conservative
divergences in the neighborhood, :math:`N_{\bf i}`, of :math:`{\bf i}`.

.. math:: D^{nc}(F)_{\bf i}= \frac{\sum_{{\bf j}\in N_{\bf i}}\kappa_{\bf j}D(F)_{\bf j}}{\sum_{{\bf j}\in N_{\bf i}}\kappa_{\bf j}}

Incorporating this form, the solution can be updated using a *hybrid
divergence*, :math:`D^H(F) = \kappa D^c(F) + (1-\kappa)D^{nc}`:

.. math:: U^{n+1,*} = U^n - \delta t D^H(F)

However, we would like our finite-volume scheme to strictly conserve the field
quantities over the domain. To enforce this, we calculate :math:`\delta M`, the
mass gained or lost by not using :math:`D^c` directly,

.. math:: \delta M_{\bf i}= \kappa (1-\kappa)(D^c(F)_{\bf i}- D^{nc}(F)_{\bf i})

This “excess material” (mass, if :math:`U=\rho`) can be *redistributed* in a
time-explicit fashion to neighboring cells, :math:`{\bf j}\in N_{\bf i}`:

.. math:: \delta M_{\bf i}= \sum_{{\bf j}\in N_{\bf i}} \delta M_{{\bf j}, {\bf i}}.

in order to preserve strict conservation over :math:`N_{\bf i}`.

Note that the physics at hand may impact the optimal choice of precisely how
the excess mass is distributed in this fashion. We introduce a weighting for
redistribution, :math:`W`,

.. math::
  :label: eqn::massweight

   \delta M_{{\bf j}, {\bf i}} =  \frac{\delta M_{\bf i}\kappa_{\bf j}
     W_{\bf j}}{\sum_{{\bf k}\in N_{\bf i}} \kappa_{\bf k}W_{\bf k}}

For all :math:`{\bf j}\in N_{\bf i}`,

.. math::

   U^{n+1}_{\bf j}= U^{n+1,*}_{\bf j}+
    \frac{\delta M_{\bf i}
     W_{\bf j}}{\sum_{{\bf k}\in N_{\bf i}} \kappa_{\bf k}W_{\bf k}}.

Typically, the redistribution neighborhood for each cell is one that can be
reached via a monotonic path in each coordinate direction of unit length (see,
e.g., :numref:`fig::redistribution`)

.. raw:: latex

   \centering

.. _fig::redistribution:

.. figure:: ./EB/redist.png
   :width: 50.0%

   : Redistribution illustration. Excess mass due to using a hybrid divergence
   :math:`D^H` instead of the conservative divergence :math:`D^C` is
   distributed to neighbor cells.


.. _sec:EB:ebinit:

Initializing EBIndexSpace, the Geometric Database
=================================================

In AMReX the geometric information is stored in a distributed database class,
:cpp:`EBIndexSpace`, which must be initialized at the start of the calculation.
The procedure for this goes as follows:

-  Define function of position which describes the surface and use it define a
   :cpp:`GeometryShop` object (see § :ref:`sec:EB:geometryshop`) – specifically,
   the scalar value returned by this function takes on a negative value inside
   the fluid, a positive value in the body, and identically zero at the EB.

-  Construct an :cpp:`EBIndexSpace` with the :cpp:`GeometryShop` object. This
   will fill the underlying database of geometric information, specifically
   tailored to the actual meshes that will be used. Thus, the construction
   requires one to specify the actual mesh resolution that will be used in a
   calculation.

To facilitate the first step, AMReX defines a virtual class, an “implicit
function”, :cpp:`BaseIF`, which encapsulates this functionality.  An instance
of a :cpp:`BaseIF` object is required for the construction of a
:cpp:`GeometryShop` object.

.. highlight:: c++

::

        GeometryShop(const BaseIF& a_localGeom)

Although the user is free to define their own instance of this class, AMReX
provides a number of preconfigured useful ones. This are listed in the next
section.

Example: Spherical EB
---------------------

The spherical implicit function, :cpp:`SphereIF`, derives from :cpp:`BaseIF`,
and defines the function

.. math:: S({\bf x}) = x^2 + y^2 + z^2 - R^2,

In this case, the solution domain is defined as the interior of a sphere of
radius :math:`R`. If the sign of :math:`S` is reversed, the solution domain is
the exterior of the sphere. The following example illustrates how to use the
SphereIF class to define a GeometryShop object:

.. highlight:: c++

::


      int nx = 1024;
      Box domain(IntVect::Zero, (nx-1)*IntVect::Unit);
      Real dx = 1.0/nx;
      Real radius = 0.1;
      RealVect center = 0.5*RealVect::Unit;
      bool insideRegular = true;
      //this is the implicit function
      SphereIF sphere(radius, center, insideRegular);

      //this is worker object that creates geometric information given an IF
      GeometryShop workshop(sphere)

      //this is the global, distributed database being initialized
      EBIndexSpace*  ebis = AMReX_EBIS::instance();
      ebis->define(domain, RealVect::Zero, dx, workshop);

In this case, we construct an :math:`r=0.1` sphere, centered within a unit
cube. The mesh resolution is :math:`1024^3`.  The :cpp:`GeometryShop` object
based on this sphere is then used to construct the :cpp:`EBIndexSpace`, as
shown.

Other basic shapes:
~~~~~~~~~~~~~~~~~~~

-  Planes are made using the class :cpp:`PlaneIF` which given a normal
   :math:`{\bf n}` and a center :math:`{\bf c}` gives the implicit function

   .. math:: I({\bf x}) = \sum_{1<=d<=D} n_d (x_d - c_d).

   .. highlight:: c++

   ::

       RealVect normal;
       RealVect center;
       // ...fill in values for normal and centre...

       PlaneIF plane(normal, centre, true);
       GeometryShop workshop(plane)

       EBIndexSpace*  ebis = AMReX_EBIS::instance();
       ebis->define(domain, RealVect::Zero, dx, workshop);

-  Polynomials of any form can be made using the class
   :cpp:`PolynomialIF`. Here is an example that makes a parabola of
   the form :math:`I({\bf x}) = x - y^2 - z^2`.

   .. highlight:: c++

   ::


       Vector<PolyTerm> poly;
       PolyTerm mono;
       Real coef;
       IntVect powers;
       Real amplitude = 1;

       // y^2 term
       coef = amplitude;
       powers = IntVect::Zero;
       powers[1] = 2;

       mono.coef   = coef;
       mono.powers = powers;
       poly.push_back(mono);

       // z^2 term
       coef = amplitude;
       RealVect translation;

       for(int idir = 0; idir < SpaceDim; idir++)
       {
           int finesize = finest_domain.size()[idir];
           translation[idir] = 0.5*finesize*fine_dx;
       }
       translation[0] = 0;

       TransformIF implicit(mirror);
       implicit.translate(translation);
       impfunc.reset(implicit.newImplicitFunction());

       powers = IntVect::Zero;
       powers[2] = 2;
       mono.coef   = coef;
       mono.powers = powers;
       poly.push_back(mono);

       // x term
       coef = -1.0;
       powers = IntVect::Zero;
       powers[0] = 1;
       mono.coef   = coef;
       mono.powers = powers;

       poly.push_back(mono);

       PolynomialIF mirror(poly,false);
       GeometryShop workshop(mirror)
       EBIndexSpace*  ebis = AMReX_EBIS::instance();
       ebis->define(domain, RealVect::Zero, dx, workshop);

Implicit Function Transformation Tools
--------------------------------------

More complex domains can be constructed by composing these fundamental shapes.
AMReX contains the following classes to compose implicit functions:

-  :cpp:`TransformIF` allows for translations and rotations of an implicit function.

-  :cpp:`UnionIF` produces the union of two implicit functions.

-  :cpp:`IntersectionIF` produces the intersection of two implicit functions.

-  :cpp:`LatheIF` creates a 3D implicit function as the surface of
   revolution of a 2D implicit function.

Multi-sphere example
--------------------

The following example creates a geometry using multiple spheres:

.. highlight:: c++

::


    vector<Real>     radius(numSpheres);
    vector<RealVect> center(numSpheres);
    //...

    //create an implicit function for each sphere
    vector<BaseIF*>  spheres(numSpheres);

    for(int isphere = 0; isphere < numSpheres; isphere++)
    {
      // Create sphere at each origin and translate
      SphereIF sphereAtZero(radius[isphere], RealVect::Zero, false);
      TransformIF* movedSphere = new TransformIF(sphereAtZero);
      movedSphere->translate(center[isphere]);
      spheres[isphere] = static_cast<BaseIF*>(movedSphere);
    }
    // Create implicit function as intersection of spheres
    IntersectionIF impMultisphere(spheres);

    // Fluid will in the complement space outside the sphere
    ComplementIF sideImpMultisphere(impMultisphere, false);

    // Construct the geometryshop
    GeometryShop workshop(sideImpMultisphere)

Geometric example 2 – Surface of revolution
-------------------------------------------

Here is an example that creates a geometric construction using a surface of
revolution of a set of polygons. This particular example only makes sense in
three dimensions. With the right polygons, it creates the surface shown in
:numref:`fig::revolution`.

.. highlight:: c++

::


    /// define EBIndexSpace from the surface of revolution of a set of polygons
    void
    defineGeometry(const Real& fine_dx, const  Box& finest_domain, int max_grid_size)
    {
      amrex::Print() << "creating geometry from polygon surfaces of revolution" << endl;

      // These  the polygons that get built around the z axis
      Vector<Vector<RealVect> > polygons;
      //....fill the polygons any way you like//

      // Make the Vector of (convex) polygons (Vectors of points) into a union
      // of convex polygons, each made from the intersection of a set of half
      // planes/spaces - all represented by implicit functions.

      // A list of all the polygons as implicit functions
      Vector<BaseIF*> polytopes;
      polytopes.resize(0);
      int numPolys = polygons.size();
      // Process each polygon
      for (int p = 0; p < numPolys; p++)
      {
        // All the half planes/spaces used to make a polygon
        Vector<BaseIF*> planes;
        planes.resize(0);

        // Get the current polygon (as a Vector of points)
        const Vector<RealVect>& polygon = polygons[p];

        // Get the number of points in the polygon
        int numPts = polygon.size();

        // Process each pair of points
        for (int n = 0; n < numPts; n++)
        {
          // The normal and point is space used to specify each half plane/space
          RealVect normal(RealVect::Zero);
          RealVect point;

          // Set the normal remembering that the last point connects to the first
          // point.
          normal[0] = -(polygon[(n+1) % numPts][1] - polygon[n][1]);
          normal[1] =  (polygon[(n+1) % numPts][0] - polygon[n][0]);

          point = polygon[n];

          // Generate the appropriate half plane/space (as an implicit function)
          PlaneIF* plane;
          plane = new PlaneIF(normal,point,true);

          // Save the result
          planes.push_back(plane);
        }

        // Intersect all the half planes/spaces to create an implicit function
        // that represents the polygon
        IntersectionIF* polygonIF = new IntersectionIF(planes);

        polytopes.push_back(polygonIF);
      }

      //this makes the cross section the union of all the polygons (around
      //z-axis, recall)
      UnionIF crossSection(polytopes);

      // In 3D rotate about the z-axis
      LatheIF lathe(crossSection, false);

      //we are starting around the z axis so we need to translate
      //over to the center of the x-y plane

      RealVect translation;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        translation[idir] = 0.5*finest_domain.size()[idir]*fine_dx;
      }
      translation[2] = 0;
      TransformIF implicit(lathe);
      implicit.translate(translation);

      //create a workshop from translated surface of revolution
      GeometryShop gshop(implicit, false);
      //define the geometric database
      AMReX_EBIS::instance()->define(finest_domain, RealVect::Zero,
                                     fine_dx, gshop, max_grid_size);
    }

.. raw:: latex

   \centering

.. _fig::revolution:

.. figure:: ./EB/revolution.png
   :width: 50.0%

   : Zero surface of an implicit function made using a surface of revolution.

Geometric example 3 – A Sphere Inside a Parabola
------------------------------------------------

Here is an example that creates a geometry of a sphere contained within a
parabola. This code creates the surface shown in :numref:`fig::parabolasphere`.

.. highlight:: c++

::

    Vector<PolyTerm> poly;

    PolyTerm mono;
    Real coef;
    IntVect powers;
    Real amplitude = 1;

    // y^2 term
    coef = amplitude;
    powers = IntVect::Zero;
    powers[1] = 2;

    mono.coef   = coef;
    mono.powers = powers;

    poly.push_back(mono);

    // z^2 term
    coef = amplitude;
    powers = IntVect::Zero;
    powers[2] = 2;
    mono.coef   = coef;
    mono.powers = powers;
    poly.push_back(mono);

    // x term
    coef = -1.0;
    powers = IntVect::Zero;
    powers[0] = 1;
    mono.coef   = coef;
    mono.powers = powers;

    poly.push_back(mono);

    PolynomialIF mirror(poly,false);
    RealVect translation;

    for(int idir = 0; idir < SpaceDim; idir++)
    {
      int finesize = finest_domain.size()[idir];
      translation[idir] = 0.5*finesize*fine_dx;
    }
    RealVect center = translation;
    translation[0] = 0;

    TransformIF transform(mirror);
    transform.translate(translation);

    Real radius = 0.2*center[0];
    SphereIF sphere(radius, center, true);
    Vector<BaseIF*> funcs(2);
    funcs[0] = &transform;
    funcs[1] = &sphere;
    UnionIF implicit(funcs);
    impfunc.reset(implicit.newImplicitFunction());
    GeometryShop gshop(impfunc, false);
    //define the geometric database
    AMReX_EBIS::instance()->define(finest_domain, RealVect::Zero,
                                     fine_dx, gshop, max_grid_size);

.. raw:: latex

   \centering

.. _fig::parabolasphere:

.. figure:: ./EB/parabsphere.png
   :width: 50.0%

   : Zero surface of an implicit function made the above code.

EBFarrayBox
===========

The fundamental data structure for embedded boundary calculations is
:cpp:`EBFArrayBox`. :cpp:`EBFArrayBox` is an a :cpp:`FArrayBox` with two extra
data members.

-  :cpp:`EBFArrayBox::getEBISBox` returns an :cpp:`EBISBox`, a data structure
   that contains the geometric information of an :cpp:`EBIndexSpace` but
   restricted to a given box.

-  :cpp:`EBFArrayBox::getEBCellFlagFab` is a :cpp:`BaseFab<EBCellFlag>`, where
   :cpp:`EBCellFlag` is a class which is a class with tools that compactly
   specifies local cell connectivities on a box.

If one compiles with ``AMREX_USE_EB = TRUE``, the state data managed by the
:cpp:`Amr` class is automatically of type :cpp:`EBFArrayBox` (typically the
data is exposed explicitly as a :cpp:`MultiFab`, but the additional
functionality may be accessed through a C++ type cast. The :cpp:`EBCellFlagFab`
can be used down in Fortran, e.g., to choose locally whether EB-specific
operations and data are required for constructing discretizations. In the next
section, we show examples of this workflow.

EBFarrayBox Usage Example
-------------------------

In order to make these EB concepts more concrete, we discuss here sample code
that appears in the AMReX tutorial, Tutorial/EB/CNS. This code implements a
time-explicit second-order method of lines integrator for hyperbolic and
parabolic transport based on a gamma-law gas EOS and constant transport
properties. This example also demonstrates how to avoid the more
complex/expensive EB-related logic if the tile under consideration has no cut
cells.

.. highlight:: c++

::

    void
    CNS::compute_dSdt (const MultiFab& S, MultiFab& dSdt, Real dt,
                       EBFluxRegister* fr_as_crse, EBFluxRegister* fr_as_fine)
    {
        BL_PROFILE("CNS::compute_dSdt()");

        const Real* dx = geom.CellSize();
        const int ncomp = dSdt.nComp();

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
         {
            //fluxes for the advance
            std::array<FArrayBox,AMREX_SPACEDIM> flux;

            for (MFIter mfi(S, MFItInfo().EnableTiling(hydro_tile_size).SetDynamic(true));
                            mfi.isValid(); ++mfi)
            {
                //this tile is the subset of the box over which we are computing
                const Box& bx = mfi.tilebox();

                //because S was created with the EBFArrayBoxFactory, we can do this cast
                const EBFArrayBox& sfab
                    = dynamic_cast<EBFArrayBox const&>(S[mfi]);

                //here we are getting the collection of flags so we know
                //kind of grid this is and if it is an EB grid, we have
                //the connectivity info
                const EBCellFlagFab & flag = sfab.getEBCellFlagFab();

                if (flag.getType(bx) == FabType::covered)
                {
                  //this tile is covered so there are no meaningful data here
                    dSdt[mfi].setVal(0.0, bx, 0, ncomp);
                }
                else
                {
                  //create the flux holders for this tile
                  for (int idim=0; idim < AMREX_SPACEDIM; ++idim)
                  {
                    flux[idim].resize(amrex::surroundingNodes(bx,idim),ncomp);
                  }

                  if (flag.getType(amrex::grow(bx,1)) == FabType::regular)
                  {
                    //this tile has no cut cells so we can just proceed
                    //with a (cheaper) non-eb call

                    cns_compute_dudt(BL_TO_FORTRAN_BOX(bx),
                    BL_TO_FORTRAN_ANYD(dSdt[mfi]),
                    BL_TO_FORTRAN_ANYD(S[mfi]),
                    BL_TO_FORTRAN_ANYD(flux[0]),
                    BL_TO_FORTRAN_ANYD(flux[1]),
                    BL_TO_FORTRAN_ANYD(flux[2]),
                    dx, &dt);

                  }
                  else
                  {
                    //this tile has cut cells so we have to send into Fortran
                    //EBCellFlagFAB as well as lots of geometric
                    //information
                    //the areafrac and facecent objects are member data
                    //filled using EBISBox
                    cns_eb_compute_dudt(BL_TO_FORTRAN_BOX(bx),
                    BL_TO_FORTRAN_ANYD(dSdt[mfi]),
                    BL_TO_FORTRAN_ANYD(S[mfi]),
                    BL_TO_FORTRAN_ANYD(flux[0]),
                    BL_TO_FORTRAN_ANYD(flux[1]),
                    BL_TO_FORTRAN_ANYD(flux[2]),
                    BL_TO_FORTRAN_ANYD(flag),
                    BL_TO_FORTRAN_ANYD(volfrac[mfi]),
                    BL_TO_FORTRAN_ANYD(bndrycent[mfi]),
                    BL_TO_FORTRAN_ANYD(areafrac[0][mfi]),
                    BL_TO_FORTRAN_ANYD(areafrac[1][mfi]),
                    BL_TO_FORTRAN_ANYD(areafrac[2][mfi]),
                    BL_TO_FORTRAN_ANYD(facecent[0][mfi]),
                    BL_TO_FORTRAN_ANYD(facecent[1][mfi]),
                    BL_TO_FORTRAN_ANYD(facecent[2][mfi]),
                    dx, &dt);
                  }
                }
              }
            }

This is the main loop in the routine to advance the state. The state,
:cpp:`MultiFab S`, comes into this routine with grow cells properly filled, and
this routine features a :cpp:`MultiFab` iterator loop to step through this
data, tile-by-tile and compute :cpp:`dSdt`. Here, we see that the definition of
:cpp:`EBFarrayBox sfab` incorporates the aforementioned type cast, enabling
queries about the EB nature of the data. Of the two possiblities handled, the
“regular” type without cut cells has a much simpler interface. The EB version
takes all the same data, but additionally requires (dense) data to specify the
volume and face area fractions, centroid information, and the
:cpp:`EBCellFlagFab flag` structure that will be queried pointwise for the
local cell connectivity.

Fortran code Snippets
---------------------

Much of the code to compute these fluxes and their divergence in this example
is too detailed to step through in this context. There are however a few
salient features worth pointing out.

The data is cell-centered, even cut cells
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to simplify the construction second-order discretizations, we can base
all the numerical operations on the assumption that all cell-based data lives
at the center of the *full*  cell containing the cut cells. This means that when
we take a standard centered difference between cell data at :math:`(i,j,k)` and
:math:`(i+1,j,k)`, e.g., we get a gradient value that is second-order and
centered on the full face at :math:`i+1/2`, regardless of the aperature.

Many EB operations can be organized as post-processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Recall that a second-order finite-volume scheme requires that fluxes be
centered on the face *centroid*. This can be accomplished by post-processing
face-centered fluxes with a linear interpolation of adjacent face values. The
resulting centroid-based fluxes are second-order, and can be used to construct
the conservative divergence we seek. Note that this operation requires the
location of the face centroids, and increases the grow cell requirement of the
flux operators, as does the necessity to form the *hybrid divergence* operator
discussed above.

The :cpp:`flag` data
~~~~~~~~~~~~~~~~~~~~~~~~~~

AMReX provides functions that query the :cpp:`flag` data in order to infer the
local connectivity of cells. For example, the cell itself or its neighbors may
be covered or cut. If cut, the data is centered at the center of the full cell.
If covered, the data is invalid and should not be involved in the fluid
advance. An example of such a call is:

.. highlight:: fortran

::

       call get_neighbor_cells(cellflag(i,j,k),nbr)

Here, for the :cpp:`flag` at :math:`(i,j,k)` is used to fill a local
:math:`3^3` array of integers with the value :math:`1` if connected to
:math:`(i,j,k)`, and :math:`0` if not. Similar queries:

.. highlight:: fortran

::

       is_covered_cell(cellflag(i,j,k))
       is_single_valued_cell(cellflag(i,j,k)

can be used to gather additional detail.

Below, we show a partial listing of the :fortran:`cns_eb_compute_dudt` code,
specifically after the face-centered fluxes have been computed, and showing
part of the work necessary to interpolate them to face centroids (while
appropriately handling covered data).

.. highlight:: fortran

::

        do n = 1, ncomp

           !
           ! First, we compute conservative divergence on (lo-2,hi+2)
           !
           iwall = 0
           do       k = lo(3)-2, hi(3)+2
              do    j = lo(2)-2, hi(2)+2
                 do i = lo(1)-2, hi(1)+2
                    divc(i,j,k) = (fluxx(i,j,k,n)-fluxx(i+1,j,k,n))*dxinv(1) &
                         +        (fluxy(i,j,k,n)-fluxy(i,j+1,k,n))*dxinv(2) &
                         +        (fluxz(i,j,k,n)-fluxz(i,j,k+1,n))*dxinv(3)
                 end do

                 do i = lo(1)-2, hi(1)+2
                    if (is_covered_cell(cellflag(i,j,k))) then
                       divc(i,j,k) = 0.d0
                    else if (is_single_valued_cell(cellflag(i,j,k))) then

                       call get_neighbor_cells(cellflag(i,j,k),nbr)

                       ! x-direction lo face
                       if (apx(i,j,k).lt.1.d0) then
                          if (centx_y(i,j,k).le.0.d0) then
                             fracy = -centx_y(i,j,k)*nbr(0,-1,0)
                             if(centx_z(i,j,k).le. 0.0d0)then
                                fracz = - centx_z(i,j,k)*nbr(0,0,-1)
                                fxm = (1.d0-fracz)*(     fracy *fluxx(i,j-1,k  ,n)  + &
                                     &             (1.d0-fracy)*fluxx(i,j  ,k  ,n)) + &
                                     &      fracz *(     fracy *fluxx(i,j-1,k-1,n)  + &
                                     &             (1.d0-fracy)*fluxx(i,j  ,k-1,n))
                             else
                                fracz =  centx_z(i,j,k)*nbr(0,0,1)
                                fxm = (1.d0-fracz)*(     fracy *fluxx(i,j-1,k  ,n)  + &
                                     &             (1.d0-fracy)*fluxx(i,j  ,k  ,n)) + &
                                     &      fracz *(     fracy *fluxx(i,j-1,k+1,n)  + &
                                     &             (1.d0-fracy)*fluxx(i,j  ,k+1,n))
                             endif
                          else
                             fracy = centx_y(i,j,k)*nbr(0,1,0)
                             if(centx_z(i,j,k).le. 0.0d0)then
                                fracz = -centx_z(i,j,k)*nbr(0,0,-1)
                                fxm = (1.d0-fracz)*(     fracy *fluxx(i,j+1,k  ,n)  + &
                                     &             (1.d0-fracy)*fluxx(i,j  ,k  ,n)) + &
                                     &      fracz *(     fracy *fluxx(i,j+1,k-1,n)  + &
                                     &             (1.d0-fracy)*fluxx(i,j  ,k-1,n))
                             else
                                fracz = centx_z(i,j,k)*nbr(0,0,1)
                                fxm = (1.d0-fracz)*(     fracy *fluxx(i,j+1,k  ,n)  + &
                                     &             (1.d0-fracy)*fluxx(i,j  ,k  ,n)) + &
                                     &      fracz *(     fracy *fluxx(i,j+1,k+1,n)  + &
                                     &             (1.d0-fracy)*fluxx(i,j  ,k+1,n))
                             endif
                          end if
                       else
                          fxm = fluxx(i,j,k,n)
                       end if

               <..... similar code for other fluxes removed ....>

                       iwall = iwall + 1
                       if (n .eq. 1) then
                          call compute_hyp_wallflux(divhyp(:,iwall), i,j,k, q(i,j,k,qrho), &
                               q(i,j,k,qu), q(i,j,k,qv), q(i,j,k,qw), q(i,j,k,qp), &
                               apx(i,j,k), apx(i+1,j,k), &
                               apy(i,j,k), apy(i,j+1,k), &
                               apz(i,j,k), apz(i,j,k+1))
                          call compute_diff_wallflux(divdiff(:,iwall), dxinv, i,j,k, &
                               q, qlo, qhi, &
                               lam, mu, xi, clo, chi, &
                               bcent, blo, bhi, &
                               apx, axlo, axhi, &
                               apy, aylo, ayhi, &
                               apz, azlo, azhi)
                       end if

                       divwn = divhyp(n,iwall) + divdiff(n,iwall)

                       ! we assume dx == dy == dz
                       divc(i,j,k) = -((apx(i+1,j,k)*fxp - apx(i,j,k)*fxm) * dxinv(1) &
                            +          (apy(i,j+1,k)*fyp - apy(i,j,k)*fym) * dxinv(2) &
                            +          (apz(i,j,k+1)*fzp - apz(i,j,k)*fzm) * dxinv(3) &
                            +          divwn * dxinv(1)) / vfrac(i,j,k)
                    end if
                 end do
              end do
           end do

One can easily identify the logic and portions of the code devoted toward the
EB corrections. Note, in particular, that diffusive fluxes into the EB need
only be computed on cut cells.

There are many approaches
~~~~~~~~~~~~~~~~~~~~~~~~~

The “fixes” that need to occur in these EB algorithms can be managed in a
number of ways, depending on the needs of the application and programming
style. In this example, the geometrical data is used to fill dense data
structures so that the sparse geometry information is available uniformally
over the entire box.  Also, the cell types are queried point-by-point in order
to form the appropriate stencil. Obviously then there is a performance penalty
if many of the cells in tile are not actually cut. There is clearly a trade-off
in such designs. Alternatively, one might build sparse data structures similar
to those AMReX uses to manage particles, and apply the EB corrections on this
sparse set directly. Future releases of AMReX will feature an expanded set of
EB tutorials to demonstrate an evolving set of tools provided.
