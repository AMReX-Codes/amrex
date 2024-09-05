.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _sec:EB:ebinit:

Initializing the Geometric Database
===================================

In AMReX geometric information is stored in a distributed database
class that must be initialized at the start of the calculation. The
procedure for this goes as follows:

- Define an implicit function of position which describes the surface of the
  embedded object. Specifically, the function class must have a public member
  function that takes a position and returns a negative value if that position
  is inside the fluid, a positive value in the body, and identically zero at the
  embedded boundary.

.. highlight:: c++

::

   Real operator() (const Array<Real,AMREX_SPACEDIM>& p) const;

- Make a :cpp:`EB2::GeometryShop` object using the implicit function.

- Build an :cpp:`EB2::IndexSpace` with the :cpp:`EB2::GeometryShop` object and a
  :cpp:`Geometry` object that contains the information about the domain and the
  mesh.

Here is a simple example of initialize the database for an embedded sphere.

.. highlight:: c++

::

    Real radius = 0.5;
    Array<Real,AMREX_SPACEDIM> center{0., 0., 0.}; //Center of the sphere
    bool inside = false;  // Is the fluid inside the sphere?
    EB2::SphereIF sphere(radius, center, inside);

    auto shop = EB2::makeShop(sphere);

    Geometry geom(...);
    EB2::Build(shop, geom, 0, 0);

Alternatively, the EB information can be initialized from an STL file
specified by a :cpp:`ParmParse` parameter ``eb2.stl_file``. (This also requires setting ``eb2.geom_type = stl``.) The
initialization is done by calling

.. highlight:: c++

::

   EB2::Build (const Geometry& geom,
               int required_coarsening_level,
               int max_coarsening_level,
               int ngrow = 4,
               bool build_coarse_level_by_coarsening = true);

Additionally one can use ``eb2.stl_scale``, ``eb2.stl_center`` and
``eb2.stl_reverse_normal`` to scale, translate and reverse the object,
respectively.

.. _sec:EB:ebinit:IF:

Implicit Function
-----------------

In ``amrex/Src/EB/``, there are a number of predefined implicit function classes
for basic shapes. One can use these directly or as template for their own
classes.

- :cpp:`AllRegularIF`:  No embedded boundaries at all.

- :cpp:`BoxIF`: Box.

- :cpp:`CylinderIF`: Cylinder.

- :cpp:`EllipsoidIF`: Ellipsoid.

- :cpp:`PlaneIF`: Half-space plane.

- :cpp:`SphereIF`: Sphere.

AMReX also provides a number of transformation operations to apply to an object.

- :cpp:`makeComplement`: Complement of an object. E.g. a sphere with fluid on
  outside becomes a sphere with fluid inside.

- :cpp:`makeIntersection`: Intersection of two or more objects.

- :cpp:`makeUnion`: Union of two or more objects.

- :cpp:`Translate`: Translates an object.

- :cpp:`scale`: Scales an object.

- :cpp:`rotate`: Rotates an object.

- :cpp:`lathe`: Creates a surface of revolution by rotating a 2D object around an axis.

Here are some examples of using these functions.

.. highlight: c++

::

    EB2::SphereIF sphere1(...);
    EB2::SphereIF sphere2(...);
    EB2::BoxIF box(...);
    EB2::CylinderIF cylinder(...);
    EB2::PlaneIF plane(...);

    // union of two spheres
    auto twospheres = EB2::makeUnion(sphere1, sphere2);

    // intersection of a rotated box, a plane and the union of two spheres
    auto box_plane = EB2::makeIntersection(amrex::rotate(box,...),
                                           plane,
                                           twospheres);

    // scale a cylinder by a factor of 2 in x and y directions, and 3 in z-direction.
    auto scylinder = EB2::scale(cylinder, {2., 2., 3.});

:cpp:`EB2::GeometryShop`
------------------------

Given an implicit function object, say :cpp:`f`, we can make a
:cpp:`GeometryShop` object with

.. highlight: c++

::

    auto shop = EB2::makeShop(f);

:cpp:`EB2::IndexSpace`
----------------------

We build :cpp:`EB2::IndexSpace` with a template function

.. highlight: c++

::

    template <typename G>
    void EB2::Build (const G& gshop, const Geometry& geom,
                     int required_coarsening_level,
                     int max_coarsening_level,
                     int ngrow = 4);

Here the template parameter is a :cpp:`EB2::GeometryShop`. :cpp:`Geometry` (see
section :ref:`sec:basics:geom`) describes the rectangular problem domain and the
mesh on the finest AMR level. Coarse level EB data is generated from coarsening
the original fine data. The :cpp:`int required_coarsening_level` parameter
specifies the number of coarsening levels required. This is usually set to
:math:`N-1`, where :math:`N` is the total number of AMR levels. The :cpp:`int
max_coarsening_levels` parameter specifies the number of coarsening levels AMReX
should try to have. This is usually set to a big number, say 20 if multigrid
solvers are used. This essentially tells the build to coarsen as much as it can.
If there are no multigrid solvers, the parameter should be set to the same as
:cpp:`required_coarsening_level`. It should be noted that coarsening could
create multi-valued cells even if the fine level does not have any multi-valued
cells. This occurs when the embedded boundary cuts a cell in such a way that
there is fluid on multiple sides of the boundary within that cell. Because
multi-valued cells are not supported, it will cause a runtime error if the
required coarsening level generates multi-valued cells. The optional :cpp:`int
ngrow` parameter specifies the number of ghost cells outside the domain on
required levels. For levels coarser than the required level, no EB data are
generated for ghost cells outside the domain.

The newly built :cpp:`EB2::IndexSpace` is pushed on to a stack. Static function
:cpp:`EB2::IndexSpace::top()` returns a :cpp:`const &` to the new
:cpp:`EB2::IndexSpace` object. We usually only need to build one
:cpp:`EB2::IndexSpace` object. However, if your application needs multiple
:cpp:`EB2::IndexSpace` objects, you can save the pointers for later use. For
simplicity, we assume there is only one `EB2::IndexSpace` object for the rest of
this chapter.

EBFArrayBoxFactory
==================

After the EB database is initialized, the next thing we build is
:cpp:`EBFArrayBoxFactory`. This object provides access to the EB database in the
format of basic AMReX objects such as :cpp:`BaseFab`, :cpp:`FArrayBox`,
:cpp:`FabArray`, and :cpp:`MultiFab`. We can construct it with

.. highlight: c++

::

    EBFArrayBoxFactory (const Geometry& a_geom,
                        const BoxArray& a_ba,
                        const DistributionMapping& a_dm,
                        const Vector<int>& a_ngrow,
                        EBSupport a_support);

or

.. highlight: c++

::

    std::unique_ptr<EBFArrayBoxFactory>
    makeEBFabFactory (const Geometry& a_geom,
                      const BoxArray& a_ba,
                      const DistributionMapping& a_dm,
                      const Vector<int>& a_ngrow,
                      EBSupport a_support);

Argument :cpp:`Vector<int> const& a_ngrow` specifies the number of
ghost cells we need for EB data at various :cpp:`EBSupport` levels,
and argument :cpp:`EBSupport a_support` specifies the level of support
needed.

- :cpp:`EBSupport:basic`:  basic flags for cell types
- :cpp:`EBSupport:volume`: basic plus volume fraction and centroid
- :cpp:`EBSupport:full`: volume plus area fraction, boundary centroid
  and face centroid

:cpp:`EBFArrayBoxFactory` is derived from :cpp:`FabFactory<FArrayBox>`.
:cpp:`MultiFab` constructors have an optional argument :cpp:`const
FabFactory<FArrayBox>&`.  We can use :cpp:`EBFArrayBoxFactory` to
build :cpp:`MultiFab`\ s that carry EB data.  Member function of
:cpp:`FabArray`

.. highlight: c++

::

    const FabFactory<FAB>& Factory () const;

can then be used to return a reference to the :cpp:`EBFArrayBoxFactory` used for
building the :cpp:`MultiFab`. Using :cpp:`dynamic_cast`, we can test whether a
:cpp:`MultiFab` is built with an :cpp:`EBFArrayBoxFactory`.

.. highlight: c++

::

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(&(mf.Factory()));
    if (factory) {
        // this is EBFArrayBoxFactory
    } else {
        // regular FabFactory<FArrayBox>
    }

Embedded Boundary Data
======================

Through member functions of :cpp:`EBFArrayBoxFactory`, we have access to the
following data:

.. highlight: c++

::

    // see section on EBCellFlagFab
    const FabArray<EBCellFlagFab>& getMultiEBCellFlagFab () const;

    // volume fraction
    const MultiFab& getVolFrac () const;

    // volume centroid
    const MultiCutFab& getCentroid () const;

    // embedded boundary centroid
    const MultiCutFab& getBndryCent () const;

    // embedded boundary normal direction
    const MultiCutFab& getBndryNormal () const;

    // embedded boundary surface area
    const MultiCutFab& getBndryArea () const;

    // area fractions
    Array<const MultiCutFab*,AMREX_SPACEDIM> getAreaFrac () const;

    // face centroid
    Array<const MultiCutFab*,AMREX_SPACEDIM> getFaceCent () const;

- **Volume fraction** is in a single-component :cpp:`MultiFab`. Data are in the range
  of :math:`[0,1]` with zero representing covered cells and one for regular
  cells.

- **Volume centroid** (also called cell centroid) is
  in a :cpp:`MultiCutFab` with ``AMREX_SPACEDIM`` components. Each component
  of the data is in the range of :math:`[-0.5,0.5]`, based on each
  cell's local coordinates with respect to the regular cell's center.

- **Boundary centroid** is also in a :cpp:`MultiCutFab` with
  ``AMREX_SPACEDIM`` components.  Each component
  of the data is in the range of :math:`[-0.5,0.5]`, based on each
  cell's local coordinates with respect to the regular cell's center.

- **Boundary normal** is in a :cpp:`MultiCutFab` with ``AMREX_SPACEDIM``
  components representing the unit vector pointing toward the covered part.

- **Boundary area** is in a :cpp:`MultiCutFab` with a single component
  representing the dimensionless boundary area. When the cell is isotropic
  (i.e., :math:`\Delta x = \Delta y = \Delta z`), it's trivial to convert it
  to physical units. If the cell size is anisotropic, the conversion
  requires multiplying by a factor of :math:`\sqrt{(n_x \Delta y \Delta
  z)^2 + (n_y \Delta x \Delta z)^2 + (n_z \Delta x \Delta y)^2}`, where
  :math:`n` is the boundary normal vector.

- **Face centroid** is in a :cpp:`MultiCutFab` with ``AMREX_SPACEDIM`` components.
  Each component of the data is in the range of :math:`[-0.5,0.5]`, based on
  each cell's local coordinates with respect to the embedded boundary.

- **Area fractions** are returned in an :cpp:`Array` of :cpp:`MultiCutFab`
  pointers. For each direction, area fraction is for the face of that direction.
  Data are in the range of :math:`[0,1]` with zero representing a covered face
  and one an un-cut face.

- **Face centroids** are returned in an :cpp:`Array` of :cpp:`MultiCutFab`
  pointers. There are two components for each direction and the
  ordering is always the same as the original ordering of the coordinates. For
  example, for :math:`y` face, the component 0 is for :math:`x` coordinate and 1
  for :math:`z`. The coordinates are in each face's local frame normalized to the
  range of :math:`[-0.5,0.5]`.


Embedded Boundary Data Structures
=================================

A :cpp:`MultiCutFab` is very similar to a :cpp:`MultiFab`. Its data can be
accessed with subscript operator

.. highlight: c++

::

    const CutFab& operator[] (const MFIter& mfi) const;

Here :cpp:`CutFab` is derived from :cpp:`FArrayBox` and can be passed to Fortran
just like :cpp:`FArrayBox`. The difference between :cpp:`MultiCutFab` and
:cpp:`MultiFab` is that to save memory :cpp:`MultiCutFab` only has data on boxes
that contain cut cells. It is an error to call :cpp:`operator[]` if that box
does not have cut cells. Thus the call must be in a :cpp:`if` test block (see
section :ref:`sec:EB:flag`).

.. _sec:EB:flag:

:cpp:`EBCellFlagFab`
--------------------

:cpp:`EBCellFlagFab` contains information on cell types.  We can use
it to determine if a box contains cut cells.

.. highlight: c++

::

    auto const& flags = factory->getMultiEBCellFlagFab();
    MultiCutFab const& centroid = factory->getCentroid();

    for (MFIter mfi ...) {
        const Box& bx = mfi.tilebox();
        FabType t = flags[mfi].getType(bx);
        if (FabType::regular == t) {
            // This box is regular
        } else if (FabType::covered == t) {
            // This box is covered
        } else if (FabType::singlevalued == t) {
            // This box has cut cells
            // Getting cutfab is safe
            const auto& centroid_fab = centroid[mfi];
        }
    }

:cpp:`EBCellFlagFab` is derived from :cpp:`BaseFab`. Its data are stored in an
array of 32-bit integers, and can be used in C++ or passed to Fortran just like
an :cpp:`IArrayBox` (section :ref:`sec:basics:fab`). AMReX provides a Fortran
module called ``amrex_ebcellflag_module``. This module contains procedures for
testing cell types and getting neighbor information. For example

.. highlight:: fortran

::

    use amrex_ebcellflag_module, only : is_regular_cell, is_single_valued_cell, is_covered_cell

    integer, intent(in) :: flags(...)

    integer :: i,j,k

    do k = ...
        do j = ...
            do i = ...
                if (is_covered_cell(flags(i,j,k))) then
                    ! this is a completely covered cells
                else if (is_regular_cell(flags(i,j,k))) then
                    ! this is a regular cell
                else if (is_single_valued_cell(flags(i,j,k))) then
                    ! this is a cut cell
                end if
            end do
        end do
    end do


.. _sec:EB:redistribution:

Small Cell Problem and Redistribution
=====================================

First, we review finite volume discretizations with embedded boundaries as used by
AMReX-based applications. Then we illustrate the small cell problem.

Finite Volume Discretizations
-----------------------------

Consider a system of PDEs to advance a conserved quantity :math:`U` with fluxes
:math:`F`:

.. math:: \frac{\partial U}{\partial t} + \nabla \cdot F = 0.
  :label: eqn::hypsys

A conservative, finite volume discretization starts with the divergence theorm

.. math:: \int_V \nabla \cdot F dV = \int_{\partial V} F \cdot n dA.

In an embedded boundary cell, the "conservative divergence" is discretized (as
:math:`D^c(F)`) as follows

.. math::
  :label: eqn::ebdiv

   D^c(F) = \frac{1}{\kappa h} \left( \sum^D_{d = 1}
     (F_{d, \mathrm{hi}} \, \alpha_{d, \mathrm{hi}} - F_{d, \mathrm{lo}}\, \alpha_{d, \mathrm{lo}})
     + F^{EB} \alpha^{EB} \right).

Geometry is discretely represented by volumes (:math:`V = \kappa h^d`) and
apertures (:math:`A= \alpha h^{d-1}`), where :math:`h` is the (uniform) mesh
spacing at that AMR level, :math:`\kappa` is the volume fraction and
:math:`\alpha` are the area fractions. Without multivalued cells the volume
fractions, area fractions and cell and face centroids (see
:numref:`fig::volume`) are the only geometric information needed to compute
second-order fluxes centered at the face centroids, and to infer the
connectivity of the cells. Cells are connected if adjacent on the Cartesian
mesh, and only via coordinate-aligned faces on the mesh. If an aperture,
:math:`\alpha = 0`, between two cells, they are not directly connected to each
other.

.. raw:: latex

   \begin{center}

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
   | | faces (labelled with A) through which fluxes      | |                                                    |
   | | flow are the "uncovered" regions of the full      | |                                                    |
   | | cell faces. The volume (labelled V) is the        | |                                                    |
   | | uncovered region of the interior.                 | |                                                    |
   +-----------------------------------------------------+------------------------------------------------------+

.. raw:: latex

   \end{center}


Small Cells And Stability
-------------------------

In the context of time-explicit advance methods for, say hyperbolic
conservation laws, a naive discretization in time of :eq:`eqn::hypsys` using
:eq:`eqn::ebdiv`,

.. math:: U^{n+1} = U^{n} - \delta t D^c(F)

would have a time step constraint :math:`\delta t \sim h \kappa^{1/D}/V_m`,
which goes to zero as the size of the smallest volume fraction :math:`\kappa` in
the calculation. Since EB volume fractions can be arbitrarily small, this presents an
unacceptable constraint. This is the so-called "small cell problem," and AMReX-based
applications address it with redistribution methods.

Flux Redistribution
-----------------------------

Consider a conservative update in the form:

.. math:: (\rho \phi)_t + \nabla \cdot ( \rho \phi u) = RHS

For each valid cell in the domain, compute the conservative divergence, :math:`(\nabla \cdot F)^c` ,
of the convective fluxes, :math:`F`

.. math:: (\nabla \cdot {F})^c_i = \dfrac{1}{\mathcal{V}_i} \sum_{f=1}^{N_f} ({F}_f\cdot{n}_f) A_f

Here :math:`N_f` is the number of faces of cell :math:`i`, :math:`\vec{n}_f` and :math:`A_f`
are the unit normal and area of the :math:`f` -th face respectively,
and :math:`\mathcal{V}_i` is the volume of cell :math:`i` given by

.. math:: \mathcal{V}_i = (\Delta x \Delta y \Delta z)\cdot \mathcal{K}_i

where :math:`\mathcal{K}_i` is the volume fraction of cell :math:`i` .

Now, a conservative update can be written as

.. math:: \frac{ \rho^{n+1} \phi^{n+1} - \rho^{n} \phi^{n} }{\Delta t} = - \nabla \cdot{F}^c

For each cell cut by the EB geometry, compute the non-conservative update, :math:`\nabla \cdot {F}^{nc}` ,

.. math:: \nabla\cdot{F}^{nc}_i = \dfrac{\sum\limits_{j\in N(i) } \mathcal{K}_j\nabla \cdot {F}^c_j} {\sum\limits_{j\in N(i) } {\mathcal{K}}_j}

where :math:`N(i)` is the index set of cell :math:`i` and its neighbors.

For each cell cut by the EB geometry, compute the convective update :math:`\nabla \cdot{F}^{EB}` follows:

.. math:: \nabla \cdot{F}^{EB}_i = \mathcal{K}_i\nabla \cdot{F}^{c}_i +(1-\mathcal{K}_i) \nabla \cdot \mathcal{F}^{nc}_i

For each cell cut by the EB geometry, redistribute its mass loss, :math:`\delta M_i` , to its neighbors:

.. math::  \nabla \cdot {F}^{EB}_j :=   \nabla \cdot {F}^{EB}_j + w_{ij}\delta M_i\, \qquad \forall j\in N(i)\setminus i

where the mass loss in cell :math:`i` , :math:`\delta M_i` , is given by

.. math:: \delta M_i =  \mathcal{K}_i(1- \mathcal{K}_i)[ \nabla \cdot {F}^c_i-  \nabla \cdot {F}^{nc}_i]

and the weights, :math:`w_{ij}` , are

.. math:: w_{ij} = \dfrac{1}{\sum\limits_{j\in N(i)\setminus i}  \mathcal{K}_j}

Note that :math:`\nabla \cdot{F}_i^{EB}` gives an update for :math:`\rho \phi` ; i.e.,

.. math:: \frac{(\rho \phi_i)^{n+1} - (\rho \phi_i)^{n} }{\Delta t} = - \nabla \cdot{F}^{EB}_i

Typically, the redistribution neighborhood for each cell is one that can be
reached via a monotonic path in each coordinate direction of unit length (see,
e.g., :numref:`fig::redistribution`)

.. raw:: latex

   \begin{center}

.. _fig::redistribution:

.. figure:: ./EB/redist.png
   :width: 50.0%

   : Redistribution illustration. Excess update distributed to neighbor cells.

.. raw:: latex

   \end{center}


State Redistribution
-----------------------------

For state redistribution we implement the weighted state
redistribution algorithm as described in Guiliani et al (2021),
which is available on `arxiv  <https://arxiv.org/abs/2112.12360>`_ .
This is an extension of the original state redistribution algorithm
of Berger and Guiliani (2020).


Linear Solvers
==============

Linear solvers for the canonical form (equation :eq:`eqn::abeclap`)
have been discussed in chapter :ref:`Chap:LinearSolvers`.

AMReX supports multi-level
1) cell-centered solvers with homogeneous Neumann, homogeneous Dirichlet,
or inhomogeneous Dirichlet boundary conditions on the EB faces, and
2) nodal solvers with homogeneous Neumann boundary conditions,
or inflow velocity conditions on the EB faces.

To use a cell-centered solver with EB, one builds a linear operator
:cpp:`MLEBABecLap` with :cpp:`EBFArrayBoxFactory` (instead of a :cpp:`MLABecLaplacian`)

.. highlight:: c++

::

    MLEBABecLap (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory);

The usage of this EB-specific class is essentially the same as
:cpp:`MLABecLaplacian`.

The default boundary condition on EB faces is homogeneous Neumann.

To set homogeneous Dirichlet boundary conditions, call

.. highlight:: c++

::

    ml_ebabeclap->setEBHomogDirichlet(lev, coeff);

where coeff can be a real number (i.e. the value is the same at every cell)
or is the MultiFab holding the coefficient of the gradient at each cell with an EB face.

To set inhomogeneous Dirichlet boundary conditions, call

.. highlight:: c++

::

    ml_ebabeclap->setEBDirichlet(lev, phi_on_eb, coeff);

where phi_on_eb is the MultiFab holding the Dirichlet values in every cut cell,
and coeff again is a real number (i.e. the value is the same at every cell)
or a MultiFab holding the coefficient of the gradient at each cell with an EB face.

Currently there are options to define the face-based coefficients on
face centers vs face centroids, and to interpret the solution variable
as being defined on cell centers vs cell centroids.

The default is for the solution variable to be defined at cell centers;
to tell the solver to interpret the solution variable as living
at cell centroids, you must set

.. highlight:: c++

::

    ml_ebabeclap->setPhiOnCentroid();

The default is for the face-based coefficients to be defined at face centers;
to tell the that the face-based coefficients should be interpreted
as living at face centroids, modify the setBCoeffs command to be

.. highlight:: c++

::

    ml_ebabeclap->setBCoeffs(lev, beta, MLMG::Location::FaceCentroid);

Tutorials
=========

`EB/CNS`_ is an AMR code for solving compressible
Navier-Stokes equations with the embedded boundary approach.

`EB/Poisson`_ is a single-level code that is a proxy for
solving the electrostatic Poisson equation for a grounded sphere with a point
charge inside.

`EB/MacProj`_ is a single-level code that computes a divergence-free
flow field around a sphere.  A MAC projection is performed on an initial velocity
field of (1,0,0).

.. _`EB/CNS`: https://amrex-codes.github.io/amrex/tutorials_html/EB_Tutorial.html

.. _`EB/Poisson`: https://amrex-codes.github.io/amrex/tutorials_html/EB_Tutorial.html

.. _`EB/MacProj`: https://amrex-codes.github.io/amrex/tutorials_html/EB_Tutorial.html
