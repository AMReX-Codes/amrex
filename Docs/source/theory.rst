========
Theory
========

Introduction
============

Computer simulations have had a profound impact on the design and
understanding of past and present plasma acceleration experiments
(**???**; **???**; **???**; **???**), with accurate modeling of wake
formation, electron self-trapping and acceleration requiring fully
kinetic methods (usually Particle-In-Cell) using large computational
resources due to the wide range of space and time scales involved.
Numerical modeling complements and guides the design and analysis of
advanced accelerators, and can reduce development costs significantly.
Despite the major recent experimental successesLeemans et al. (2014;
Blumenfeld et al. 2007; Bulanov S V and Wilkens J J and Esirkepov T Zh
and Korn G and Kraft G and Kraft S D and Molls M and Khoroshkov V S
2014; Steinke et al. 2016), the various advanced acceleration concepts
need significant progress to fulfill their potential. To this end,
large-scale simulations will continue to be a key component toward
reaching a detailed understanding of the complex interrelated physics
phenomena at play.

For such simulations, the most popular algorithm is the Particle-In-Cell
(or PIC) technique, which represents electromagnetic fields on a grid
and particles by a sample of macroparticles. However, these simulations
are extremely computationally intensive, due to the need to resolve the
evolution of a driver (laser or particle beam) and an accelerated beam
into a structure that is orders of magnitude longer and wider than the
accelerated beam. Various techniques or reduced models have been
developed to allow multidimensional simulations at manageable
computational costs: quasistatic approximation (**???**; **???**;
**???**; Mora and Antonsen 1997; Huang et al. 2006), ponderomotive
guiding center (PGC) models (**???**; **???**; Huang et al. 2006;
**???**; **???**), simulation in an optimal Lorentz boosted frame
(**???**; Bruhwiler et al. 2009; Vay et al. 2009; Vay
:math:`\backslash`\ it Et Al. 2009; Martins :math:`\backslash`\ It Et
Al. 2009; Vay et al. 2010; **???**; **???**; **???**; **???**; Vay et
al. 2011; **???**; Yu et al. 2016), expanding the fields into a
truncated series of azimuthal modes (**???**; Lifschitz et al. 2009;
Davidson et al. 2015; Lehe et al. 2016; Andriyash, Lehe, and Lifschitz
2016), fluid approximation (**???**; **???**; **???**) and scaled
parameters (**???**; **???**). Many codes have been developed and are
used for the modeling of plasma accelerators. A list of such codes is
given in table [table\_codes], with the name of the code, its main
characteristics, the web site if existing or a reference, and the
availability and license, if known.

[]

| @llll@ Code & Type & Web site/reference & Availability/License
| ALaDyn/PICCANTE & EM-PIC 3D & http://aladyn.github.io/piccante &
  Open/GPLv3+
| Architect & EM-PIC RZ & https://github.com/albz/Architect & Open/GPL
| Calder & EM-PIC 3D &
  http://iopscience.iop.org/article/10.1088/0029-5515/43/7/317 &
  Collaborators/Proprietary
| Calder-Circ & EM-PIC RZ\ :math:`^{+}` &
  http://dx.doi.org/10.1016/j.jcp.2008.11.017 & Upon Request/Proprietary
| CHIMERA & EM-PIC RZ+ & https://github.com/hightower8083/chimera &
  Open/GPLv3
| ELMIS & EM-PIC 3D &
  http://www.diva-portal.org/smash/record.jsf?pid=diva2%3A681092&dswid=-8610
  & Collaborators/Proprietary
| EPOCH & EM-PIC 3D& http://www.ccpp.ac.uk/codes.html &
  Collaborators/GPL
| FBPIC & EM-PIC RZ\ :math:`^{+}` & https://fbpic.github.io &
  Open/modified BSD
| HiPACE & QS-PIC 3D & http://dx.doi.org/10.1088/0741-3335/56/8/084012 &
  Collaborators/Proprietary
| INF&RNO & QS/EM-PIC RZ & http://dx.doi.org/10.1063/1.3520323 &
  Collaborators/Proprietary
| LCODE & QS-PIC RZ & http://www.inp.nsk.su/~lotov/lcode & Open/None
| LSP & EM-PIC 3D/RZ & http://www.lspsuite.com/LSP/index.html &
  Commercial/Proprietary
| MAGIC & EM-PIC 3D & http://www.mrcwdc.com/magic/index.html &
  Commercial/Proprietary
| Osiris & EM-PIC 3D/RZ\ :math:`^{+}` &
  http://picksc.idre.ucla.edu/software/software-production-codes/osiris
  & Collaborators/Proprietary
| PHOTON-PLASMA & EM-PIC 3D & https://bitbucket.org/thaugboelle/ppcode &
  Open/GPLv2
| PICADOR & EM-PIC 3D &
  http://hpc-education.unn.ru/en/research/overview/laser-plasma &
  Collaborators/Proprietary
| PIConGPU & EM-PIC 3D & http://picongpu.hzdr.de & Open/GPLv3+
| PICLS & EM-PIC 3D & http://dx.doi.org/10.1016/j.jcp.2008.03.043 &
  Collaborators/Proprietary
| PSC & EM-PIC 3D &
  http://www.sciencedirect.com/science/article/pii/S0021999116301413 &
  Open/GPLv3
| QuickPIC & QS-PIC 3D &
  http://picksc.idre.ucla.edu/software/software-production-codes/quickpic
  & Collaborators/Proprietary
| REMP & EM-PIC 3D & http://dx.doi.org/10.1016/S0010-4655(00)00228-9 &
  Collaborators/Proprietary
| Smilei & EM-PIC 2D &
  http://www.maisondelasimulation.fr/projects/Smilei/html/licence.html &
  Open/CeCILL
| TurboWave & EM-PIC 3D/RZ & http://dx.doi.org/10.1109/27.893300 &
  Collaborators/Proprietary
| UPIC-EMMA & EM-PIC 3D &
  http://picksc.idre.ucla.edu/software/software-production-codes/upic-emma
  & Collaborators/Proprietary
| VLPL & EM/QS-PIC 3D & http://www.tp1.hhu.de/~pukhov/ &
  Collaborators/Proprietary
| VPIC & EM-PIC 3D & http://github.com/losalamos/vpic & Open/BSD
  clause-3 license
| VSim (Vorpal) & EM-PIC 3D & https://txcorp.com/vsim &
  Commercial/Proprietary
| Wake & QS-PIC RZ & http://dx.doi.org/10.1063/1.872134 &
  Collaborators/Proprietary
| Warp & EM-PIC 3D/RZ\ :math:`^{+}` & http://warp.lbl.gov &
  Open/modified BSD

| EM=electromagnetic, QS=quasi-static, PIC=Particle-In-Cell,
  3D=three-dimensional, RZ=axi-symmetric, RZ\ :math:`^+`\ =axi-symmetric
  with azimuthal Fourier decomposition.

[table\_codes]

In Section 2 of this chapter, we review the standard methods employed in
relativistic electromagnetic Particle-In-Cell (PIC) simulations of
plasma accelerators, including the core PIC loop steps (particle push,
fields update, current deposition from the particles to the grid and
fields gathering from the grid to the particles positions), the use of
moving window and Lorentz boosted frame, the numerical Cherenkov
instability and its mitigation. The electromagnetic quasistatic
approximation is presented in section 3, the ponderomotive guiding
center approximation in section 4, and azimuthal Fourier decomposition
in section 5. Additional considerations such as filtering and
inputs/outputs are discussed respectively in sections 6 and 7.

The electromagnetic Particle-In-Cell method
===========================================

In the electromagnetic Particle-In-Cell method (**???**), the
electromagnetic fields are solved on a grid, usually using Maxwell’s
equations

.. math::

   \begin{aligned}
   \frac{\mathbf{\partial B}}{\partial t} & = & -\nabla\times\mathbf{E}\label{Eq:Faraday-1}\\
   \frac{\mathbf{\partial E}}{\partial t} & = & \nabla\times\mathbf{B}-\mathbf{J}\label{Eq:Ampere-1}\\
   \nabla\cdot\mathbf{E} & = & \rho\label{Eq:Gauss-1}\\
   \nabla\cdot\mathbf{B} & = & 0\label{Eq:divb-1}\end{aligned}

given here in natural units (:math:`\epsilon_0=\mu_0=c=1`), where
:math:`t` is time, :math:`\mathbf{E}` and :math:`\mathbf{B}` are the
electric and magnetic field components, and :math:`\rho` and
:math:`\mathbf{J}` are the charge and current densities. The charged
particles are advanced in time using the Newton-Lorentz equations of
motion

.. math::

   \begin{aligned}
   \frac{d\mathbf{x}}{dt}= & \mathbf{v},\label{Eq:Lorentz_x-1}\\
   \frac{d\left(\gamma\mathbf{v}\right)}{dt}= & \frac{q}{m}\left(\mathbf{E}+\mathbf{v}\times\mathbf{B}\right),\label{Eq:Lorentz_v-1}\end{aligned}

where :math:`m`, :math:`q`, :math:`\mathbf{x}`, :math:`\mathbf{v}` and
:math:`\gamma=1/\sqrt{1-v^{2}}` are respectively the mass, charge,
position, velocity and relativistic factor of the particle given in
natural units (:math:`c=1`). The charge and current densities are
interpolated on the grid from the particles’ positions and velocities,
while the electric and magnetic field components are interpolated from
the grid to the particles’ positions for the velocity update.

Particle push
-------------

A centered finite-difference discretization of the Newton-Lorentz
equations of motion is given by

.. math::

   \begin{aligned}
   \frac{\mathbf{x}^{i+1}-\mathbf{x}^{i}}{\Delta t}= & \mathbf{v}^{i+1/2},\label{Eq:leapfrog_x}\\
   \frac{\gamma^{i+1/2}\mathbf{v}^{i+1/2}-\gamma^{i-1/2}\mathbf{v}^{i-1/2}}{\Delta t}= & \frac{q}{m}\left(\mathbf{E}^{i}+\mathbf{\bar{v}}^{i}\times\mathbf{B}^{i}\right).\label{Eq:leapfrog_v}\end{aligned}

In order to close the system, :math:`\bar{\mathbf{v}}^{i}` must be
expressed as a function of the other quantities. The two implementations
that have become the most popular are presented below.

Boris relativistic velocity rotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The solution proposed by Boris Boris (1970) is given by

.. math::

   \begin{aligned}
   \mathbf{\bar{v}}^{i}= & \frac{\gamma^{i+1/2}\mathbf{v}^{i+1/2}+\gamma^{i-1/2}\mathbf{v}^{i-1/2}}{2\bar{\gamma}^{i}}.\label{Eq:boris_v}\end{aligned}

 where :math:`\bar{\gamma}^{i}` is defined by
:math:`\bar{\gamma}^{i} \equiv (\gamma^{i+1/2}+\gamma^{i-1/2} )/2`.

The system ([Eq:leapfrog\_v],[Eq:boris\_v]) is solved very efficiently
following Boris’ method, where the electric field push is decoupled from
the magnetic push. Setting :math:`\mathbf{u}=\gamma\mathbf{v}`, the
velocity is updated using the following sequence:

.. math::

   \begin{aligned}
   \mathbf{u^{-}}= & \mathbf{u}^{i-1/2}+\left(q\Delta t/2m\right)\mathbf{E}^{i}\\
   \mathbf{u'}= & \mathbf{u}^{-}+\mathbf{u}^{-}\times\mathbf{t}\\
   \mathbf{u}^{+}= & \mathbf{u}^{-}+\mathbf{u'}\times2\mathbf{t}/(1+t^{2})\\
   \mathbf{u}^{i+1/2}= & \mathbf{u}^{+}+\left(q\Delta t/2m\right)\mathbf{E}^{i}\end{aligned}

where :math:`\mathbf{t}=\left(q\Delta
  t/2m\right)\mathbf{B}^{i}/\bar{\gamma}^{i}` and where
:math:`\bar{\gamma}^{i}` can be calculated as
:math:`\bar{\gamma}^{i}=\sqrt{1+(\mathbf{u}^-/c)^2}`.

The Boris implementation is second-order accurate, time-reversible and
fast. Its implementation is very widespread and used in the vast
majority of PIC codes.

Lorentz-invariant formulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It was shown in (**???**) that the Boris formulation is not Lorentz
invariant and can lead to significant errors in the treatment of
relativistic dynamics. A Lorentz invariant formulation is obtained by
considering the following velocity average

.. math::

   \begin{aligned}
   \mathbf{\bar{v}}^{i}= & \frac{\mathbf{v}^{i+1/2}+\mathbf{v}^{i-1/2}}{2},\label{Eq:new_v}\end{aligned}

 This gives a system that is solvable analytically (see (**???**) for a
detailed derivation), giving the following velocity update:

.. math::

   \begin{aligned}
   \mathbf{u^{*}}= & \mathbf{u}^{i-1/2}+\frac{q\Delta t}{m}\left(\mathbf{E}^{i}+\frac{\mathbf{v}^{i-1/2}}{2}\times\mathbf{B}^{i}\right),\label{pusher_gamma}\\
   \mathbf{u}^{i+1/2}= & \left[\mathbf{u^{*}}+\left(\mathbf{u^{*}}\cdot\mathbf{t}\right)\mathbf{t}+\mathbf{u^{*}}\times\mathbf{t}\right]/\left(1+t^{2}\right),\label{pusher_upr}\end{aligned}

where :math:`\mathbf{t}=\bm{\tau}/\gamma^{i+1/2}`,
:math:`\bm{\tau}=\left(q\Delta t/2m\right)\mathbf{B}^{i}`,
:math:`\gamma^{i+1/2}=\sqrt{\sigma+\sqrt{\sigma^{2}+\left(\tau^{2}+w^{2}\right)}}`,
:math:`w=\mathbf{u^{*}}\cdot\bm{\tau}`,
:math:`\sigma=\left(\gamma'^{2}-\tau^{2}\right)/2` and
:math:`\gamma'=\sqrt{1+(\mathbf{u}^{*}/c)^{2}}`. This Lorentz invariant
formulation is particularly well suited for the modeling of
ultra-relativistic charged particle beams, where the accurate account of
the cancellation of the self-generated electric and magnetic fields is
essential, as shown in (**???**).

Field solve
-----------

Various methods are available for solving Maxwell’s equations on a grid,
based on finite-differences, finite-volume, finite-element, spectral, or
other discretization techniques that apply most commonly on single
structured or unstructured meshes and less commonly on multiblock
multiresolution grid structures. In this chapter, we summarize the
widespread second order finite-difference time-domain (FDTD) algorithm,
its extension to non-standard finite-differences as well as the
pseudo-spectral analytical time-domain (PSATD) and pseudo-spectral
time-domain (PSTD) algorithms. Extension to multiresolution (or mesh
refinement) PIC is described in, e.g. Vay et al. (2012; Vay, Adam, and
Heron 2004).

Finite-Difference Time-Domain (FDTD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most popular algorithm for electromagnetic PIC codes is the
Finite-Difference Time-Domain (or FDTD) solver

.. math::

   \begin{aligned}
   D_{t}\mathbf{B} & = & -\nabla\times\mathbf{E}\label{Eq:Faraday-2}\\
   D_{t}\mathbf{E} & = & \nabla\times\mathbf{B}-\mathbf{J}\label{Eq:Ampere-2}\\
   \left[\nabla\cdot\mathbf{E}\right. & = & \left.\rho\right]\label{Eq:Gauss-2}\\
   \left[\nabla\cdot\mathbf{B}\right. & = & \left.0\right].\label{Eq:divb-2}\end{aligned}

The differential operator is defined as
:math:`\nabla=D_{x}\mathbf{\hat{x}}+D_{y}\mathbf{\hat{y}}+D_{z}\mathbf{\hat{z}}`
and the finite-difference operators in time and space are defined
respectively as
:math:` `\ :math:`D_{t}G|_{i,j,k}^{n}=\left(G|_{i,j,k}^{n+1/2}-G|_{i,j,k}^{n-1/2}\right)/\Delta t`\ :math:` `
and
:math:`D_{x}G|_{i,j,k}^{n}=\left(G|_{i+1/2,j,k}^{n}-G|_{i-1/2,j,k}^{n}\right)/\Delta x`,
where :math:`\Delta t` and :math:`\Delta x` are respectively the time
step and the grid cell size along :math:`x`, :math:`n` is the time index
and :math:`i`, :math:`j` and :math:`k` are the spatial indices along
:math:`x`, :math:`y` and :math:`z` respectively. The difference
operators along :math:`y` and :math:`z` are obtained by circular
permutation. The equations in brackets are given for completeness, as
they are often not actually solved, thanks to the usage of a so-called
charge conserving algorithm, as explained below. As shown in Figure
[fig:yee\_grid], the quantities are given on a staggered (or “Yee”) grid
Yee (1966), where the electric field components are located between
nodes and the magnetic field components are located in the center of the
cell faces. Knowing the current densities at half-integer steps, the
electric field components are updated alternately with the magnetic
field components at integer and half-integer steps respectively.

Non-Standard Finite-Difference Time-Domain (NSFDTD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In (**???**; **???**), Cole introduced an implementation of the
source-free Maxwell’s wave equations for narrow-band applications based
on non-standard finite-differences (NSFD). In (**???**), Karkkainen *et
al.* adapted it for wideband applications. At the Courant limit for the
time step and for a given set of parameters, the stencil proposed in
(**???**) has no numerical dispersion along the principal axes, provided
that the cell size is the same along each dimension (i.e. cubic cells in
3D). The “Cole-Karkkainnen” (or CK) solver uses the non-standard finite
difference formulation (based on extended stencils) of the
Maxwell-Ampere equation and can be implemented as follows (**???**):

.. math::

   \begin{aligned}
   D_{t}\mathbf{B} & = & -\nabla^{*}\times\mathbf{E}\label{Eq:Faraday}\\
   D_{t}\mathbf{E} & = & \nabla\times\mathbf{B}-\mathbf{J}\label{Eq:Ampere}\\
   \left[\nabla\cdot\mathbf{E}\right. & = & \left.\rho\right]\label{Eq:Gauss}\\
   \left[\nabla^{*}\cdot\mathbf{B}\right. & = & \left.0\right]\label{Eq:divb}\end{aligned}

Eq. [Eq:Gauss] and [Eq:divb] are not being solved explicitly but
verified via appropriate initial conditions and current deposition
procedure. The NSFD differential operators is given by
:math:`\nabla^{*}=D_{x}^{*}\mathbf{\hat{x}}+D_{y}^{*}\mathbf{\hat{y}}+D_{z}^{*}\mathbf{\hat{z}}`
where
:math:`D_{x}^{*}=\left(\alpha+\beta S_{x}^{1}+\xi S_{x}^{2}\right)D_{x}`
with
:math:`S_{x}^{1}G|_{i,j,k}^{n}=G|_{i,j+1,k}^{n}+G|_{i,j-1,k}^{n}+G|_{i,j,k+1}^{n}+G|_{i,j,k-1}^{n}`,
:math:`S_{x}^{2}G|_{i,j,k}^{n}=G|_{i,j+1,k+1}^{n}+G|_{i,j-1,k+1}^{n}+G|_{i,j+1,k-1}^{n}+G|_{i,j-1,k-1}^{n}`.
:math:`G` is a sample vector component, while :math:`\alpha`,
:math:`\beta` and :math:`\xi` are constant scalars satisfying
:math:`\alpha+4\beta+4\xi=1`. As with the FDTD algorithm, the quantities
with half-integer are located between the nodes (electric field
components) or in the center of the cell faces (magnetic field
components). The operators along :math:`y` and :math:`z`, i.e.
:math:`D_{y}`, :math:`D_{z}`, :math:`D_{y}^{*}`, :math:`D_{z}^{*}`,
:math:`S_{y}^{1}`, :math:`S_{z}^{1}`, :math:`S_{y}^{2}`, and
:math:`S_{z}^{2}`, are obtained by circular permutation of the indices.

Assuming cubic cells (:math:`\Delta x=\Delta y=\Delta z`), the
coefficients given in (**???**) (:math:`\alpha=7/12`, :math:`\beta=1/12`
and :math:`\xi=1/48`) allow for the Courant condition to be at
:math:`\Delta t=\Delta x`, which equates to having no numerical
dispersion along the principal axes. The algorithm reduces to the FDTD
algorithm with :math:`\alpha=1` and :math:`\beta=\xi=0`. An extension to
non-cubic cells is provided by Cowan, *et al.* in 3-D in Cowan et al.
(2013) and was given by Pukhov in 2-D in Pukhov (1999). An alternative
NSFDTD implementation that enables superluminous waves is also given by
Lehe *et al.* in Lehe et al. (2013).

As mentioned above, a key feature of the algorithms based on NSFDTD is
that some implementations (**???**; Cowan et al. 2013) enable the time
step :math:`\Delta t=\Delta x` along one or more axes and no numerical
dispersion along those axes. However, as shown in (**???**), an
instability develops at the Nyquist wavelength at (or very near) such a
timestep. It is also shown in the same paper that removing the Nyquist
component in all the source terms using a bilinear filter (see
description of the filter below) suppresses this instability.

Pseudo Spectral Analytical Time Domain (PSATD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Maxwell’s equations in Fourier space are given by

.. math::

   \begin{aligned}
   \frac{\partial{\mathbf{\tilde{E}}}}{\partial t} & = & i{\mathbf{k}}\times{\mathbf{\tilde{B}}}-{\mathbf{\tilde{J}}}\\
   \frac{\partial{\mathbf{\tilde{B}}}}{\partial t} & = & -i{\mathbf{k}}\times{\mathbf{\tilde{E}}}\\
   {}[i{\mathbf{k}}\cdot{\mathbf{\tilde{E}}}& = & \tilde{\rho}]\\
   {}[i{\mathbf{k}}\cdot{\mathbf{\tilde{B}}}& = & 0]\end{aligned}

where :math:`\tilde{a}` is the Fourier Transform of the quantity
:math:`a`. As with the real space formulation, provided that the
continuity equation
:math:`\partial\tilde{\rho}/\partial t+i{\mathbf{k}}\cdot{\mathbf{\tilde{J}}}=0`
is satisfied, then the last two equations will automatically be
satisfied at any time if satisfied initially and do not need to be
explicitly integrated.

Decomposing the electric field and current between longitudinal and
transverse components
:math:`{\mathbf{\tilde{E}}}={\mathbf{\tilde{E}}}_{L}+{\mathbf{\tilde{E}}}_{T}={\mathbf{\hat{k}}}({\mathbf{\hat{k}}}\cdot{\mathbf{\tilde{E}}})-{\mathbf{\hat{k}}}\times({\mathbf{\hat{k}}}\times{\mathbf{\tilde{E}}})`
and
:math:`{\mathbf{\tilde{J}}}={\mathbf{\tilde{J}}}_{L}+{\mathbf{\tilde{J}}}_{T}={\mathbf{\hat{k}}}({\mathbf{\hat{k}}}\cdot{\mathbf{\tilde{J}}})-{\mathbf{\hat{k}}}\times({\mathbf{\hat{k}}}\times{\mathbf{\tilde{J}}})`
gives

.. math::

   \begin{aligned}
   \frac{\partial{\mathbf{\tilde{E}}}_{T}}{\partial t} & = & i{\mathbf{k}}\times{\mathbf{\tilde{B}}}-\mathbf{\tilde{J}_{T}}\\
   \frac{\partial{\mathbf{\tilde{E}}}_{L}}{\partial t} & = & -\mathbf{\tilde{J}_{L}}\\
   \frac{\partial{\mathbf{\tilde{B}}}}{\partial t} & = & -i{\mathbf{k}}\times{\mathbf{\tilde{E}}}\end{aligned}

with :math:`{\mathbf{\hat{k}}}={\mathbf{k}}/k`.

If the sources are assumed to be constant over a time interval
:math:`\Delta t`, the system of equations is solvable analytically and
is given by (see (**???**) for the original formulation and Vay, Haber,
and Godfrey (2013) for a more detailed derivation):

[Eq:PSATD]

.. math::

   \begin{aligned}
   {\mathbf{\tilde{E}}}_{T}^{n+1} & = & C{\mathbf{\tilde{E}}}_{T}^{n}+iS{\mathbf{\hat{k}}}\times{\mathbf{\tilde{B}}}^{n}-\frac{S}{k}{\mathbf{\tilde{J}}}_{T}^{n+1/2}\label{Eq:PSATD_transverse_1}\\
   {\mathbf{\tilde{E}}}_{L}^{n+1} & = & {\mathbf{\tilde{E}}}_{L}^{n}-\Delta t{\mathbf{\tilde{J}}}_{L}^{n+1/2}\\
   {\mathbf{\tilde{B}}}^{n+1} & = & C{\mathbf{\tilde{B}}}^{n}-iS{\mathbf{\hat{k}}}\times{\mathbf{\tilde{E}}}^{n}\\
   &+&i\frac{1-C}{k}{\mathbf{\hat{k}}}\times{\mathbf{\tilde{J}}}^{n+1/2}\label{Eq:PSATD_transverse_2}\end{aligned}

with :math:`C=\cos\left(k\Delta t\right)` and
:math:`S=\sin\left(k\Delta t\right)`.

Combining the transverse and longitudinal components, gives

.. math::

   \begin{aligned}
   {\mathbf{\tilde{E}}}^{n+1} & = & C{\mathbf{\tilde{E}}}^{n}+iS{\mathbf{\hat{k}}}\times{\mathbf{\tilde{B}}}^{n}-\frac{S}{k}{\mathbf{\tilde{J}}}^{n+1/2}\\
    & + &(1-C){\mathbf{\hat{k}}}({\mathbf{\hat{k}}}\cdot{\mathbf{\tilde{E}}}^{n})\nonumber \\
    & + & {\mathbf{\hat{k}}}({\mathbf{\hat{k}}}\cdot{\mathbf{\tilde{J}}}^{n+1/2})\left(\frac{S}{k}-\Delta t\right),\label{Eq_PSATD_1}\\
   {\mathbf{\tilde{B}}}^{n+1} & = & C{\mathbf{\tilde{B}}}^{n}-iS{\mathbf{\hat{k}}}\times{\mathbf{\tilde{E}}}^{n}\\
   &+&i\frac{1-C}{k}{\mathbf{\hat{k}}}\times{\mathbf{\tilde{J}}}^{n+1/2}.\label{Eq_PSATD_2}\end{aligned}

For fields generated by the source terms without the self-consistent
dynamics of the charged particles, this algorithm is free of numerical
dispersion and is not subject to a Courant condition. Furthermore, this
solution is exact for any time step size subject to the assumption that
the current source is constant over that time step.

As shown in Vay, Haber, and Godfrey (2013), by expanding the
coefficients :math:`S_{h}` and :math:`C_{h}` in Taylor series and
keeping the leading terms, the PSATD formulation reduces to the perhaps
better known pseudo-spectral time-domain (PSTD) formulation Dawson
(1983; **???**):

.. math::

   \begin{aligned}
   {\mathbf{\tilde{E}}}^{n+1} & = & {\mathbf{\tilde{E}}}^{n}+i\Delta t{\mathbf{k}}\times{\mathbf{\tilde{B}}}^{n+1/2}-\Delta t{\mathbf{\tilde{J}}}^{n+1/2},\\
   {\mathbf{\tilde{B}}}^{n+3/2} & = & {\mathbf{\tilde{B}}}^{n+1/2}-i\Delta t{\mathbf{k}}\times{\mathbf{\tilde{E}}}^{n+1}.\end{aligned}

The dispersion relation of the PSTD solver is given by
:math:`\sin(\frac{\omega\Delta t}{2})=\frac{k\Delta t}{2}.` In contrast
to the PSATD solver, the PSTD solver is subject to numerical dispersion
for a finite time step and to a Courant condition that is given by
:math:`\Delta t\leq \frac{2}{\pi}\left(\frac{1}{\Delta x^{2}}+\frac{1}{\Delta y^{2}}+\frac{1}{\Delta x^{2}}\right)^{-1/2}.`

The PSATD and PSTD formulations that were just given apply to the field
components located at the nodes of the grid. As noted in Ohmura and
Okamura (2010), they can also be easily recast on a staggered Yee grid
by multiplication of the field components by the appropriate phase
factors to shift them from the collocated to the staggered locations.
The choice between a collocated and a staggered formulation is
application-dependent.

Spectral solvers used to be very popular in the years 1970s to early
1990s, before being replaced by finite-difference methods with the
advent of parallel supercomputers that favored local methods. However,
it was shown recently that standard domain decomposition with Fast
Fourier Transforms that are local to each subdomain could be used
effectively with PIC spectral methods Vay, Haber, and Godfrey (2013), at
the cost of truncation errors in the guard cells that could be
neglected. A detailed analysis of the effectiveness of the method with
exact evaluation of the magnitude of the effect of the truncation error
is given in Vincenti and Vay (2016) for stencils of arbitrary order
(up-to the infinite “spectral” order).

Current deposition
------------------

The current densities are deposited on the computational grid from the
particle position and velocities, employing splines of various orders
(**???**).

.. math::

   \begin{aligned}
   \rho & = & \frac{1}{\Delta x \Delta y \Delta z}\sum_nq_nS_n\\
   \mathbf{J} & = & \frac{1}{\Delta x \Delta y \Delta z}\sum_nq_n\mathbf{v_n}S_n\end{aligned}

In most applications, it is essential to prevent the accumulation of
errors resulting from the violation of the discretized Gauss’ Law. This
is accomplished by providing a method for depositing the current from
the particles to the grid that preserves the discretized Gauss’ Law, or
by providing a mechanism for “divergence cleaning” (**???**; **???**;
**???**; **???**; Munz et al. 2000). For the former, schemes that allow
a deposition of the current that is exact when combined with the Yee
solver is given in (**???**) for linear splines and in Esirkepov (2001)
for splines of arbitrary order.

The NSFDTD formulations given above and in Pukhov (1999; **???**; Cowan
et al. 2013; Lehe et al. 2013) apply to the Maxwell-Faraday equation,
while the discretized Maxwell-Ampere equation uses the FDTD formulation.
Consequently, the charge conserving algorithms developed for current
deposition (**???**; Esirkepov 2001) apply readily to those NSFDTD-based
formulations. More details concerning those implementations, including
the expressions for the numerical dispersion and Courant condition are
given in Pukhov (1999; **???**; Cowan et al. 2013; Lehe et al. 2013).

In the case of the pseudospectral solvers, the current deposition
algorithm generally does not satisfy the discretized continuity equation
in Fourier space
:math:`\tilde{\rho}^{n+1}=\tilde{\rho}^{n}-i\Delta t{\mathbf{k}}\cdot\mathbf{\tilde{J}}^{n+1/2}`.
In this case, a Boris correction (**???**) can be applied in :math:`k`
space in the form
:math:`{\mathbf{\tilde{E}}}_{c}^{n+1}={\mathbf{\tilde{E}}}^{n+1}-\left({\mathbf{k}}\cdot{\mathbf{\tilde{E}}}^{n+1}+i\tilde{\rho}^{n+1}\right){\mathbf{\hat{k}}}/k`,
where :math:`{\mathbf{\tilde{E}}}_{c}` is the corrected field.
Alternatively, a correction to the current can be applied (with some
similarity to the current deposition presented by Morse and Nielson in
their potential-based model in (**???**)) using
:math:`{\mathbf{\tilde{J}}}_{c}^{n+1/2}={\mathbf{\tilde{J}}}^{n+1/2}-\left[{\mathbf{k}}\cdot{\mathbf{\tilde{J}}}^{n+1/2}-i\left(\tilde{\rho}^{n+1}-\tilde{\rho}^{n}\right)/\Delta t\right]{\mathbf{\hat{k}}}/k`,
where :math:`{\mathbf{\tilde{J}}}_{c}` is the corrected current. In this
case, the transverse component of the current is left untouched while
the longitudinal component is effectively replaced by the one obtained
from integration of the continuity equation, ensuring that the corrected
current satisfies the continuity equation. The advantage of correcting
the current rather than the electric field is that it is more local and
thus more compatible with domain decomposition of the fields for
parallel computation J. L. Vay, Haber, and Godfrey (2013).

Alternatively, an exact current deposition can be written for the
pseudospectral solvers, following the geometrical interpretation of
existing methods in real space (**???**; **???**; Esirkepov 2001),
thereby averaging the currents of the paths following grid lines between
positions :math:`(x^n,y^n)` and :math:`(x^{n+1},y^{n+1})`, which is
given in 2D (extension to 3D follows readily) for :math:`k\neq0` by J.
L. Vay, Haber, and Godfrey (2013):

.. math::

   \begin{aligned}
   {\mathbf{\tilde{J}}}^{k\neq0}=\frac{i\mathbf{\tilde{D}}}{{\mathbf{k}}}
   \label{Eq_Jdep_1}\end{aligned}

 with

.. math::

   \begin{aligned}
   D_x   =  \frac{1}{2\Delta t}\sum_i q_i
     [\Gamma(x_i^{n+1},y_i^{n+1})-\Gamma(x_i^{n},y_i^{n+1}) \nonumber\\
   +\Gamma(x_i^{n+1},y_i^{n})-\Gamma(x_i^{n},y_i^{n})],\\
   D_y   =  \frac{1}{2\Delta t}\sum_i q_i
     [\Gamma(x_i^{n+1},y_i^{n+1})-\Gamma(x_i^{n+1},y_i^{n}) \nonumber \\
   +\Gamma(x_i^{n},y_i^{n+1})-\Gamma(x_i^{n},y_i^{n})],\end{aligned}

 where :math:`\Gamma` is the macro-particle form factor. The
contributions for :math:`k=0` are integrated directly in real space J.
L. Vay, Haber, and Godfrey (2013).

Field gather
------------

In general, the field is gathered from the mesh onto the macroparticles
using splines of the same order as for the current deposition
:math:`\mathbf{S}=\left(S_{x},S_{y},S_{z}\right)`. Three variations are
considered:

-  “momentum conserving”: fields are interpolated from the grid nodes to
   the macroparticles using
   :math:`\mathbf{S}=\left(S_{nx},S_{ny},S_{nz}\right)` for all field
   components (if the fields are known at staggered positions, they are
   first interpolated to the nodes on an auxiliary grid),

-  “energy conserving (or Galerkin)”: fields are interpolated from the
   staggered Yee grid to the macroparticles using
   :math:`\left(S_{nx-1},S_{ny},S_{nz}\right)` for :math:`E_{x}`,
   :math:`\left(S_{nx},S_{ny-1},S_{nz}\right)` for :math:`E_{y}`,
   :math:`\left(S_{nx},S_{ny},S_{nz-1}\right)` for :math:`E_{z}`,
   :math:`\left(S_{nx},S_{ny-1},S_{nz-1}\right)` for :math:`B_{x}`,
   :math:`\left(S_{nx-1},S_{ny},S_{nz-1}\right)` for :math:`B{}_{y}`
   and\ :math:`\left(S_{nx-1},S_{ny-1},S_{nz}\right)` for :math:`B_{z}`
   (if the fields are known at the nodes, they are first interpolated to
   the staggered positions on an auxiliary grid),

-  “uniform”: fields are interpolated directly form the Yee grid to the
   macroparticles using
   :math:`\mathbf{S}=\left(S_{nx},S_{ny},S_{nz}\right)` for all field
   components (if the fields are known at the nodes, they are first
   interpolated to the staggered positions on an auxiliary grid).

As shown in (**???**; Hockney and Eastwood 1988; Lewis 1972), the
momentum and energy conserving schemes conserve momentum and energy
respectively at the limit of infinitesimal time steps and generally
offer better conservation of the respective quantities for a finite time
step. The uniform scheme does not conserve momentum nor energy in the
sense defined for the others but is given for completeness, as it has
been shown to offer some interesting properties in the modeling of
relativistically drifting plasmas Godfrey and Vay (2013).

Moving window and optimal Lorentz boosted frame
-----------------------------------------------

The simulations of plasma accelerators from first principles are
extremely computationally intensive, due to the need to resolve the
evolution of a driver (laser or particle beam) and an accelerated
particle beam into a plasma structure that is orders of magnitude longer
and wider than the accelerated beam. As is customary in the modeling of
particle beam dynamics in standard particle accelerators, a moving
window is commonly used to follow the driver, the wake and the
accelerated beam. This results in huge savings, by avoiding the meshing
of the entire plasma that is orders of magnitude longer than the other
length scales of interest.

Even using a moving window, however, a full PIC simulation of a plasma
accelerator can be extraordinarily demanding computationally, as many
time steps are needed to resolve the crossing of the short driver beam
with the plasma column. As it turns out, choosing an optimal frame of
reference that travels close to the speed of light in the direction of
the laser or particle beam (as opposed to the usual choice of the
laboratory frame) enables speedups by orders of magnitude (**???**;
**???**). This is a result of the properties of Lorentz contraction and
dilation of space and time. In the frame of the laboratory, a very short
driver (laser or particle) beam propagates through a much longer plasma
column, necessitating millions to tens of millions of time steps for
parameters in the range of the BELLA or FACET-II experiments. As
sketched in Fig. [fig:PIC], in a frame moving with the driver beam in
the plasma at velocity :math:`v=\beta c` (where :math:`c` is the speed
of light in vacuum), the beam length is now elongated by
:math:`\approx(1+\beta)\gamma` while the plasma contracts by
:math:`\gamma` (where :math:`\gamma=1/\sqrt{1-\beta^2}` is the
relativistic factor associated with the frame velocity). The number of
time steps that is needed to simulate a “longer” beam through a
“shorter” plasma is now reduced by up to
:math:`\approx(1+\beta) \gamma^2` (a detailed derivation of the speedup
is given below).

The modeling of a plasma acceleration stage in a boosted frame involves
the fully electromagnetic modeling of a plasma propagating at near the
speed of light, for which Numerical Cerenkov (**???**; **???**) is a
potential issue, as explained in more details below. In addition, for a
frame of reference moving in the direction of the accelerated beam (or
equivalently the wake of the laser), waves emitted by the plasma in the
forward direction expand while the ones emitted in the backward
direction contract, following the properties of the Lorentz
transformation. If one had to resolve both forward and backward
propagating waves emitted from the plasma, there would be no gain in
selecting a frame different from the laboratory frame. However, the
physics of interest for a laser wakefield is the laser driving the wake,
the wake, and the accelerated beam. Backscatter is weak in the
short-pulse regime, and does not interact as strongly with the beam as
do the forward propagating waves which stay in phase for a long period.
It is thus often assumed that the backward propagating waves can be
neglected in the modeling of plasma accelerator stages. The accuracy of
this assumption has been demonstrated by comparison between explicit
codes which include both forward and backward waves and envelope or
quasistatic codes which neglect backward waves (**???**; **???**;
**???**).

Theoretical speedup dependency with the frame boost
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The derivation that is given here reproduces the one given in (**???**),
where the obtainable speedup is derived as an extension of the formula
that was derived earlier(**???**), taking in addition into account the
group velocity of the laser as it traverses the plasma.

Assuming that the simulation box is a fixed number of plasma periods
long, which implies the use (which is standard) of a moving window
following the wake and accelerated beam, the speedup is given by the
ratio of the time taken by the laser pulse and the plasma to cross each
other, divided by the shortest time scale of interest, that is the laser
period. To first order, the wake velocity :math:`v_w` is set by the 1D
group velocity of the laser driver, which in the linear (low intensity)
limit, is given by (**???**):

.. math:: v_w/c=\beta_w=\left(1-\frac{\omega_p^2}{\omega^2}\right)^{1/2}

where :math:`\omega_p=\sqrt{(n_e e^2)/(\epsilon_0 m_e)}` is the plasma
frequency, :math:`\omega=2\pi c/\lambda` is the laser frequency,
:math:`n_e` is the plasma density, :math:`\lambda` is the laser
wavelength in vacuum, :math:`\epsilon_0` is the permittivity of vacuum,
:math:`c` is the speed of light in vacuum, and :math:`e` and :math:`m_e`
are respectively the charge and mass of the electron.

In practice, the runs are typically stopped when the last electron beam
macro-particle exits the plasma, and a measure of the total time of the
simulation is then given by

.. math:: T=\frac{L+\eta \lambda_p}{v_w-v_p}

 where :math:`\lambda_p\approx 2\pi c/\omega_p` is the wake wavelength,
:math:`L` is the plasma length, :math:`v_w` and :math:`v_p=\beta_p c`
are respectively the velocity of the wake and of the plasma relative to
the frame of reference, and :math:`\eta` is an adjustable parameter for
taking into account the fraction of the wake which exited the plasma at
the end of the simulation. For a beam injected into the :math:`n^{th}`
bucket, :math:`\eta` would be set to :math:`n-1/2`. If positrons were
considered, they would be injected half a wake period ahead of the
location of the electrons injection position for a given period, and one
would have :math:`\eta=n-1`. The numerical cost :math:`R_t` scales as
the ratio of the total time to the shortest timescale of interest, which
is the inverse of the laser frequency, and is thus given by

.. math:: R_t=\frac{T c}{\lambda}=\frac{\left(L+\eta \lambda_p\right)}{\left(\beta_w-\beta_p\right) \lambda}

 In the laboratory, :math:`v_p=0` and the expression simplifies to

.. math:: R_{lab}=\frac{T c}{\lambda}=\frac{\left(L+\eta \lambda_p\right)}{\beta_w \lambda}

 In a frame moving at :math:`\beta c`, the quantities become

.. math::

   \begin{aligned}
   \lambda_p^*&=&\lambda_p/\left[\gamma \left(1-\beta_w \beta\right)\right] \\
   L^*&=&L/\gamma \\
   \lambda^*&=& \gamma\left(1+\beta\right) \lambda\\
   \beta_w^*&=&\left(\beta_w-\beta\right)/\left(1-\beta_w\beta\right) \\
   v_p^*&=&-\beta c \\
   T^*&=&\frac{L^*+\eta \lambda_p^*}{v_w^*-v_p^*} \\
   R_t^*&=&\frac{T^* c}{\lambda^*} = \frac{\left(L^*+\eta \lambda_p^*\right)}{\left(\beta_w^*+\beta\right) \lambda^*}\end{aligned}

 where :math:`\gamma=1/\sqrt{1-\beta^2}`.

The expected speedup from performing the simulation in a boosted frame
is given by the ratio of :math:`R_{lab}` and :math:`R_t^*`

.. math::

   S=\frac{R_{lab}}{R_t^*}=\frac{\left(1+\beta\right)\left(L+\eta \lambda_p\right)}{\left(1-\beta\beta_w\right)L+\eta \lambda_p}
   \label{Eq_scaling1d0}

We note that assuming that :math:`\beta_w\approx1` (which is a valid
approximation for most practical cases of interest) and that
:math:`\gamma<<\gamma_w`, this expression is consistent with the
expression derived earlier (**???**) for the laser-plasma acceleration
case, which states that :math:`R_t^*=\alpha R_t/\left(1+\beta\right)`
with :math:`\alpha=\left(1-\beta+l/L\right)/\left(1+l/L\right)`, where
:math:`l` is the laser length which is generally proportional to
:math:`\eta \lambda_p`, and :math:`S=R_t/R_T^*`. However, higher values
of :math:`\gamma` are of interest for maximum speedup, as shown below.

For intense lasers (:math:`a\sim 1`) typically used for acceleration,
the energy gain is limited by dephasing (**???**), which occurs over a
scale length :math:`L_d \sim \lambda_p^3/2\lambda^2`. Acceleration is
compromised beyond :math:`L_d` and in practice, the plasma length is
proportional to the dephasing length, i.e. :math:`L= \xi L_d`. In most
cases, :math:`\gamma_w^2>>1`, which allows the approximations
:math:`\beta_w\approx1-\lambda^2/2\lambda_p^2`, and
:math:`L=\xi \lambda_p^3/2\lambda^2\approx \xi \gamma_w^2 \lambda_p/2>>\eta \lambda_p`,
so that Eq.([Eq\_scaling1d0]) becomes

.. math::

   S=\left(1+\beta\right)^2\gamma^2\frac{\xi\gamma_w^2}{\xi\gamma_w^2+\left(1+\beta\right)\gamma^2\left(\xi\beta/2+2\eta\right)}
   \label{Eq_scaling1d}

 For low values of :math:`\gamma`, i.e. when :math:`\gamma<<\gamma_w`,
Eq.([Eq\_scaling1d]) reduces to

.. math::

   S_{\gamma<<\gamma_w}=\left(1+\beta\right)^2\gamma^2
   \label{Eq_scaling1d_simpl2}

 Conversely, if :math:`\gamma\rightarrow\infty`, Eq.([Eq\_scaling1d])
becomes

.. math::

   S_{\gamma\rightarrow\infty}=\frac{4}{1+4\eta/\xi}\gamma_w^2
   \label{Eq_scaling_gamma_inf}

 Finally, in the frame of the wake, i.e. when :math:`\gamma=\gamma_w`,
assuming that :math:`\beta_w\approx1`, Eq.([Eq\_scaling1d]) gives

.. math::

   S_{\gamma=\gamma_w}\approx\frac{2}{1+2\eta/\xi}\gamma_w^2
   \label{Eq_scaling_gamma_wake}

 Since :math:`\eta` and :math:`\xi` are of order unity, and the
practical regimes of most interest satisfy :math:`\gamma_w^2>>1`, the
speedup that is obtained by using the frame of the wake will be near the
maximum obtainable value given by Eq.([Eq\_scaling\_gamma\_inf]).

Note that without the use of a moving window, the relativistic effects
that are at play in the time domain would also be at play in the spatial
domain (**???**), and the :math:`\gamma^2` scaling would transform to
:math:`\gamma^4`. Hence, it is important to use a moving window even in
simulations in a Lorentz boosted frame. For very high values of the
boosted frame, the optimal velocity of the moving window may vanish
(i.e. no moving window) or even reverse.

Numerical Stability and alternate formulation in a Galilean frame
-----------------------------------------------------------------

The numerical Cherenkov instability (NCI) (**???**) is the most serious
numerical instability affecting multidimensional PIC simulations of
relativistic particle beams and streaming plasmas (**???**; Vay et al.
2010; **???**; Sironi and Spitkovsky 2011; Godfrey and Vay 2013; Xu et
al. 2013). It arises from coupling between possibly numerically
distorted electromagnetic modes and spurious beam modes, the latter due
to the mismatch between the Lagrangian treatment of particles and the
Eulerian treatment of fields Godfrey (1975).

In recent papers the electromagnetic dispersion relations for the
numerical Cherenkov instability were derived and solved for both FDTD
Godfrey and Vay (2013; Godfrey and Vay 2014) and PSATD Godfrey, Vay, and
Haber (2014b; Godfrey, Vay, and Haber 2014c) algorithms.

Several solutions have been proposed to mitigate the NCI Godfrey, Vay,
and Haber (2014a; Godfrey, Vay, and Haber 2014c; Godfrey, Vay, and Haber
2014b; Godfrey and Vay 2015; Yu, Xu, Decyk, et al. 2015; Yu, Xu,
Tableman, et al. 2015). Although these solutions efficiently reduce the
numerical instability, they typically introduce either strong smoothing
of the currents and fields, or arbitrary numerical corrections, which
are tuned specifically against the NCI and go beyond the natural
discretization of the underlying physical equation. Therefore, it is
sometimes unclear to what extent these added corrections could impact
the physics at stake for a given resolution.

For instance, NCI-specific corrections include periodically smoothing
the electromagnetic field components (**???**), using a special time
step Vay et al. (2010; **???**) or applying a wide-band smoothing of the
current components Vay et al. (2010; **???**; Vay et al. 2011). Another
set of mitigation methods involve scaling the deposited currents by a
carefully-designed wavenumber-dependent factor Godfrey and Vay (2014;
Godfrey, Vay, and Haber 2014c) or slightly modifying the ratio of
electric and magnetic fields (:math:`E/B`) before gathering their value
onto the macroparticles Godfrey, Vay, and Haber (2014b; Godfrey and Vay
2015). Yet another set of NCI-specific corrections Yu, Xu, Decyk, et al.
(2015; Yu, Xu, Tableman, et al. 2015) consists in combining a small
timestep :math:`\Delta t`, a sharp low-pass spatial filter, and a
spectral or high-order scheme that is tuned so as to create a small,
artificial “bump” in the dispersion relation Yu, Xu, Decyk, et al.
(2015). While most mitigation methods have only been applied to
Cartesian geometry, this last set of methods (Yu, Xu, Decyk, et al.
(2015; Yu, Xu, Tableman, et al. 2015)) has the remarkable property that
it can be applied Yu, Xu, Tableman, et al. (2015) to both Cartesian
geometry and quasi-cylindrical geometry (i.e. cylindrical geometry with
azimuthal Fourier decomposition Lifschitz et al. (2009; Davidson et al.
2015; Lehe et al. 2016)). However, the use of a small timestep
proportionally slows down the progress of the simulation, and the
artificial “bump” is again an arbitrary correction that departs from the
underlying physics.

A new scheme was recently proposed, in (**???**; **???**), which
completely eliminates the NCI for a plasma drifting at a uniform
relativistic velocity – with no arbitrary correction – by simply
integrating the PIC equations in *Galilean coordinates* (also known as
*comoving coordinates*). More precisely, in the new method, the Maxwell
equations *in Galilean coordinates* are integrated analytically, using
only natural hypotheses, within the PSATD framework
(Pseudo-Spectral-Analytical-Time-Domain (**???**; J. L. Vay, Haber, and
Godfrey 2013)).

The idea of the proposed scheme is to perform a Galilean change of
coordinates, and to carry out the simulation in the new coordinates:

.. math::

   \label{eq:change-var}
   {\boldsymbol{x}}' = {\boldsymbol{x}} - {{\boldsymbol{v}}_{gal}}t

 where
:math:`{\boldsymbol{x}} = x\,{\boldsymbol{u}}_x + y\,{\boldsymbol{u}}_y + z\,{\boldsymbol{u}}_z`
and
:math:`{\boldsymbol{x}}' = x'\,{\boldsymbol{u}}_x + y'\,{\boldsymbol{u}}_y + z'\,{\boldsymbol{u}}_z`
are the position vectors in the standard and Galilean coordinates
respectively.

When choosing :math:`{{\boldsymbol{v}}_{gal}}= {\boldsymbol{v}}_0`,
where :math:`{\boldsymbol{v}}_0` is the speed of the bulk of the
relativistic plasma, the plasma does not move with respect to the grid
in the Galilean coordinates :math:`{\boldsymbol{x}}'` – or,
equivalently, in the standard coordinates :math:`{\boldsymbol{x}}`, the
grid moves along with the plasma. The heuristic intuition behind this
scheme is that these coordinates should prevent the discrepancy between
the Lagrangian and Eulerian point of view, which gives rise to the NCI
Godfrey (1975).

An important remark is that the Galilean change of coordinates
([eq:change-var]) is a simple translation. Thus, when used in the
context of Lorentz-boosted simulations, it does of course preserve the
relativistic dilatation of space and time which gives rise to the
characteristic computational speedup of the boosted-frame technique.

Another important remark is that the Galilean scheme is *not* equivalent
to a moving window (and in fact the Galilean scheme can be independently
*combined* with a moving window). Whereas in a moving window, gridpoints
are added and removed so as to effectively translate the boundaries, in
the Galilean scheme the gridpoints *themselves* are not only translated
but in this case, the physical equations are modified accordingly. Most
importantly, the assumed time evolution of the current
:math:`{\boldsymbol{J}}` within one timestep is different in a standard
PSATD scheme with moving window and in a Galilean PSATD scheme
(**???**).

In the Galilean coordinates :math:`{\boldsymbol{x}}'`, the equations of
particle motion and the Maxwell equations take the form

.. math::

   \begin{aligned}
   \frac{d{\boldsymbol{x}}'}{dt} &= \frac{{\boldsymbol{p}}}{\gamma m} - {{\boldsymbol{v}}_{gal}}\label{eq:motion1} \\
   \frac{d{\boldsymbol{p}}}{dt} &= q \left( {\boldsymbol{E}} +
   \frac{{\boldsymbol{p}}}{\gamma m} \times {\boldsymbol{B}} \right) \label{eq:motion2}\\
   \left( { \frac{\partial \;}{\partial t}} - {{\boldsymbol{v}}_{gal}}\cdot{{\boldsymbol{\nabla'}}}\right){\boldsymbol{B}} &= -{{\boldsymbol{\nabla'}}}\times{\boldsymbol{E}} \label{eq:maxwell1}\\
   \frac{1}{c^2}\left( { \frac{\partial \;}{\partial t}} - {{\boldsymbol{v}}_{gal}}\cdot{{\boldsymbol{\nabla'}}}\right){\boldsymbol{E}} &= {{\boldsymbol{\nabla'}}}\times{\boldsymbol{B}} - \mu_0{\boldsymbol{J}} \label{eq:maxwell2}\end{aligned}

where :math:`{{\boldsymbol{\nabla'}}}` denotes a spatial derivative with
respect to the Galilean coordinates :math:`{\boldsymbol{x}}'`.

Integrating these equations from :math:`t=n\Delta
t` to :math:`t=(n+1)\Delta t` results in the following update equations
(see (**???**) for the details of the derivation):

.. math::

   \begin{aligned}
   {\mathbf{\tilde{B}}}^{n+1} &= \theta^2 C {\mathbf{\tilde{B}}}^n
    -\frac{\theta^2 S}{ck}i{\boldsymbol{k}}\times {\mathbf{\tilde{E}}}^n \nonumber \\
   & + \;\frac{\theta \chi_1}{\epsilon_0c^2k^2}\;i{\boldsymbol{k}} \times
                        {\mathbf{\tilde{J}}}^{n+1/2} \label{eq:disc-maxwell1}\\
   {\mathbf{\tilde{E}}}^{n+1} &=  \theta^2 C  {\mathbf{\tilde{E}}}^n
    +\frac{\theta^2 S}{k} \,c i{\boldsymbol{k}}\times {\mathbf{\tilde{B}}}^n \nonumber \\
   & +\frac{i\nu \theta \chi_1 - \theta^2S}{\epsilon_0 ck} \; {\mathbf{\tilde{J}}}^{n+1/2}\nonumber \\
   & - \frac{1}{\epsilon_0k^2}\left(\; \chi_2\;{\hat{\mathcal{\rho}}}^{n+1} -
     \theta^2\chi_3\;{\hat{\mathcal{\rho}}}^{n} \;\right) i{\boldsymbol{k}} \label{eq:disc-maxwell2}

where we used the short-hand notations
:math:`{\mathbf{\tilde{E}}}^n \equiv
{\mathbf{\tilde{E}}}({\boldsymbol{k}}, n\Delta t)`,
:math:`{\mathbf{\tilde{B}}}^n \equiv
{\mathbf{\tilde{B}}}({\boldsymbol{k}}, n\Delta t)` as well as:

.. math::

   \begin{aligned}
   &C = \cos(ck\Delta t) \quad S = \sin(ck\Delta t) \quad k
   = |{\boldsymbol{k}}| \label{eq:def-C-S}\\&
   \nu = \frac{{\boldsymbol{k}}\cdot{{\boldsymbol{v}}_{gal}}}{ck} \quad \theta =
     e^{i{\boldsymbol{k}}\cdot{{\boldsymbol{v}}_{gal}}\Delta t/2} \quad \theta^* =
     e^{-i{\boldsymbol{k}}\cdot{{\boldsymbol{v}}_{gal}}\Delta t/2} \label{eq:def-nu-theta}\\&
   \chi_1 =  \frac{1}{1 -\nu^2} \left( \theta^* -  C \theta + i
     \nu \theta S \right) \label{eq:def-chi1}\\&
   \chi_2 = \frac{\chi_1 - \theta(1-C)}{\theta^*-\theta} \quad
   \chi_3 = \frac{\chi_1-\theta^*(1-C)}{\theta^*-\theta} \label{eq:def-chi23}\end{aligned}

Note that, in the limit
:math:`{{\boldsymbol{v}}_{gal}}={\boldsymbol{0}}`, ([eq:disc-maxwell1])
and ([eq:disc-maxwell2]) reduce to the standard PSATD equations
(**???**), as expected. As shown in (**???**; **???**), the elimination
of the NCI with the new Galilean integration is verified empirically via
PIC simulations of uniform drifting plasmas and laser-driven plasma
acceleration stages, and confirmed by a theoretical analysis of the
instability.

The electromagnetic quasi-static method
=======================================

The electromagnetic quasi-static method was developed earlier than the
Lorentz boosted frame method, as a way to tackle the large separation of
length and time scales between the plasma and the driver. The
quasi-static approximation (**???**) takes advantage of the facts that
(a) the laser or particle beam driver is moving close to the speed of
light, and is hence very rigid with a slow time response, and (b) the
plasma response is extremely fast, in comparison to the driver’s. The
separation of the driver and plasma time responses enables a separation
in the treatment of the two components as follows.

Assuming the driver at a given time and position, its high rigidity
enables the approximation that it is quasi-static during the time that
it takes for traversing a transverse slice of the plasma (assumed to be
unperturbed by the driver ahead of it). The response of the plasma can
thus be computed by following the evolution of the plasma slice as the
driver propagates through it (See Fig. [fig:quasistatic]). The
reconstruction of the longitudinal and transverse structure of the wake
from the succession of transverse slices gives the full electric and
magnetic field map for evolving the beam momenta and positions on a time
scale that is commensurate with its rigidity.

Most formulations use the speed-of-light frame, defined as
:math:`\zeta=z-ct`, to follow the evolution of the plasma slices.
Assuming a slice initialized ahead of the driver, the evolution of the
plasma particles inside the slice is given by:

.. math::

   \begin{aligned}
   \frac{d\mathbf{x}_p}{d\zeta} & = & \frac{d\mathbf{x}_p}{dt}\frac{dt}{d\zeta}=\frac{\mathbf{v}_p}{v_{pz}-c},\\
   \frac{d\mathbf{p}_p}{d\zeta} & = & \frac{q}{v_{pz}-c}\left(\mathbf{E}+\mathbf{v}_p\times\mathbf{B} \right).\end{aligned}

The plasma charge and current densities are computed by accumulating the
contributions of each plasma macro-particle :math:`i`, corrected by the
time taken by the particle to cross an interval of :math:`\zeta`:

.. math::

   \begin{aligned}
   \rho_p&=&\frac{1}{\delta x \delta y \delta \zeta}\sum_i \frac{q_i}{1-v_{iz}/c}, \\
   \mathbf{J}_p&=&\frac{1}{\delta x \delta y \delta \zeta}\sum_i \frac{q_i \mathbf{v_i}}{1-v_{iz}/c}.\end{aligned}

In contrast, the evolution of a charged particle driver or witness beam
(assumed to propagate near the speed of light), is given using the
standard equations of motion:

.. math::

   \begin{aligned}
   \frac{d\mathbf{x}_{d/w}}{dt} & = &\mathbf{v}_{d/w},\\
   \frac{d\mathbf{p}_{d/w}}{dt} & = & q_{d/w}\left(\mathbf{E}+\mathbf{v}_{d/w}\times\mathbf{B} \right),\end{aligned}

while their contributions to the charge and current densities are

.. math::

   \begin{aligned}
   \rho_{d/w}&=&\frac{1}{\delta x \delta y \delta z}\sum_i q_i, \\
   \mathbf{J}_{d/w}&=&\frac{1}{\delta x \delta y \delta z}\sum_i q_i \mathbf{v_i}.\end{aligned}

The electric and magnetic fields are obtained by either solving the
equations of the scalar and vector potentials in the Coulomb or Lorentz
gauge Mora and Antonsen (1997; Huang et al. 2006) or directly the
Maxwell’s equations (**???**; Lotov 2003; Mehrling et al. 2014; An et
al. 2013) which, under the quasi-static assumption

.. math:: \partial/\partial \zeta = \partial/\partial z = -\partial/\partial ct

become

.. math::

   \begin{aligned}
   &&\nabla_\bot \times \mathbf{E}_\bot = \frac{\partial B_z}{\partial\zeta}\hat{\mathbf{z}}, \\
   &&\nabla_\bot \times E_z\hat{\mathbf{z}} = \frac{\partial \left(\mathbf{B}_\bot-\hat{\mathbf{z}}\times \mathbf{E}_\bot \right)}{\partial\zeta},\\
   &&\nabla_\bot \times \mathbf{B}_\bot -J_z \hat{\mathbf{z}} = -\frac{\partial E_z}{\partial\zeta}\hat{\mathbf{z}}, \\
   &&\nabla_\bot \times B_z\hat{\mathbf{z}} -\mathbf{J}_\bot = -\frac{\partial \left(\mathbf{E}_\bot+\hat{\mathbf{z}}\times \mathbf{B}_\bot \right)}{\partial\zeta}, \\
   &&\nabla_\bot \cdot \mathbf{E}_\bot - \rho = -\frac{\partial E_z}{\partial\zeta}, \\
   &&\nabla_\bot \cdot \mathbf{B}_\bot = -\frac{\partial B_z}{\partial\zeta}. \end{aligned}

The set of equations on the potentials or the fields can then be
rearranged in a set of 2-D Poisson-like equations that are solved
iteratively with the particle motion equations. Unlike the
Particle-In-Cell method, there is no single way of marching the set of
equations together and the reader should refer to the descriptions of
implementations in the various codes for more specific details Mora and
Antonsen (1997; Huang et al. 2006; **???**; Lotov 2003; Mehrling et al.
2014; An et al. 2013).

The Ponderomotive Guiding Center approximation
==============================================

For laser pulses with envelopes that are long compared to the laser
oscillations, it is advantageous to average over the fast laser
oscillations and solve the laser evolution with an envelope equation
Mora and Antonsen (1997; Gordon, Mori, and Antonsen 2000; Huang et al.
2006; **???**; Benedetti et al. 2012). Assuming a laser pulse in the
form of an envelope modulating a plane wave traveling at the speed of
light,

.. math::

   \begin{aligned}
   \tilde{A}_\bot = \hat{A}_\bot\left(z,\mathbf{x}_\bot,t\right)\exp{ik_0\zeta}+c.c.,\end{aligned}

the average response of a plasma to the fast laser oscillations can be
described by a ponderomotive force that inserts into a modified equation
of motion:

.. math::

   \begin{aligned}
   \frac{d\mathbf{p}}{dt} & = & q\left(\mathbf{E}+\mathbf{v}\times\mathbf{B} \right)-\frac{q^2}{\gamma mc^2}\nabla |\hat{A}_\bot|^2,\end{aligned}

with

.. math::

   \begin{aligned}
   \gamma = \sqrt{1+\frac{|\mathbf{p}|^2}{m^2c^2}+\frac{2|q\hat{A}_\bot|^2}{m^2c^4}}.\end{aligned}

Most codes Mora and Antonsen (1997; Gordon, Mori, and Antonsen 2000;
Huang et al. 2006; **???**) solve the approximate envelope equation

.. math::

   \begin{aligned}
   \left[\frac{2}{c}\frac{\partial}{\partial t}\left(ik_0+\frac{\partial}{\partial \zeta}\right) + \nabla^2_\bot\right] \hat{A}_\bot \nonumber \\
   = \frac{q^2}{mc^2} \Bigg \langle \frac{n}{\gamma}\Bigg \rangle \hat{A}_\bot \end{aligned}

while the more complete envelope equation

.. math::

   \begin{aligned}
   \left[\frac{2}{c}\frac{\partial}{\partial t}\left(ik_0+\frac{\partial}{\partial \zeta}\right) + \nabla^2_\bot - \frac{\partial^2}{\partial t^2}\right] \hat{A}_\bot \nonumber \\
   = \frac{q^2}{mc^2} \Bigg \langle \frac{n}{\gamma}\Bigg \rangle \hat{A}_\bot \end{aligned}

that retains the second time derivative is solved in the code INF&RNO
Benedetti et al. (2012). The latter equation is more exact, enabling the
accurate simulation of the laser depletion into strongly depleted
stages. As noted in Benedetti et al. (2012), in order to avoid numerical
inaccuracies, or having to grid the simulation very finely in the
longitudinal direction, it is advantageous to use the polar
representation of the laser complex field, namely
:math:`\hat{A}_\bot(\zeta)=A_\bot(\zeta)\exp[i\theta(\zeta)]`, rather
than the more common Cartesian splitting between the real and imaginary
parts
:math:`\hat{A}_\bot(\zeta)=\Re[A_\bot(\zeta)]+\imath\Im[A_\bot(\zeta)]`.
As it turns out, the functions :math:`A_\bot(\zeta)` and
:math:`\theta(\zeta)` are much smoother functions with respect to
:math:`\zeta` than :math:`\Re[A_\bot(\zeta)]` and
:math:`\Im[A_\bot(\zeta)]`, which both exhibit very short wavelength
oscillations in :math:`\zeta`, leading to more accurate numerical
differentiation along :math:`\zeta` of the polar representation at a
given longitudinal resolution.

Axi-symmetry and azimuthal Fourier decomposition
================================================

Although full PIC codes are powerful tools, which capture a wide range
of physical phenomena, they also require large computational ressources.
This is partly due to the use of a 3D Cartesian grid, which leads to a
very large number of grid cells. (Typical 3D simulations of
laser-wakefield acceleration require :math:`\sim 10^6`–:math:` 10^8`
grid cells.) For this reason, these algorithms need to be highly
parallelized, and high-resolution simulations can only be run on costly
large-scale computer facilities. However, when the driver is
cylindrically-symmetric, it is possible to take advantage of the
symmetry of the problem to reduce the computational cost of the
algorithm (**???**; Lifschitz et al. 2009; Davidson et al. 2015; Lehe et
al. 2016).

Azimuthal decomposition
-----------------------

Let us consider the fields :math:`{\boldsymbol{E}}`,
:math:`{\boldsymbol{B}}`, :math:`{\boldsymbol{J}}` and :math:`\rho` in
cylindral coordinates :math:`(r,\theta,z)`, expressed as a Fourier
series in :math:`\theta`:

.. math::

   F(r,\theta,z) = \mathrm{Re}\left[ \sum_{\ell=0}^\infty
     \tilde{F}_{\ell}(r,z) e^{-i\ell\theta} \right]
   \label{eq:chap2:azimuthal}

.. math::

   \mathrm{with} \qquad \tilde{F}_{\ell} = C_\ell \int_0^{2\pi} d\theta
   \,F(r,\theta,z)e^{i\ell\theta} \qquad
   \label{eq:chap2:Fourier-coeffs}

.. math::

   \mathrm{and} \;
   \left \{ \begin{array}{l l}
   C_{0} = 1/2\pi &\\
   C_\ell = 1/\pi &\mathrm{for}\,\ell > 0
   \end{array} \right.

 where :math:`F` represents any of the quantities :math:`E_r`,
:math:`E_\theta`, :math:`E_z`, :math:`B_r`, :math:`B_\theta`,
:math:`B_z`, :math:`J_r`, :math:`J_\theta`, :math:`J_z` are
:math:`\rho`, and where the :math:`\tilde{F}_\ell` are the associated
Fourier components (:math:`\ell` is the index of the corresponding
azimuthal mode). In the general case, this azimuthal decomposition does
not simplify the problem, since an infinity of modes have to be
considered in ([eq:chap2:azimuthal]). However, in the case of a
cylindrically-symmetric laser pulse, only the very first modes have
non-zero components. For instance, the wakefield is represented
exclusively by the mode :math:`\ell = 0`. (This is because the
quantities :math:`E_r`, :math:`E_\theta`, :math:`E_z`, :math:`B_r`,
:math:`B_\theta`, :math:`B_z`, :math:`J_r`, :math:`J_\theta`,
:math:`J_z` and :math:`\rho` associated with the wakefield are
independent of :math:`\theta`.) On the other hand, the field of the
laser pulse *does* depend on :math:`\theta`, in cylindrical coordinates.
For example, for a cylindrically-symmetric pulse propagating along
:math:`z` and polarized along
:math:`{\boldsymbol{e}}_\alpha = \cos(\alpha){\boldsymbol{e}}_x + \sin(\alpha){\boldsymbol{e}}_y`:

.. math::

   \begin{aligned}
   {\boldsymbol{E}} &= E_0(r,z){\boldsymbol{e}}_\alpha \\
   & = E_0(r,z) [\; \cos(\alpha)(\cos(\theta){\boldsymbol{e}}_r - \sin(\theta){\boldsymbol{e}}_\theta) \; \nonumber \\
   & + \; \sin(\alpha)(\sin(\theta){\boldsymbol{e}}_r + \cos(\theta){\boldsymbol{e}}_\theta) \; ]\\
   & = \mathrm{Re}[ \; E_0(r,z) e^{i\alpha} e^{-i\theta} \; ]{\boldsymbol{e}}_r \; \nonumber \\
   & + \; \mathrm{Re}[ \; -i E_0(r,z) e^{i\alpha} e^{-i\theta} \; ]{\boldsymbol{e}}_\theta.\end{aligned}

 Here the amplitude :math:`E_0` does not depend on :math:`\theta`
because the pulse was assumed to be cylindrically symmetric. In this
case, the above relation shows that the fields :math:`E_r` and
:math:`E_\theta` of the laser are represented exclusively by the mode
:math:`\ell = 1`. A similar calculation shows that the same holds for
:math:`B_r` and :math:`B_\theta`. On the whole, only the modes
:math:`\ell = 0` and :math:`\ell = 1` are a priori necessary to model
laser-wakefield acceleration. Under those conditions, the infinite sum
in ([eq:chap2:azimuthal]) is truncated at a chosen :math:`\ell_{max}`.
In principle, :math:`\ell_{max} = 1` is sufficient for laser-wakefield
acceleration. However, :math:`\ell_{max}` is kept as a free parameter in
the algorithm, in order to verify that higher modes are negligible, as
well as to allow for less-symmetric configurations. Because codes based
on this algorithm are able to take into account the modes with
:math:`\ell > 0`, they are said to be “quasi-cylindrical” (or “quasi-3D”
by some authors Davidson et al. (2015)), in contrast to cylindrical
codes, which assume that all fields are independent of :math:`\theta`,
and thus only consider the mode :math:`\ell = 0`.

Discretized Maxwell equations
-----------------------------

When the Fourier expressions of the fields are injected into the Maxwell
equations (written in cylindrical coordinates), the different azimuthal
modes decouple. In this case, the Maxwell-Ampère and Maxwell-Faraday
equations – which are needed to update the fields in the PIC cycle – can
be written separately for each azimuthal mode :math:`\ell`:

.. math::

   \begin{aligned}
   \frac{\partial \tilde{B}_{r,\ell} }{\partial t} &=
   \frac{i\ell}{r}\tilde{E}_{z,\ell} + \frac{\partial
     \tilde{E}_{\theta,\ell}}{\partial z} \\[3mm]
   \frac{\partial \tilde{B}_{\theta,\ell} }{\partial t} &=
    - \frac{\partial \tilde{E}_{r,\ell}}{\partial z} + \frac{\partial
     \tilde{E}_{z,\ell}}{\partial r} \\[3mm]
   \frac{\partial \tilde{B}_{z,\ell} }{\partial t} &=
   - \frac{1}{r} \frac{\partial (r\tilde{E}_{\theta,\ell})}{\partial r} - \frac{i\ell}{r}\tilde{E}_{r,\ell} \\[3mm]
   \frac{1}{c^2} \frac{\partial \tilde{E}_{r,\ell} }{\partial t} &=
   -\frac{i\ell}{r}\tilde{B}_{z,\ell} - \frac{\partial
     \tilde{B}_{\theta,\ell}}{\partial z} - \mu_0 \tilde{J}_{r,\ell} \\[3mm]
   \frac{1}{c^2}\frac{\partial \tilde{E}_{\theta,\ell} }{\partial t} &=
    \frac{\partial \tilde{B}_{r,\ell}}{\partial z} - \frac{\partial
     \tilde{B}_{z,\ell}}{\partial r} - \mu_0 \tilde{J}_{\theta,\ell} \\[3mm]
   \frac{1}{c^2}\frac{\partial \tilde{E}_{z,\ell} }{\partial t} &=
    \frac{1}{r} \frac{\partial (r\tilde{B}_{\theta,\ell})}{\partial r} +
    \frac{i\ell}{r}\tilde{B}_{r,\ell} - \mu_0 \tilde{J}_{z,\ell}\end{aligned}

In order to discretize these equations, each azimuthal mode is
represented on a two-dimensional grid, on which the discretized
Maxwell-Ampère and Maxwell-Faraday equations are given by

.. math::

   \begin{aligned}
   D_{t}\tilde{B}_r|_{j,\ell,k+{\frac{1}{2}}}^{n} \nonumber
   =& \frac{i\,\ell}{j\Delta r}{\tilde{E_z}^{n}_{j,\ell,k+{\frac{1}{2}}}} \\
   & + D_z \tilde{E}_{\theta}|^n_{j,\ell,k+{\frac{1}{2}}} \\
   D_{t}\tilde{B}_\theta|_{j+{\frac{1}{2}},\ell,k+{\frac{1}{2}}}^{n} \nonumber
   =& -D_z \tilde{E}_r|^n_{j+{\frac{1}{2}},\ell,k+{\frac{1}{2}}} \\
   & + D_r \tilde{E}_z|^{n}_{j+{\frac{1}{2}},\ell,k+{\frac{1}{2}}} \\
   D_{t}\tilde{B}_z|_{j+{\frac{1}{2}},\ell,k}^{n} =& \nonumber
    -\frac{(j+1){\tilde{E_\theta}^{n}_{j+1,\ell,k}} }{(j+{\frac{1}{2}})\Delta r} \\ \nonumber
    & +\frac{ j{\tilde{E_\theta}^{n}_{j,\ell,k}}}{(j+{\frac{1}{2}})\Delta r} \\
    & - \frac{i\,\ell}{(j+{\frac{1}{2}}) \Delta r}{\tilde{E_r}^{n}_{j+{\frac{1}{2}},\ell,k}} \end{aligned}

for the magnetic field components, and

.. math::

   \begin{aligned}
   \frac{1}{c^2}D_{t}\tilde{E}_r|_{j+{\frac{1}{2}},\ell,k}^{n+{\frac{1}{2}}} \nonumber
   =& -\frac{i\,\ell}{(j+{\frac{1}{2}})\Delta r}{\tilde{B_z}^{n+{\frac{1}{2}}}_{j+{\frac{1}{2}},\ell,k}} \\
   & - D_z \tilde{B}_{\theta}|^{n+{\frac{1}{2}}}_{j+{\frac{1}{2}},\ell,k} \nonumber\\
   & - \mu_0{\tilde{J_r}^{n+{\frac{1}{2}}}_{j+{\frac{1}{2}},\ell,k}} \\
   \frac{1}{c^2}D_{t}\tilde{E}_\theta|_{j,\ell,k}^{n+{\frac{1}{2}}} \nonumber
   =& D_z \tilde{B}_r|^{n+{\frac{1}{2}}}_{j,\ell,k} - D_r \tilde{B}_z|^{n+{\frac{1}{2}}}_{j,\ell,k} \\
   & - \mu_0{\tilde{J_\theta}^{n+{\frac{1}{2}}}_{j,\ell,k}} \\
   \frac{1}{c^2}D_{t}\tilde{E}_z|_{j,\ell,k+{\frac{1}{2}}}^{n+{\frac{1}{2}}} \nonumber
   =&   \frac{\left(j+{\frac{1}{2}}\right){\tilde{B_\theta}^{n+{\frac{1}{2}}}_{j+{\frac{1}{2}},\ell,k+{\frac{1}{2}}}} }{j\Delta r} \\
   =&   -\frac{\left(j-{\frac{1}{2}}\right){\tilde{B_\theta}^{n+{\frac{1}{2}}}_{j-{\frac{1}{2}},\ell,k+{\frac{1}{2}}}}}{j\Delta r} \nonumber\\
   & + \frac{i\,\ell}{j\Delta r}{\tilde{B_r}^{n+{\frac{1}{2}}}_{j,\ell,k+{\frac{1}{2}}}} \nonumber\\
   & - \mu_0{\tilde{J_z}^{n+{\frac{1}{2}}}_{j,\ell,k+{\frac{1}{2}}}}\end{aligned}

for the electric field components.

The numerical operator :math:`D_r` and :math:`D_z` are defined by

.. math::

   \begin{aligned}
   (D_r F)_{j',\ell,k'} = \frac{F_{j'+{\frac{1}{2}},\ell,k'}-F_{j'-{\frac{1}{2}},\ell,k'} }{\Delta r} \\
   (D_z F)_{j',\ell,k'} = \frac{F_{j',\ell,k'+{\frac{1}{2}}}-F_{j',\ell,k'-{\frac{1}{2}}} }{\Delta z} \\ \end{aligned}

 where :math:`j'` and :math:`k'` can be integers or half-integers.
Notice that these discretized Maxwell equations are not valid on-axis
(i.e. for :math:`j=0`), due to singularities in some of the terms.
Therefore, on the axis, these equations are replaced by specific
boundary conditions, which are based on the symmetry properties of the
fields (see Lifschitz et al. (2009) for details).

Compared to a 3D Cartesian calculation with
:math:`n_x\times n_y \times n_z` grid cells, a quasi-cylindrical
calculation with two modes (:math:`l=0` and :math:`l=1`) will require
only :math:`3 \,n_r \times n_z` grid cells. Assuming
:math:`n_x=n_y=n_r=100` as a typical transverse resolution, a
quasi-cylindrical calculation is typically over an order of magnitude
less computationally demanding than its 3D Cartesian equivalent.

Filtering
=========

It is common practice to apply digital filtering to the charge or
current density in Particle-In-Cell simulations as a complement or an
alternative to using higher order splines (**???**). A commonly used
filter in PIC simulations is the three points filter
:math:`\phi_{j}^{f}=\alpha\phi_{j}+\left(1-\alpha\right)\left(\phi_{j-1}+\phi_{j+1}\right)/2`
where :math:`\phi^{f}` is the filtered quantity. This filter is called a
bilinear filter when :math:`\alpha=0.5`. Assuming :math:`\phi=e^{jkx}`
and :math:`\phi^{f}=g\left(\alpha,k\right)e^{jkx}`, the filter gain
:math:`g` is given as a function of the filtering coefficient
:math:`\alpha` and the wavenumber :math:`k` by
:math:`g\left(\alpha,k\right)=\alpha+\left(1-\alpha\right)\cos\left(k\Delta x\right)\approx1-\left(1-\alpha\right)\frac{\left(k\Delta x\right)^{2}}{2}+O\left(k^{4}\right)`.
The total attenuation :math:`G` for :math:`n` successive applications of
filters of coefficients :math:`\alpha_{1}`...\ :math:`\alpha_{n}` is
given by
:math:`G=\prod_{i=1}^{n}g\left(\alpha_{i},k\right)\approx1-\left(n-\sum_{i=1}^{n}\alpha_{i}\right)\frac{\left(k\Delta x\right)^{2}}{2}+O\left(k^{4}\right)`.
A sharper cutoff in :math:`k` space is provided by using
:math:`\alpha_{n}=n-\sum_{i=1}^{n-1}\alpha_{i}`, so that
:math:`G\approx1+O\left(k^{4}\right)`. Such step is called a
“compensation” step (**???**). For the bilinear filter
(:math:`\alpha=1/2`), the compensation factor is
:math:`\alpha_{c}=2-1/2=3/2`. For a succession of :math:`n` applications
of the bilinear factor, it is :math:`\alpha_{c}=n/2+1`.

It is sometimes necessary to filter on a relatively wide band of
wavelength, necessitating the application of a large number of passes of
the bilinear filter or on the use of filters acting on many points. The
former can become very intensive computationally while the latter is
problematic for parallel computations using domain decomposition, as the
footprint of the filter may eventually surpass the size of subdomains. A
workaround is to use a combination of filters of limited footprint. A
solution based on the combination of three point filters with various
strides was proposed in (**???**) and operates as follows.

The bilinear filter provides complete suppression of the signal at the
grid Nyquist wavelength (twice the grid cell size). Suppression of the
signal at integer multiples of the Nyquist wavelength can be obtained by
using a stride :math:`s` in the filter
:math:`\phi_{j}^{f}=\alpha\phi_{j}+\left(1-\alpha\right)\left(\phi_{j-s}+\phi_{j+s}\right)/2`
for which the gain is given by
:math:`g\left(\alpha,k\right)=\alpha+\left(1-\alpha\right)\cos\left(sk\Delta x\right)\approx1-\left(1-\alpha\right)\frac{\left(sk\Delta x\right)^{2}}{2}+O\left(k^{4}\right)`.
For a given stride, the gain is given by the gain of the bilinear filter
shifted in k space, with the pole :math:`g=0` shifted from the
wavelength :math:`\lambda=2/\Delta x` to :math:`\lambda=2s/\Delta x`,
with additional poles, as given by
:math:`sk\Delta x=\arccos\left(\frac{\alpha}{\alpha-1}\right)\pmod{2\pi}`.
The resulting filter is pass band between the poles, but since the poles
are spread at different integer values in k space, a wide band low pass
filter can be constructed by combining filters using different strides.
As shown in (**???**), the successive application of 4-passes +
compensation of filters with strides 1, 2 and 4 has a nearly equivalent
fall-off in gain as 80 passes + compensation of a bilinear filter. Yet,
the strided filter solution needs only 15 passes of a three-point
filter, compared to 81 passes for an equivalent n-pass bilinear filter,
yielding a gain of 5.4 in number of operations in favor of the
combination of filters with stride. The width of the filter with stride
4 extends only on 9 points, compared to 81 points for a single pass
equivalent filter, hence giving a gain of 9 in compactness for the
stride filters combination in comparison to the single-pass filter with
large stencil, resulting in more favorable scaling with the number of
computational cores for parallel calculations.

Inputs and outputs
==================

Initialization of the plasma columns and drivers (laser or particle
beam) is performed via the specification of multidimensional functions
that describe the initial state with, if needed, a time dependence, or
from reconstruction of distributions based on experimental data. Care is
needed when initializing quantities in parallel to avoid double counting
and ensure smoothness of the distributions at the interface of
computational domains. When the sum of the initial distributions of
charged particles is not charge neutral, initial fields are computed
using generally a static approximation with Poisson solves accompanied
by proper relativistic scalings (**???**; Cowan et al. 2013).

Outputs include dumps of particle and field quantities at regular
intervals, histories of particle distributions moments, spectra, etc,
and plots of the various quantities. In parallel simulations, the
diagnostic subroutines need to handle additional complexity from the
domain decomposition, as well as large amount of data that may
necessitate data reduction in some form before saving to disk.

Simulations using the quasistatic method or in a Lorentz boosted frame
require additional considerations, as described below.

Outputs in parallel quasi-static simulations
--------------------------------------------

When running in parallel with domain decomposition in the direction of
propagation, the slices of the driver and wake are known at various
positions (or “stations”) in the plasma column, complicating the output
of data. The parallelization in the transverse direction (perpendicular
to the direction of propagation) uses standard domain decomposition of
the particles and fields, and does not add complexity compared to
standard PIC simulations. The parallelization in the direction of
propagation uses pipelining Feng et al. (2009), as follows.

Assuming that the grid has :math:`n_z` cells along the longitudinal
dimension and that :math:`N` processors are used for a run,
:math:`n_z/N` consecutive slices are assigned to each computational
core, as sketched in figure [Fig\_QSpipelining]. During the first
iteration, the computational core :math:`N` evolves the slices at the
front of the grid that it contains, then passes the data marching
through :math:`\zeta` to the previous core :math:`N-1` while pushing the
driver beam particle to the next station in the plasma column. During
the second iteration, both cores :math:`N` and :math:`N-1` evolve the
slices that they contain at two subsequent stations. After :math:`N`
steps, all processors are active and the procedure is repeated until the
slices on processor 1 reach the last station at the exit of the plasma
column (See Fig. [Fig\_QSpipelining]). Because each core contains a
slice that is at a different station in the plasma column, the output
modules must write data at different steps for each computational core.

Inputs and outputs in a boosted frame simulation
------------------------------------------------

|image| |image|

The input and output data are often known from, or compared to,
experimental data. Thus, calculating in a frame other than the
laboratory entails transformations of the data between the calculation
frame and the laboratory frame. This section describes the procedures
that have been implemented in the Particle-In-Cell framework Warp Grote
et al. (2005) to handle the input and output of data between the frame
of calculation and the laboratory frame (**???**). Simultaneity of
events between two frames is valid only for a plane that is
perpendicular to the relative motion of the frame. As a result, the
input/output processes involve the input of data (particles or fields)
through a plane, as well as output through a series of planes, all of
which are perpendicular to the direction of the relative velocity
between the frame of calculation and the other frame of choice.

Input in a boosted frame simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Particles -
^^^^^^^^^^^^

Particles are launched through a plane using a technique that is generic
and applies to Lorentz boosted frame simulations in general, including
plasma acceleration, and is illustrated using the case of a positively
charged particle beam propagating through a background of cold electrons
in an assumed continuous transverse focusing system, leading to a
well-known growing transverse “electron cloud” instability (**???**). In
the laboratory frame, the electron background is initially at rest and a
moving window is used to follow the beam progression. Traditionally, the
beam macroparticles are initialized all at once in the window, while
background electron macroparticles are created continuously in front of
the beam on a plane that is perpendicular to the beam velocity. In a
frame moving at some fraction of the beam velocity in the laboratory
frame, the beam initial conditions at a given time in the calculation
frame are generally unknown and one must initialize the beam
differently. However, it can be taken advantage of the fact that the
beam initial conditions are often known for a given plane in the
laboratory, either directly, or via simple calculation or projection
from the conditions at a given time in the labortory frame. Given the
position and velocity :math:`\{x,y,z,v_x,v_y,v_z\}` for each beam
macroparticle at time :math:`t=0` for a beam moving at the average
velocity :math:`v_b=\beta_b c` (where :math:`c` is the speed of light)
in the laboratory, and using the standard synchronization
(:math:`z=z'=0` at :math:`t=t'=0`) between the laboratory and the
calculation frames, the procedure for transforming the beam quantities
for injection in a boosted frame moving at velocity :math:`\beta c` in
the laboratory is as follows (the superscript :math:`'` relates to
quantities known in the boosted frame while the superscript :math:`^*`
relates to quantities that are know at a given longitudinal position
:math:`z^*` but different times of arrival):

#. project positions at :math:`z^*=0` assuming ballistic propagation

   .. math::

      \begin{aligned}
          t^* &=& \left(z-\bar{z}\right)/v_z \label{Eq:t*}\\
          x^* &=& x-v_x t^* \label{Eq:x*}\\
          y^* &=& y-v_y t^* \label{Eq:y*}\\
          z^* &=& 0 \label{Eq:z*}\end{aligned}

    the velocity components being left unchanged,

#. apply Lorentz transformation from laboratory frame to boosted frame

   .. math::

      \begin{aligned}
          t'^* &=& -\gamma t^* \label{Eq:tp*}\\
          x'^* &=& x^* \label{Eq:xp*}\\
          y'^* &=& y^* \label{Eq:yp*}\\
          z'^* &=& \gamma\beta c t^* \label{Eq:zp*}\\
          v'^*_x&=&\frac{v_x^*}{\gamma\left(1-\beta \beta_b\right)} \label{Eq:vxp*}\\
          v'^*_y&=&\frac{v_y^*}{\gamma\left(1-\beta \beta_b\right)} \label{Eq:vyp*}\\
          v'^*_z&=&\frac{v_z^*-\beta c}{1-\beta \beta_b} \label{Eq:vzp*}\end{aligned}

    where :math:`\gamma=1/\sqrt{1-\beta^2}`. With the knowledge of the
   time at which each beam macroparticle crosses the plane into
   consideration, one can inject each beam macroparticle in the
   simulation at the appropriate location and time.

#. synchronize macroparticles in boosted frame, obtaining their
   positions at a fixed :math:`t'=0` (before any particle is injected)

   .. math::

      \begin{aligned}
          z' &=& z'^*-\bar{v}'^*_z t'^* \label{Eq:zp}\end{aligned}

    This additional step is needed for setting the electrostatic or
   electromagnetic fields at the plane of injection. In a
   Particle-In-Cell code, the three-dimensional fields are calculated by
   solving the Maxwell equations (or static approximation like Poisson,
   Darwin or other (**???**)) on a grid on which the source term is
   obtained from the macroparticles distribution. This requires
   generation of a three-dimensional representation of the beam
   distribution of macroparticles at a given time before they cross the
   injection plane at :math:`z'^*`. This is accomplished by expanding
   the beam distribution longitudinally such that all macroparticles (so
   far known at different times of arrival at the injection plane) are
   synchronized to the same time in the boosted frame. To keep the beam
   shape constant, the particles are “frozen” until they cross that
   plane: the three velocity components and the two position components
   perpendicular to the boosted frame velocity are kept constant, while
   the remaining position component is advanced at the average beam
   velocity. As particles cross the plane of injection, they become
   regular “active” particles with full 6-D dynamics.

Figure [Fig\_inputoutput] (top) shows a snapshot of a beam that has
passed partly through the injection plane. As the frozen beam
macroparticles pass through the injection plane (which moves opposite to
the beam in the boosted frame), they are converted to “active"
macroparticles. The charge or current density is accumulated from the
active and the frozen particles, thus ensuring that the fields at the
plane of injection are consistent.

Laser -
^^^^^^^^

Similarly to the particle beam, the laser is injected through a plane
perpendicular to the axis of propagation of the laser (by default
:math:`z`). The electric field :math:`E_\perp` that is to be emitted is
given by the formula

.. math:: E_\perp\left(x,y,t\right)=E_0 f\left(x,y,t\right) \sin\left[\omega t+\phi\left(x,y,\omega\right)\right]

 where :math:`E_0` is the amplitude of the laser electric field,
:math:`f\left(x,y,t\right)` is the laser envelope, :math:`\omega` is the
laser frequency, :math:`\phi\left(x,y,\omega\right)` is a phase function
to account for focusing, defocusing or injection at an angle, and
:math:`t` is time. By default, the laser envelope is a three-dimensional
gaussian of the form

.. math:: f\left(x,y,t\right)=e^{-\left(x^2/2 \sigma_x^2+y^2/2 \sigma_y^2+c^2t^2/2 \sigma_z^2\right)}

 where :math:`\sigma_x`, :math:`\sigma_y` and :math:`\sigma_z` are the
dimensions of the laser pulse; or it can be defined arbitrarily by the
user at runtime. If :math:`\phi\left(x,y,\omega\right)=0`, the laser is
injected at a waist and parallel to the axis :math:`z`.

If, for convenience, the injection plane is moving at constant velocity
:math:`\beta_s c`, the formula is modified to take the Doppler effect on
frequency and amplitude into account and becomes

.. math::

   \begin{aligned}
   E_\perp\left(x,y,t\right)&=&\left(1-\beta_s\right)E_0 f\left(x,y,t\right)\nonumber \\
   &\times& \sin\left[\left(1-\beta_s\right)\omega t+\phi\left(x,y,\omega\right)\right].\end{aligned}

The injection of a laser of frequency :math:`\omega` is considered for a
simulation using a boosted frame moving at :math:`\beta c` with respect
to the laboratory. Assuming that the laser is injected at a plane that
is fixed in the laboratory, and thus moving at :math:`\beta_s=-\beta` in
the boosted frame, the injection in the boosted frame is given by

.. math::

   \begin{aligned}
   E_\perp\left(x',y',t'\right)&=&\left(1-\beta_s\right)E'_0 f\left(x',y',t'\right)\nonumber \\
   &\times&\sin\left[\left(1-\beta_s\right)\omega' t'+\phi\left(x',y',\omega'\right)\right]\\
   &=&\left(E_0/\gamma\right) f\left(x',y',t'\right) \nonumber\\
   &\times&\sin\left[\omega t'/\gamma+\phi\left(x',y',\omega'\right)\right]\end{aligned}

 since :math:`E'_0/E_0=\omega'/\omega=1/\left(1+\beta\right)\gamma`.

The electric field is then converted into currents that get injected via
a 2D array of macro-particles, with one positive and one dual negative
macro-particle for each array cell in the plane of injection, whose
weights and motion are governed by :math:`E_\perp\left(x',y',t'\right)`.
Injecting using this dual array of macroparticles offers the advantage
of automatically including the longitudinal component that arises from
emitting into a boosted frame, and to automatically verify the discrete
Gauss’ law thanks to using charge conserving (e.g. Esirkepov) current
deposition scheme Esirkepov (2001).

Output in a boosted frame simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some quantities, e.g. charge or dimensions perpendicular to the boost
velocity, are Lorentz invariant. Those quantities are thus readily
available from standard diagnostics in the boosted frame calculations.
Quantities that do not fall in this category are recorded at a number of
regularly spaced “stations", immobile in the laboratory frame, at a
succession of time intervals to record data history, or averaged over
time. A visual example is given on Fig. [Fig\_inputoutput] (bottom).
Since the space-time locations of the diagnostic grids in the laboratory
frame generally do not coincide with the space-time positions of the
macroparticles and grid nodes used for the calculation in a boosted
frame, some interpolation is performed at runtime during the data
collection process. As a complement or an alternative, selected particle
or field quantities can be dumped at regular intervals and quantities
are reconstructed in the laboratory frame during a post-processing
phase. The choice of the methods depends on the requirements of the
diagnostics and particular implementations.

Outlook
=======

The development of plasma-based accelerators depends critically on
high-performance, high-fidelity modeling to capture the full complexity
of acceleration processes that develop over a large range of space and
time scales. The field will continue to be a driver for pushing the
state-of-the-art in the detailed modeling of relativistic plasmas. The
modeling of tens of multi-GeV stages, as envisioned for plasma-based
high-energy physics colliders, will require further advances in
algorithmic, coupled to preparing the codes to take full advantage of
the upcoming generation of exascale supercomputers.

Acknowledgments
===============

This work was supported by US-DOE Contract DE-AC02-05CH11231.

This document was prepared as an account of work sponsored in part by
the United States Government. While this document is believed to contain
correct information, neither the United States Government nor any agency
thereof, nor The Regents of the University of California, nor any of
their employees, nor the authors makes any warranty, express or implied,
or assumes any legal responsibility for the accuracy, completeness, or
usefulness of any information, apparatus, product, or process disclosed,
or represents that its use would not infringe privately owned rights.
Reference herein to any specific commercial product, process, or service
by its trade name, trademark, manufacturer, or otherwise, does not
necessarily constitute or imply its endorsement, recommendation, or
favoring by the United States Government or any agency thereof, or The
Regents of the University of California. The views and opinions of
authors expressed herein do not necessarily state or reflect those of
the United States Government or any agency thereof or The Regents of the
University of California.

.. raw:: html

   <div id="refs" class="references">

.. raw:: html

   <div id="ref-Quickpic2">

An, Weiming, Viktor K. Decyk, Warren B. Mori, and Thomas M. Antonsen.
2013. “An improved iteration loop for the three dimensional quasi-static
particle-in-cell algorithm: QuickPIC.” *Journal of Computational
Physics* 250: 165–77.
doi:\ `10.1016/j.jcp.2013.05.020 <https://doi.org/10.1016/j.jcp.2013.05.020>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-AndriyashPoP2016">

Andriyash, Igor A., Remi Lehe, and Agustin Lifschitz. 2016.
“Laser-Plasma Interactions with a Fourier-Bessel Particle-in-Cell
Method.” *Physics of Plasmas* 23 (3).
doi:\ `http://dx.doi.org/10.1063/1.4943281 <https://doi.org/http://dx.doi.org/10.1063/1.4943281>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-INFERNO">

Benedetti, Carlo, Carl B. Schroeder, Eric Esarey, and Wim P. Leemans.
2012. “Efficient Modeling of Laser-plasma Accelerators Using the
Ponderomotive-based Code INF&RNO.” In *ICAP*, THAAI2.
Rostock-Warnemünde, Germany: Jacow.
http://accelconf.web.cern.ch/AccelConf/ICAP2012/papers/thaai2.pdf.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Blumenfeld2007">

Blumenfeld, Ian, Christopher E Clayton, Franz-Josef Decker, Mark J
Hogan, Chengkun Huang, Rasmus Ischebeck, Richard Iverson, et al. 2007.
“Energy doubling of 42[thinsp]GeV electrons in a metre-scale plasma
wakefield accelerator.” *Nature* 445 (7129): 741–44.
`http://dx.doi.org/10.1038/nature05538 http://www.nature.com/nature/journal/v445/n7129/suppinfo/nature05538{\\\_}S1.html <http://dx.doi.org/10.1038/nature05538 http://www.nature.com/nature/journal/v445/n7129/suppinfo/nature05538{\_}S1.html>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-BorisICNSP70">

Boris, Jp. 1970. “Relativistic Plasma Simulation-Optimization of a
Hybrid Code.” In *Proc. Fourth Conf. Num. Sim. Plasmas*, 3–67. Naval
Res. Lab., Wash., D. C.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Bruhwileraac08">

Bruhwiler, D L, J R Cary, B M Cowan, K Paul, C G R Geddes, P J
Mullowney, P Messmer, et al. 2009. “New Developments In The Simulation
Of Advanced Accelerator Concepts.” In *Aip Conference Proceedings*,
1086:29–37.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-BulanovSV2014">

Bulanov S V and Wilkens J J and Esirkepov T Zh and Korn G and Kraft G
and Kraft S D and Molls M and Khoroshkov V S. 2014. “Laser ion
acceleration for hadron therapy.” *Physics-Uspekhi* 57 (12): 1149.
http://stacks.iop.org/1063-7869/57/i=12/a=1149.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-CowanPRSTAB13">

Cowan, Benjamin M, David L Bruhwiler, John R Cary, Estelle
Cormier-Michel, and Cameron G R Geddes. 2013. “Generalized algorithm for
control of numerical dispersion in explicit time-domain electromagnetic
simulations.” *Physical Review Special Topics-Accelerators And Beams* 16
(4).
doi:\ `10.1103/PhysRevSTAB.16.041303 <https://doi.org/10.1103/PhysRevSTAB.16.041303>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-DavidsonJCP2015">

Davidson, A., A. Tableman, W. An, F.S. Tsung, W. Lu, J. Vieira, R.A.
Fonseca, L.O. Silva, and W.B. Mori. 2015. “Implementation of a hybrid
particle code with a PIC description in r–z and a gridless description
in ϕ into OSIRIS.” *Journal of Computational Physics* 281: 1063–77.
doi:\ `10.1016/j.jcp.2014.10.064 <https://doi.org/10.1016/j.jcp.2014.10.064>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-DawsonRMP83">

Dawson, J M. 1983. “Particle Simulation Of Plasmas.” *Reviews Of Modern
Physics* 55 (2): 403–47.
doi:\ `10.1103/RevModPhys.55.403 <https://doi.org/10.1103/RevModPhys.55.403>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Esirkepovcpc01">

Esirkepov, Tz. 2001. “Exact Charge Conservation Scheme For
Particle-In-Cell Simulation With An Arbitrary Form-Factor.” *Computer
Physics Communications* 135 (2): 144–53.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-QuickpicParallel">

Feng, B., C. Huang, V. Decyk, W.B. Mori, P. Muggli, and T. Katsouleas.
2009. “Enhancing parallel quasi-static particle-in-cell simulations with
a pipelining algorithm.” *Journal of Computational Physics* 228 (15):
5340–8.
doi:\ `10.1016/j.jcp.2009.04.019 <https://doi.org/10.1016/j.jcp.2009.04.019>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Godfreyjcp75">

Godfrey, Bb. 1975. “Canonical Momenta And Numerical Instabilities In
Particle Codes.” *Journal of Computational Physics* 19 (1): 58–76.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-GodfreyJCP2013">

Godfrey, Brendan B, and Jean-Luc Vay. 2013. “Numerical stability of
relativistic beam multidimensional {PIC} simulations employing the
Esirkepov algorithm.” *Journal of Computational Physics* 248 (0): 33–46.
doi:\ `http://dx.doi.org/10.1016/j.jcp.2013.04.006 <https://doi.org/http://dx.doi.org/10.1016/j.jcp.2013.04.006>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-GodfreyJCP2014">

Godfrey, Brendan B, Jean-Luc Vay, and Irving Haber. 2014a. “Numerical
stability analysis of the pseudo-spectral analytical time-domain {PIC}
algorithm.” *Journal of Computational Physics* 258 (0): 689–704.
doi:\ `http://dx.doi.org/10.1016/j.jcp.2013.10.053 <https://doi.org/http://dx.doi.org/10.1016/j.jcp.2013.10.053>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-GodfreyJCP2014_FDTD">

Godfrey, Brendan B., and Jean Luc Vay. 2014. “Suppressing the numerical
Cherenkov instability in FDTD PIC codes.” *Journal of Computational
Physics* 267: 1–6.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-GodfreyCPC2015">

———. 2015. “Improved numerical Cherenkov instability suppression in the
generalized PSTD PIC algorithm.” *Computer Physics Communications* 196.
Elsevier: 221–25.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-GodfreyJCP2014_PSATD">

Godfrey, Brendan B., Jean Luc Vay, and Irving Haber. 2014b. “Numerical
stability analysis of the pseudo-spectral analytical time-domain PIC
algorithm.” *Journal of Computational Physics* 258: 689–704.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-GodfreyIEEE2014">

———. 2014c. “Numerical stability improvements for the pseudospectral EM
PIC algorithm.” *IEEE Transactions on Plasma Science* 42 (5). Institute
of Electrical; Electronics Engineers Inc.: 1339–44.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Turbowave">

Gordon, Daniel F, W B Mori, and Thomas M Antonsen. 2000. “A
Ponderomotive Guiding Center Particle-in-Cell Code for Efficient
Modeling of Laser–Plasma Interactions.” *IEEE TRANSACTIONS ON PLASMA
SCIENCE* 28 (4).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Warp">

Grote, D P, A Friedman, J.-L. Vay, and I Haber. 2005. “The Warp Code:
Modeling High Intensity Ion Beams.” In *Aip Conference Proceedings*,
55–58. 749.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-HockneyEastwoodBook">

Hockney, R W, and J W Eastwood. 1988. *Computer simulation using
particles*. Book.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Quickpic">

Huang, C, V K Decyk, C Ren, M Zhou, W Lu, W B Mori, J H Cooley, T M
Antonsen Jr., and T Katsouleas. 2006. “Quickpic: A Highly Efficient
Particle-In-Cell Code For Modeling Wakefield Acceleration In Plasmas.”
*Journal of Computational Physics* 217 (2): 658–79.
doi:\ `10.1016/J.Jcp.2006.01.039 <https://doi.org/10.1016/J.Jcp.2006.01.039>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-LeemansPRL2014">

Leemans, W P, A J Gonsalves, H.-S. Mao, K Nakamura, C Benedetti, C B
Schroeder, Cs. Tóth, et al. 2014. “Multi-GeV Electron Beams from
Capillary-Discharge-Guided Subpetawatt Laser Pulses in the Self-Trapping
Regime.” *Phys. Rev. Lett.* 113 (24). American Physical Society: 245002.
doi:\ `10.1103/PhysRevLett.113.245002 <https://doi.org/10.1103/PhysRevLett.113.245002>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-LehePRSTAB13">

Lehe, R, A Lifschitz, C Thaury, V Malka, and X Davoine. 2013. “Numerical
growth of emittance in simulations of laser-wakefield acceleration.”
*Physical Review Special Topics-Accelerators And Beams* 16 (2).
doi:\ `10.1103/PhysRevSTAB.16.021301 <https://doi.org/10.1103/PhysRevSTAB.16.021301>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Lehe2016">

Lehe, Rémi, Manuel Kirchen, Igor A. Andriyash, Brendan B. Godfrey, and
Jean-Luc Vay. 2016. “A spectral, quasi-cylindrical and dispersion-free
Particle-In-Cell algorithm.” *Computer Physics Communications* 203:
66–82.
doi:\ `10.1016/j.cpc.2016.02.007 <https://doi.org/10.1016/j.cpc.2016.02.007>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-LewisJCP1972">

Lewis, H.Ralph. 1972. “Variational algorithms for numerical simulation
of collisionless plasma with point particles including electromagnetic
interactions.” *Journal of Computational Physics* 10 (3): 400–419.
doi:\ `http://dx.doi.org/10.1016/0021-9991(72)90044-7 <https://doi.org/http://dx.doi.org/10.1016/0021-9991(72)90044-7>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-LifschitzJCP2009">

Lifschitz, A F, X Davoine, E Lefebvre, J Faure, C Rechatin, and V Malka.
2009. “Particle-in-Cell modelling of laser{â}plasma interaction using
Fourier decomposition.” *Journal of Computational Physics* 228 (5):
1803–14.
doi:\ `http://dx.doi.org/10.1016/j.jcp.2008.11.017 <https://doi.org/http://dx.doi.org/10.1016/j.jcp.2008.11.017>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-LotovPRSTAB2003">

Lotov, K. V. 2003. “Fine wakefield structure in the blowout regime of
plasma wakefield accelerators.” *Physical Review Special Topics -
Accelerators and Beams* 6 (6). American Physical Society: 061301.
doi:\ `10.1103/PhysRevSTAB.6.061301 <https://doi.org/10.1103/PhysRevSTAB.6.061301>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Martinspac09">

Martins :math:`\backslash`\ It Et Al., S F. 2009. “Boosted Frame Pic
Simulations Of Lwfa: Towards The Energy Frontier.” In *Proc. Particle
Accelerator Conference*. Vancouver, Canada.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Hipace">

Mehrling, T, C Benedetti, C B Schroeder, and J Osterhoff. 2014. “HiPACE:
a quasi-static particle-in-cell code.” *Plasma Physics and Controlled
Fusion* 56 (8). IOP Publishing: 084012.
doi:\ `10.1088/0741-3335/56/8/084012 <https://doi.org/10.1088/0741-3335/56/8/084012>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Morapop1997">

Mora, P, and Tm Antonsen. 1997. “Kinetic Modeling Of Intense, Short
Laser Pulses Propagating In Tenuous Plasmas.” *Phys. Plasmas* 4 (1):
217–29. doi:\ `10.1063/1.872134 <https://doi.org/10.1063/1.872134>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Munzjcp2000">

Munz, Cd, P Omnes, R Schneider, E Sonnendrucker, and U Voss. 2000.
“Divergence Correction Techniques For Maxwell Solvers Based On A
Hyperbolic Model.” *Journal of Computational Physics* 161 (2): 484–511.
doi:\ `10.1006/Jcph.2000.6507 <https://doi.org/10.1006/Jcph.2000.6507>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Ohmurapiers2010">

Ohmura, Y, and Y Okamura. 2010. “Staggered Grid Pseudo-Spectral
Time-Domain Method For Light Scattering Analysis.” *Piers Online* 6 (7):
632–35.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-PukhovJPP99">

Pukhov, A. 1999. “Three-dimensional electromagnetic relativistic
particle-in-cell code VLPL (Virtual Laser Plasma Lab).” *Journal of
Plasma Physics* 61 (3): 425–33.
doi:\ `10.1017/S0022377899007515 <https://doi.org/10.1017/S0022377899007515>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Spitkovsky:Icnsp2011">

Sironi, L, and A Spitkovsky. 2011. “No Title.”

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Steinke2016">

Steinke, S, J van Tilborg, C Benedetti, C G R Geddes, C B Schroeder, J
Daniels, K K Swanson, et al. 2016. “Multistage coupling of independent
laser-plasma accelerators.” *Nature* 530 (7589). Nature Publishing
Group, a division of Macmillan Publishers Limited. All Rights Reserved.:
190–93.
`http://dx.doi.org/10.1038/nature16525 http://10.1038/nature16525 <http://dx.doi.org/10.1038/nature16525 http://10.1038/nature16525>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Vaypac09">

Vay :math:`\backslash`\ it Et Al., J.-L. 2009. “Application Of The
Reduction Of Scale Range In A Lorentz Boosted Frame To The Numerical
Simulation Of Particle Acceleration Devices.” In *Proc. Particle
Accelerator Conference*. Vancouver, Canada.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-VayAAC2010">

Vay, J -. L, C G R Geddes, C Benedetti, D L Bruhwiler, E Cormier-Michel,
B M Cowan, J R Cary, and D P Grote. 2010. “Modeling Laser Wakefield
Accelerators In A Lorentz Boosted Frame.” *Aip Conference Proceedings*
1299: 244–49.
doi:\ `10.1063/1.3520322 <https://doi.org/10.1063/1.3520322>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Vaycpc04">

Vay, J.-L., J.-C. Adam, and A Heron. 2004. “Asymmetric Pml For The
Absorption Of Waves. Application To Mesh Refinement In Electromagnetic
Particle-In-Cell Plasma Simulations.” *Computer Physics Communications*
164 (1-3): 171–77.
doi:\ `10.1016/J.Cpc.2004.06.026 <https://doi.org/10.1016/J.Cpc.2004.06.026>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Vayscidac09">

Vay, J.-L., D L Bruhwiler, C G R Geddes, W M Fawley, S F Martins, J R
Cary, E Cormier-Michel, et al. 2009. “Simulating Relativistic Beam And
Plasma Systems Using An Optimal Boosted Frame.” *Journal of Physics:
Conference Series* 180 (1): 012006 (5 Pp.).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-VayCSD12">

Vay, J.-L., D P Grote, R H Cohen, and A Friedman. 2012. “Novel methods
in the particle-in-cell accelerator code-framework warp.” Journal Paper.
*Computational Science and Discovery* 5 (1): 014019 (20 pp.).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-VayJCP2013">

Vay, Jean Luc, Irving Haber, and Brendan B. Godfrey. 2013. “A domain
decomposition method for pseudo-spectral electromagnetic simulations of
plasmas.” *Journal of Computational Physics* 243: 260–68.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-VayJCP13">

Vay, Jean-Luc, Irving Haber, and Brendan B Godfrey. 2013. “A domain
decomposition method for pseudo-spectral electromagnetic simulations of
plasmas.” *Journal of Computational Physics* 243 (June): 260–68.
doi:\ `10.1016/j.jcp.2013.03.010 <https://doi.org/10.1016/j.jcp.2013.03.010>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-VayPOPL2011">

Vay, Jl, C G R Geddes, E Cormier-Michel, and D P Grote. 2011. “Effects
Of Hyperbolic Rotation In Minkowski Space On The Modeling Of Plasma
Accelerators In A Lorentz Boosted Frame.” *Physics Of Plasmas* 18 (3):
30701. doi:\ `10.1063/1.3559483 <https://doi.org/10.1063/1.3559483>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Vincenti2016a">

Vincenti, H., and J.-L. Vay. 2016. “Detailed analysis of the effects of
stencil spatial variations with arbitrary high-order finite-difference
Maxwell solver.” *Computer Physics Communications* 200 (March). ELSEVIER
SCIENCE BV, PO BOX 211, 1000 AE AMSTERDAM, NETHERLANDS: 147–67.
doi:\ `10.1016/j.cpc.2015.11.009 <https://doi.org/10.1016/j.cpc.2015.11.009>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-XuJCP2013">

Xu, Xinlu, Peicheng Yu, Samual F Martins, Frank S Tsung, Viktor K Decyk,
Jorge Vieira, Ricardo A Fonseca, Wei Lu, Luis O Silva, and Warren B
Mori. 2013. “Numerical instability due to relativistic plasma drift in
EM-PIC simulations.” *Computer Physics Communications* 184 (11):
2503–14.
doi:\ `http://dx.doi.org/10.1016/j.cpc.2013.07.003 <https://doi.org/http://dx.doi.org/10.1016/j.cpc.2013.07.003>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Yee">

Yee, Ks. 1966. “Numerical Solution Of Initial Boundary Value Problems
Involving Maxwells Equations In Isotropic Media.” *Ieee Transactions On
Antennas And Propagation* Ap14 (3): 302–7.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Yu2016">

Yu, Peicheng, Xinlu Xu, Asher Davidson, Adam Tableman, Thamine
Dalichaouch, Fei Li, Michael D. Meyers, et al. 2016. “Enabling Lorentz
boosted frame particle-in-cell simulations of laser wakefield
acceleration in quasi-3D geometry.” *Journal of Computational Physics*.
doi:\ `10.1016/j.jcp.2016.04.014 <https://doi.org/10.1016/j.jcp.2016.04.014>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-YuCPC2015">

Yu, Peicheng, Xinlu Xu, Viktor K. Decyk, Frederico Fiuza, Jorge Vieira,
Frank S. Tsung, Ricardo A. Fonseca, Wei Lu, Luis O. Silva, and Warren B.
Mori. 2015. “Elimination of the numerical Cerenkov instability for
spectral EM-PIC codes.” *Computer Physics Communications* 192 (July).
ELSEVIER SCIENCE BV, PO BOX 211, 1000 AE AMSTERDAM, NETHERLANDS: 32–47.
doi:\ `10.1016/j.cpc.2015.02.018 <https://doi.org/10.1016/j.cpc.2015.02.018>`__.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-YuCPC2015-Circ">

Yu, Peicheng, Xinlu Xu, Adam Tableman, Viktor K. Decyk, Frank S. Tsung,
Frederico Fiuza, Asher Davidson, et al. 2015. “Mitigation of numerical
Cerenkov radiation and instability using a hybrid finite difference-FFT
Maxwell solver and a local charge conserving current deposit.” *Computer
Physics Communications* 197 (December). ELSEVIER SCIENCE BV, PO BOX 211,
1000 AE AMSTERDAM, NETHERLANDS: 144–52.
doi:\ `10.1016/j.cpc.2015.08.026 <https://doi.org/10.1016/j.cpc.2015.08.026>`__.

.. raw:: html

   </div>

.. raw:: html

   </div>

.. |image| image:: Input.pdf
   :width: 12.00000cm
.. |image| image:: Output.pdf
   :width: 12.00000cm
