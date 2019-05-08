Boundary conditions
===================

Open boundary condition for electromagnetic waves
-------------------------------------------------

For the TE case, the original Berenger’s Perfectly Matched Layer (PML) writes

.. math::

   \begin{aligned}
   \varepsilon _{0}\frac{\partial E_{x}}{\partial t}+\sigma _{y}E_{x} = & \frac{\partial H_{z}}{\partial y}\label{PML_def_1} \\
   \varepsilon _{0}\frac{\partial E_{y}}{\partial t}+\sigma _{x}E_{y} = & -\frac{\partial H_{z}}{\partial x}\label{PML_def_2} \\
   \mu _{0}\frac{\partial H_{zx}}{\partial t}+\sigma ^{*}_{x}H_{zx} = & -\frac{\partial E_{y}}{\partial x}\label{PML_def_3} \\
   \mu _{0}\frac{\partial H_{zy}}{\partial t}+\sigma ^{*}_{y}H_{zy} = & \frac{\partial E_{x}}{\partial y}\label{PML_def_4} \\
   H_{z}  = & H_{zx}+H_{zy}\label{PML_def_5}\end{aligned}

This can be generalized to

.. math::

   \begin{aligned}
   \varepsilon _{0}\frac{\partial E_{x}}{\partial t}+\sigma _{y}E_{x} = & \frac{c_{y}}{c}\frac{\partial H_{z}}{\partial y}+\overline{\sigma }_{y}H_{z}\label{APML_def_1} \\
   \varepsilon _{0}\frac{\partial E_{y}}{\partial t}+\sigma _{x}E_{y} = & -\frac{c_{x}}{c}\frac{\partial H_{z}}{\partial x}+\overline{\sigma }_{x}H_{z}\label{APML_def_2} \\
   \mu _{0}\frac{\partial H_{zx}}{\partial t}+\sigma ^{*}_{x}H_{zx} = & -\frac{c^{*}_{x}}{c}\frac{\partial E_{y}}{\partial x}+\overline{\sigma }_{x}^{*}E_{y}\label{APML_def_3} \\
   \mu _{0}\frac{\partial H_{zy}}{\partial t}+\sigma ^{*}_{y}H_{zy} = & \frac{c^{*}_{y}}{c}\frac{\partial E_{x}}{\partial y}+\overline{\sigma }_{y}^{*}E_{x}\label{APML_def_4} \\
   H_{z} = & H_{zx}+H_{zy}\label{APML_def_5}\end{aligned}

For :math:`c_{x}=c_{y}=c^{*}_{x}=c^{*}_{y}=c` and :math:`\overline{\sigma }_{x}=\overline{\sigma }_{y}=\overline{\sigma }_{x}^{*}=\overline{\sigma }_{y}^{*}=0`,
this system reduces to the Berenger PML medium, while adding the additional
constraint :math:`\sigma _{x}=\sigma _{y}=\sigma _{x}^{*}=\sigma _{y}^{*}=0`
leads to the system of Maxwell equations in vacuum.

[Sec:analytic theory, propa plane wave]Propagation of a Plane Wave in an APML Medium
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We consider a plane wave of magnitude (:math:`E_{0},H_{zx0},H_{zy0}`)
and pulsation :math:`\omega` propagating in the APML medium with an
angle :math:`\varphi` relative to the x axis

.. math::

   \begin{aligned}
   E_{x} = & -E_{0}\sin \varphi e^{i\omega \left( t-\alpha x-\beta y\right) }\label{Plane_wave_APML_def_1} \\
   E_{y} = & E_{0}\cos \varphi e^{i\omega \left( t-\alpha x-\beta y\right) }\label{Plane_wave_APML_def_2} \\
   H_{zx} = & H_{zx0}e^{i\omega \left( t-\alpha x-\beta y\right) }\label{Plane_wave_AMPL_def_3} \\
   H_{zy} = & H_{zy0}e^{i\omega \left( t-\alpha x-\beta y\right) }\label{Plane_wave_APML_def_4}\end{aligned}

where :math:`\alpha` and\ :math:`\beta` are two complex constants to
be determined.

Introducing (`[Plane_wave_APML_def_1] <#Plane_wave_APML_def_1>`__), (`[Plane_wave_APML_def_2] <#Plane_wave_APML_def_2>`__),
(`[Plane_wave_AMPL_def_3] <#Plane_wave_AMPL_def_3>`__) and (`[Plane_wave_APML_def_4] <#Plane_wave_APML_def_4>`__)
into (`[APML_def_1] <#APML_def_1>`__), (`[APML_def_2] <#APML_def_2>`__), (`[APML_def_3] <#APML_def_3>`__)
and (`[APML_def_4] <#APML_def_4>`__) gives

.. math::

   \begin{aligned}
   \varepsilon _{0}E_{0}\sin \varphi -i\frac{\sigma _{y}}{\omega }E_{0}\sin \varphi  = & \beta \frac{c_{y}}{c}\left( H_{zx0}+H_{zy0}\right) +i\frac{\overline{\sigma }_{y}}{\omega }\left( H_{zx0}+H_{zy0}\right) \label{Plane_wave_APML_1_1} \\
   \varepsilon _{0}E_{0}\cos \varphi -i\frac{\sigma _{x}}{\omega }E_{0}\cos \varphi  = & \alpha \frac{c_{x}}{c}\left( H_{zx0}+H_{zy0}\right) -i\frac{\overline{\sigma }_{x}}{\omega }\left( H_{zx0}+H_{zy0}\right) \label{Plane_wave_APML_1_2} \\
   \mu _{0}H_{zx0}-i\frac{\sigma ^{*}_{x}}{\omega }H_{zx0} = & \alpha \frac{c^{*}_{x}}{c}E_{0}\cos \varphi -i\frac{\overline{\sigma }^{*}_{x}}{\omega }E_{0}\cos \varphi \label{Plane_wave_APML_1_3} \\
   \mu _{0}H_{zy0}-i\frac{\sigma ^{*}_{y}}{\omega }H_{zy0} = & \beta \frac{c^{*}_{y}}{c}E_{0}\sin \varphi +i\frac{\overline{\sigma }^{*}_{y}}{\omega }E_{0}\sin \varphi \label{Plane_wave_APML_1_4}\end{aligned}

Defining :math:`Z=E_{0}/\left( H_{zx0}+H_{zy0}\right)` and using (`[Plane_wave_APML_1_1] <#Plane_wave_APML_1_1>`__)
and (`[Plane_wave_APML_1_2] <#Plane_wave_APML_1_2>`__), we get

.. math::

   \begin{aligned}
   \beta  = & \left[ Z\left( \varepsilon _{0}-i\frac{\sigma _{y}}{\omega }\right) \sin \varphi -i\frac{\overline{\sigma }_{y}}{\omega }\right] \frac{c}{c_{y}}\label{Plane_wave_APML_beta_of_g} \\
   \alpha  = & \left[ Z\left( \varepsilon _{0}-i\frac{\sigma _{x}}{\omega }\right) \cos \varphi +i\frac{\overline{\sigma }_{x}}{\omega }\right] \frac{c}{c_{x}}\label{Plane_wave_APML_alpha_of_g}\end{aligned}

Adding :math:`H_{zx0}` and :math:`H_{zy0}` from (`[Plane_wave_APML_1_3] <#Plane_wave_APML_1_3>`__)
and (`[Plane_wave_APML_1_4] <#Plane_wave_APML_1_4>`__) and substituting the expressions
for :math:`\alpha` and :math:`\beta` from (`[Plane_wave_APML_beta_of_g] <#Plane_wave_APML_beta_of_g>`__)
and (`[Plane_wave_APML_alpha_of_g] <#Plane_wave_APML_alpha_of_g>`__) yields

.. math::

   \begin{aligned}
   \frac{1}{Z} = & \frac{Z\left( \varepsilon _{0}-i\frac{\sigma _{x}}{\omega }\right) \cos \varphi \frac{c^{*}_{x}}{c_{x}}+i\frac{\overline{\sigma }_{x}}{\omega }\frac{c^{*}_{x}}{c_{x}}-i\frac{\overline{\sigma }^{*}_{x}}{\omega }}{\mu _{0}-i\frac{\sigma ^{*}_{x}}{\omega }}\cos \varphi \nonumber \\
    + & \frac{Z\left( \varepsilon _{0}-i\frac{\sigma _{y}}{\omega }\right) \sin \varphi \frac{c^{*}_{y}}{c_{y}}-i\frac{\overline{\sigma }_{y}}{\omega }\frac{c^{*}_{y}}{c_{y}}+i\frac{\overline{\sigma }^{*}_{y}}{\omega }}{\mu _{0}-i\frac{\sigma ^{*}_{y}}{\omega }}\sin \varphi\end{aligned}

If :math:`c_{x}=c^{*}_{x}`, :math:`c_{y}=c^{*}_{y}`, :math:`\overline{\sigma }_{x}=\overline{\sigma }^{*}_{x}`, :math:`\overline{\sigma }_{y}=\overline{\sigma }^{*}_{y}`, :math:`\frac{\sigma _{x}}{\varepsilon _{0}}=\frac{\sigma ^{*}_{x}}{\mu _{0}}` and :math:`\frac{\sigma _{y}}{\varepsilon _{0}}=\frac{\sigma ^{*}_{y}}{\mu _{0}}` then

.. math::

   \begin{aligned}
   Z = & \pm \sqrt{\frac{\mu _{0}}{\varepsilon _{0}}}\label{APML_impedance}\end{aligned}

which is the impedance of vacuum. Hence, like the PML, given some
restrictions on the parameters, the APML does not generate any reflection
at any angle and any frequency. As for the PML, this property is not
retained after discretization, as shown subsequently in this paper.

Calling :math:`\psi` any component of the field and :math:`\psi _{0}`
its magnitude, we get from (`[Plane_wave_APML_def_1] <#Plane_wave_APML_def_1>`__), (`[Plane_wave_APML_beta_of_g] <#Plane_wave_APML_beta_of_g>`__),
(`[Plane_wave_APML_alpha_of_g] <#Plane_wave_APML_alpha_of_g>`__) and (`[APML_impedance] <#APML_impedance>`__) that

.. math::

   \label{Plane_wave_absorption}
   \psi =\psi _{0}e^{i\omega \left( t\mp x\cos \varphi /c_{x}\mp y\sin \varphi /c_{y}\right) }e^{-\left( \pm \frac{\sigma _{x}\cos \varphi }{\varepsilon _{0}c_{x}}+\overline{\sigma }_{x}\frac{c}{c_{x}}\right) x}e^{-\left( \pm \frac{\sigma _{y}\sin \varphi }{\varepsilon _{0}c_{y}}+\overline{\sigma }_{y}\frac{c}{c_{y}}\right) y}

We assume that we have an APML layer of thickness :math:`\delta` (measured
along :math:`x`) and that :math:`\sigma _{y}=\overline{\sigma }_{y}=0`
and :math:`c_{y}=c.` Using (`[Plane_wave_absorption] <#Plane_wave_absorption>`__), we determine
that the coefficient of reflection given by this layer is

.. math::

   \begin{aligned}
   R_{APML}\left( \theta \right)  = & e^{-\left( \sigma _{x}\cos \varphi /\varepsilon _{0}c_{x}+\overline{\sigma }_{x}c/c_{x}\right) \delta }e^{-\left( \sigma _{x}\cos \varphi /\varepsilon _{0}c_{x}-\overline{\sigma }_{x}c/c_{x}\right) \delta }\nonumber \\
    = & e^{-2\left( \sigma _{x}\cos \varphi /\varepsilon _{0}c_{x}\right) \delta }\end{aligned}

which happens to be the same as the PML theoretical coefficient of
reflection if we assume :math:`c_{x}=c`. Hence, it follows that for
the purpose of wave absorption, the term :math:`\overline{\sigma }_{x}`
seems to be of no interest. However, although this conclusion is true
at the infinitesimal limit, it does not hold for the discretized counterpart.

Discretization
~~~~~~~~~~~~~~

.. math::

   \begin{aligned}
   \frac{E_x|^{n+1}_{j+1/2,k,l}-E_x|^{n}_{j+1/2,k,l}}{\Delta t} + \sigma_y \frac{E_x|^{n+1}_{j+1/2,k,l}+E_x|^{n}_{j+1/2,k,l}}{2} = & \frac{H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j+1/2,k-1/2,l}}{\Delta y} \\
   %
   \frac{E_y|^{n+1}_{j,k+1/2,l}-E_y|^{n}_{j,k+1/2,l}}{\Delta t} + \sigma_x \frac{E_y|^{n+1}_{j,k+1/2,l}+E_y|^{n}_{j,k+1/2,l}}{2} = & - \frac{H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j-1/2,k+1/2,l}}{\Delta x} \\
   %
   \frac{H_{zx}|^{n+3/2}_{j+1/2,k+1/2,l}-H_{zx}|^{n}_{j+1/2,k+1/2,l}}{\Delta t} + \sigma^*_x \frac{H_{zx}|^{n+3/2}_{j+1/2,k+1/2,l}+H_{zx}|^{n}_{j+1/2,k+1/2,l}}{2} = & - \frac{E_y|^{n+1}_{j+1,k+1/2,l}-E_y|^{n+1}_{j,k+1/2,l}}{\Delta x} \\
   %
   \frac{H_{zy}|^{n+3/2}_{j+1/2,k+1/2,l}-H_{zy}|^{n}_{j+1/2,k+1/2,l}}{\Delta t} + \sigma^*_y \frac{H_{zy}|^{n+3/2}_{j+1/2,k+1/2,l}+H_{zy}|^{n}_{j+1/2,k+1/2,l}}{2} = & \frac{E_x|^{n+1}_{j+1/2,k+1,l}-E_x|^{n+1}_{j+1/2,k,l}}{\Delta y} \\
   %
   H_z = & H_{zx}+H_{zy}\end{aligned}

.. math::

   \begin{aligned}
   E_x|^{n+1}_{j+1/2,k,l} = & \left(\frac{1-\sigma_y \Delta t/2}{1+\sigma_y \Delta t/2}\right) E_x|^{n}_{j+1/2,k,l} + \frac{\Delta t/\Delta y}{1+\sigma_y \Delta t/2} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j+1/2,k-1/2,l}\right) \\
   %
   E_y|^{n+1}_{j,k+1/2,l} = & \left(\frac{1-\sigma_x \Delta t/2}{1+\sigma_x \Delta t/2}\right) E_y|^{n}_{j,k+1/2,l} - \frac{\Delta t/\Delta x}{1+\sigma_x \Delta t/2} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j-1/2,k+1/2,l}\right) \\
   %
   H_{zx}|^{n+3/2}_{j+1/2,k+1/2,l} = & \left(\frac{1-\sigma^*_x \Delta t/2}{1+\sigma^*_x \Delta t/2}\right) H_{zx}|^{n}_{j+1/2,k+1/2,l} - \frac{\Delta t/\Delta x}{1+\sigma^*_x \Delta t/2} \left(E_y|^{n+1}_{j+1,k+1/2,l}-E_y|^{n+1}_{j,k+1/2,l}\right) \\
   %
   H_{zy}|^{n+3/2}_{j+1/2,k+1/2,l} = & \left(\frac{1-\sigma^*_y \Delta t/2}{1+\sigma^*_y \Delta t/2}\right) H_{zy}|^{n}_{j+1/2,k+1/2,l} + \frac{\Delta t/\Delta y}{1+\sigma^*_y \Delta t/2} \left(E_x|^{n+1}_{j+1/2,k+1,l}-E_x|^{n+1}_{j+1/2,k,l}\right) \\
   %
   H_z = & H_{zx}+H_{zy}\end{aligned}

.. math::

   \begin{aligned}
   E_x|^{n+1}_{j+1/2,k,l} = & e^{-\sigma_y\Delta t} E_x|^{n}_{j+1/2,k,l} + \frac{1-e^{-\sigma_y\Delta t}}{\sigma_y \Delta y} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j+1/2,k-1/2,l}\right) \\
   %
   E_y|^{n+1}_{j,k+1/2,l} = & e^{-\sigma_x\Delta t} E_y|^{n}_{j,k+1/2,l} - \frac{1-e^{-\sigma_x\Delta t}}{\sigma_x \Delta x} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j-1/2,k+1/2,l}\right) \\
   %
   H_{zx}|^{n+3/2}_{j+1/2,k+1/2,l} = & e^{-\sigma^*_x\Delta t} H_{zx}|^{n}_{j+1/2,k+1/2,l} - \frac{1-e^{-\sigma^*_x\Delta t}}{\sigma^*_x \Delta x} \left(E_y|^{n+1}_{j+1,k+1/2,l}-E_y|^{n+1}_{j,k+1/2,l}\right) \\
   %
   H_{zy}|^{n+3/2}_{j+1/2,k+1/2,l} = & e^{-\sigma^*_y\Delta t} H_{zy}|^{n}_{j+1/2,k+1/2,l} + \frac{1-e^{-\sigma^*_y\Delta t}}{\sigma^*_y \Delta y} \left(E_x|^{n+1}_{j+1/2,k+1,l}-E_x|^{n+1}_{j+1/2,k,l}\right) \\
   %
   H_z = & H_{zx}+H_{zy}\end{aligned}

.. math::

   \begin{aligned}
   E_x|^{n+1}_{j+1/2,k,l} = & e^{-\sigma_y\Delta t} E_x|^{n}_{j+1/2,k,l} + \frac{1-e^{-\sigma_y\Delta t}}{\sigma_y \Delta y}\frac{c_y}{c} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j+1/2,k-1/2,l}\right) \\
   %
   E_y|^{n+1}_{j,k+1/2,l} = & e^{-\sigma_x\Delta t} E_y|^{n}_{j,k+1/2,l} - \frac{1-e^{-\sigma_x\Delta t}}{\sigma_x \Delta x}\frac{c_x}{c} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j-1/2,k+1/2,l}\right) \\
   %
   H_{zx}|^{n+3/2}_{j+1/2,k+1/2,l} = & e^{-\sigma^*_x\Delta t} H_{zx}|^{n}_{j+1/2,k+1/2,l} - \frac{1-e^{-\sigma^*_x\Delta t}}{\sigma^*_x \Delta x}\frac{c^*_x}{c} \left(E_y|^{n+1}_{j+1,k+1/2,l}-E_y|^{n+1}_{j,k+1/2,l}\right) \\
   %
   H_{zy}|^{n+3/2}_{j+1/2,k+1/2,l} = & e^{-\sigma^*_y\Delta t} H_{zy}|^{n}_{j+1/2,k+1/2,l} + \frac{1-e^{-\sigma^*_y\Delta t}}{\sigma^*_y \Delta y}\frac{c^*_y}{c} \left(E_x|^{n+1}_{j+1/2,k+1,l}-E_x|^{n+1}_{j+1/2,k,l}\right) \\
   %
   H_z = & H_{zx}+H_{zy}\end{aligned}

.. math::

   \begin{aligned}
   c_x = & c e^{-\sigma_x\Delta t} \frac{\sigma_x \Delta x}{1-e^{-\sigma_x\Delta t}} \\
   c_y = & c e^{-\sigma_y\Delta t} \frac{\sigma_y \Delta y}{1-e^{-\sigma_y\Delta t}} \\
   c^*_x = & c e^{-\sigma^*_x\Delta t} \frac{\sigma^*_x \Delta x}{1-e^{-\sigma^*_x\Delta t}} \\
   c^*_y = & c e^{-\sigma^*_y\Delta t} \frac{\sigma^*_y \Delta y}{1-e^{-\sigma^*_y\Delta t}}\end{aligned}

.. math::

   \begin{aligned}
   E_x|^{n+1}_{j+1/2,k,l} = & e^{-\sigma_y\Delta t} \left[ E_x|^{n}_{j+1/2,k,l} + \frac{\Delta t}{\Delta y} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j+1/2,k-1/2,l}\right) \right] \\
   %
   E_y|^{n+1}_{j,k+1/2,l} = & e^{-\sigma_x\Delta t} \left[ E_y|^{n}_{j,k+1/2,l} - \frac{\Delta t}{\Delta x}  \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j-1/2,k+1/2,l}\right) \right] \\
   %
   H_{zx}|^{n+3/2}_{j+1/2,k+1/2,l} = & e^{-\sigma^*_x\Delta t} \left[ H_{zx}|^{n}_{j+1/2,k+1/2,l} - \frac{\Delta t}{\Delta x}  \left(E_y|^{n+1}_{j+1,k+1/2,l}-E_y|^{n+1}_{j,k+1/2,l}\right) \right] \\
   %
   H_{zy}|^{n+3/2}_{j+1/2,k+1/2,l} = & e^{-\sigma^*_y\Delta t} \left[ H_{zy}|^{n}_{j+1/2,k+1/2,l} + \frac{\Delta t}{\Delta y}  \left(E_x|^{n+1}_{j+1/2,k+1,l}-E_x|^{n+1}_{j+1/2,k,l}\right) \right] \\
   %
   H_z = & H_{zx}+H_{zy}\end{aligned}

.. math::

   \begin{aligned}
   E_x|^{n+1}_{j+1/2,k,l} = & E_x|^{n}_{j+1/2,k,l} + \frac{\Delta t}{\Delta y} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j+1/2,k-1/2,l}\right) \\
   %
   E_y|^{n+1}_{j,k+1/2,l} = & E_y|^{n}_{j,k+1/2,l} - \frac{\Delta t}{\Delta x} \left(H_z|^{n+1/2}_{j+1/2,k+1/2,l}-H_z|^{n+1/2}_{j-1/2,k+1/2,l}\right) \\
   %
   H_{zx}|^{n+3/2}_{j+1/2,k+1/2,l} = & H_{zx}|^{n}_{j+1/2,k+1/2,l} - \frac{\Delta t}{\Delta x} \left(E_y|^{n+1}_{j+1,k+1/2,l}-E_y|^{n+1}_{j,k+1/2,l}\right) \\
   %
   H_{zy}|^{n+3/2}_{j+1/2,k+1/2,l} = & H_{zy}|^{n}_{j+1/2,k+1/2,l} + \frac{\Delta t}{\Delta y} \left(E_x|^{n+1}_{j+1/2,k+1,l}-E_x|^{n+1}_{j+1/2,k,l}\right) \\
   %
   H_z = & H_{zx}+H_{zy}\end{aligned}
