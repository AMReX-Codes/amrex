.. _Chap:InputsComputeBackends:

Compute Backends
================

The following inputs must be preceded by ``amrex.`` and determine runtime options of CPU or GPU compute implementations.

+------------------------+-----------------------------------------------------------------------+-------------+------------+
| Parameter              | Description                                                           |   Type      | Default    |
+========================+=======================================================================+=============+============+
| ``omp_threads``        | If OpenMP is enabled, this can be used to set the default number of   |   String    | ``system`` |
|                        | threads. The special value ``nosmt`` can be used to avoid using       |   or Int    |            |
|                        | threads for virtual cores (aka Hyperthreading or SMT), as is default  |             |            |
|                        | in OpenMP, and instead only spawns threads equal to the number of     |             |            |
|                        | physical cores in the system.                                         |             |            |
|                        | For the values ``system`` and ``nosmt``, the environment variable     |             |            |
|                        | ``OMP_NUM_THREADS`` takes precedence. For Integer values,             |             |            |
|                        | ``OMP_NUM_THREADS`` is ignored.                                       |             |            |
+------------------------+-----------------------------------------------------------------------+-------------+------------+

For GPU-specific parameters, see also the :ref:`GPU chapter <sec:gpu:parameters>`.
