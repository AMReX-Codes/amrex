..
  COMMENT: This page uses the built-in python syntax highlighting to get things to
  highlight correctly in the final output.

.. _input_parameters:

Input Parameters
================

This page contains a list of AMReX input parameters. These values can be accessed by
either including them in an inputs file, or specifying them at the command line.
For more information see :ref:`sec:basics:parmparse`.

Plotfile Parameters
-------------------

The parameters below affect plotfile operations:

amr prefix
~~~~~~~~~~

.. py:data:: amr.plot_vars
   :type: string

   Specify list of variables to be included in plot.
   E.g. ``var1 var2 var3``. Use ``ALL`` to plot all variables, and
   ``NONE`` to clear all plot variables.

.. py:data:: amr.derive_plot_vars
   :type: string

   Specify list of derived variables to be included in plot.
   E.g. ``var1 var2 var3``. Use ``ALL`` to plot all variables, and
   ``NONE`` to clear all plot variables.

.. py:data:: amr.small_plot_vars
   :type: string

   Specify list of variables to be included in small plot.
   E.g. ``var1 var2 var3``. Use ``ALL`` to plot all variables, and
   ``NONE`` to clear all plot variables.

.. py:data:: amr.small_derive_plot_vars
   :type: string

   Specify list of derived variables to be included in small plot.
   E.g. ``var1 var2 var3``. Use ``ALL`` to plot all variables, and
   ``NONE`` to clear all plot variables.

.. py:data:: amr.plotfile_on_restart
   :type: bool

   Write plotfile when we restart ( only used if ``plot_int > 0``)

.. py:data:: amr.file_name_digits
   :type: int

   Specify how many zeros to use when padding the filename. E.g. ``amr.file_name_digits=5``
   with the 37th plot will have suffix, ``00037``.

.. py:data:: amr.plot_files_output
   :type: bool

   To enable or disable plotting.

.. py:data:: amr.plot_nfiles
   :type: int

   Number of processes to use/number of plotfiles to write simultaneously.

.. py:data:: amr.plot_file
   :type: string
   :value: "plt" (default)

   Plotfile root filename. E.g. ``amr.plot_file="plot"`` would give ``plot00000``.

.. py:data:: amr.plot_int
   :type: int

   Frequency of plotfile output.

.. py:data:: amr.plot_per
   :type: amrex::Real

   How often to plot in units of time. Set to negative value to disable.

.. py:data:: amr.plot_log_per
   :type: amrex::Real

   How often to plot in units of :math:`\log_{10}{}` time.
   Set to a negative value to disable.

.. py:data:: amr.small_plot_file
   :type: string
   :value: "smallplt" (default)

   Small plotfile root filename. E.g. ``amr.plot_file="plot"`` would give ``plot00000``.

.. py:data:: amr.small_plot_int
   :type: int

   Frequency of small plotfile output.

.. py:data:: amr.small_plot_per
   :type: amrex::Real

   How often to plot the small plotfile in units of time. Set to a negative
   value to disable.

.. py:data:: amr.small_plot_log_per
   :type: amrex::Real

   How often to plot the small plotfile in units of :math:`\log_{10}{}` time.
   Set to a negative value to disable.

.. py:data:: amr.write_plotfile_with_checkpoint
   :type: int

   Set to 1 to write a plotfile whenever a checkpoint is written.

.. py:data:: amr.stream_max_tries
   :type: int
   :value: 4 (default)

   Set maximum number of attempts to write plotfile. Does not apply to AsyncOut.

.. py:data:: amr.abort_on_stream_retry_failure
   :type: bool
   :value: False (default)

   If true, will AMReX will abort when a file cannot be written after exceeding
   the number of attempts specified.

.. py:data:: amr.precreateDirectories
   :type: bool
   :value: True (default)

   Write all directories at once.

.. py:data:: amr.prereadFAHeaders
   :type: bool
   :value: True (default)

   Preread and broadcast all FabArray headers.

.. py:data:: amr.plot_headerversion
   :type: int
   :value: 1 (default)

   Specify a version of the FabArray Header code. If the version does not match default,
   the specified version is used.

------


vismf prefix
~~~~~~~~~~~~

.. py:data:: vismf.v
   :type: int

   Set to 1 for verbose output of VisMF operations.

.. py:data:: vismf.headerversion
   :type: int

   Specify a version for the output file  header. If the version does not match default,
   the specified version is used.

.. py:data:: vismf.groupsets
   :type: bool
   :value: False (default)

   Wait to receive each group of processes, rather than each individual process.

..
  HACK: Needs verification

.. py:data:: vismf.setbuf
   :type: bool
   :value: True (default)

   Resize the butter to ``GetIOBufferSize``.

.. py:data:: vismf.usesingleread
   :type: bool
   :value: False (default)

   Read all FabArrays at one time.

.. py:data:: vismf.usesinglewrite
   :type: bool
   :value: False (default)

   Write all FabArrays at one time.

.. py:data:: vismf.checkfilepositions
   :type: bool
   :value: False (default)

   Verify that the bytes written match string stream.

.. py:data:: vismf.usepersistentifstreams
   :type: bool
   :value: False (default)

   *Description Needed.*

.. py:data:: vismf.usesyhchronousreads
   :type: bool
   :value: False (default)

   *Description Needed.*

.. py:data:: vismf.usedynamicsetselection
   :type: bool
   :value: True (default)

   Dynamically decide which proc writes to file.

.. py:data:: vismf.iobuffersize
   :type: Long
   :value: 262144*8 (default)

   Set the I/O buffer size.

.. py:data:: vismf.allowsparsewrites
   :type: Long
   :value: True (default)

   Allow writing when the MultiFab is determined to contain sparse data.

