.. _Chap:InputsCheckpoint:

Checkpoint/Restart
==================

The following inputs must be preceded by "amr" and control checkpoint/restart.

+------------------+-----------------------------------------------------------------------+-------------+-----------+
|                  | Description                                                           |   Type      | Default   |
+==================+=======================================================================+=============+===========+
| restart          | If present, then the name of file to restart from                     |    String   | None      |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
| check_int        | Frequency of checkpoint output;                                       |    Int      | -1        |
|                  | if -1 then no checkpoints will be written                             |             |           |
+------------------+-----------------------------------------------------------------------+-------------+-----------+
| check_file       | Prefix to use for checkpoint output                                   |  String     | chk       |
+------------------+-----------------------------------------------------------------------+-------------+-----------+

