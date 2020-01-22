.. _developers-workflow:

Workflow
========

Make a new Github release
-------------------------

WarpX has one release per month. To make the release, you need to:

    * Update the version number in all source files. There is a script for that, so you can do:

      .. code-block:: sh
          cd Tools/
          ./update_release.sh

    * Merge ``dev`` branch into ``master`` branch.

    * Click the `Draft a new release` button at https://github.com/ECP-WarpX/WarpX/releases and follow instructions.
