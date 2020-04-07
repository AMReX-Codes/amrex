.. _developers-contributing:

Contribute to WarpX
===================

We welcome new contributors!
Here is how to participate to the WarpX development.

Git workflow
------------

The WarpX project uses `git <https://git-scm.com>`_ for version control.
If you are new to git, you can follow one of these tutorials:

- `Learn git with bitbucket <https://www.atlassian.com/git/tutorials/learn-git-with-bitbucket-cloud>`_
- `git - the simple guide <http://rogerdudler.github.io/git-guide/>`_

Make your own fork and create a branch on it
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The basic WarpX workflow is:

1. Fork the main repo (or update it if you already created it).
2. Implement your changes and push them on a new branch ``<branch_name>`` on your fork.
3. Create a Pull Request from branch ``<branch_name>`` on your fork to branch `master` on the main WarpX repo.

First, let us setup your local git repo. Make your own fork of the main (``upstream``) WarpX repo:
on the `WarpX Github page <https://github.com/ECP-WarpX/WarpX>`_, press the fork button.
Then, you can execute:

.. code-block:: sh

   # These 4 first lines are the same as for a standard WarpX install
   mkdir warpx_directory
   cd warpx_directory
   git clone --branch master https://bitbucket.org/berkeleylab/picsar.git
   git clone --branch development https://github.com/AMReX-Codes/amrex.git

   # Clone your fork on your local computer. You can get this address on your fork's Github page.
   git clone https://github.com/<myGithubUsername>/WarpX.git
   cd warpx
   # Keep track of the main WarpX repo, to remain up-to-date.
   git remote add upstream https://github.com/ECP-WarpX/WarpX.git

Now you are free to play with your fork (for additional information, you can visit the
`Github fork help page <https://help.github.com/en/articles/fork-a-repo>`_).

.. note::

   You do not have to re-do the setup above every time.
   Instead, in the future, all you need is to update the ``master`` branch on your fork with:

   .. code-block:: sh

      git checkout master
      git pull upstream master

Make sure you are on WarpX ``master`` branch with

.. code-block:: sh

   git checkout master

in the WarpX directory.

Create a branch ``<branch_name>`` (the branch name should reflect the piece of code you want to add, like ``fix_spectral_solver``) with

.. code-block:: sh

   git checkout -b <branch_name>

and do the coding you want.
It is probably a good time to look at the `AMReX documentation <https://amrex-codes.github.io/amrex/docs_html/>`_ and at the `AMReX Doxygen reference <https://ccse.lbl.gov/pub/AMReX_Docs/index.html>`_.
Add the files you work on to the git staging area with

.. code-block:: sh

   git add <file_I_created> <and_file_I_modified>

Commit & push your changes
^^^^^^^^^^^^^^^^^^^^^^^^^^

Periodically commit your changes with

.. code-block:: sh

   git commit -m "This is a 50-char description to explain my work"

The commit message (between quotation marks) is super important in order to follow the developments and identify bugs.

For the moment, commits are on your local repo only.
You can push them to your fork with

.. code-block:: sh

   git push -u origin <branch_name>

If you want to synchronize your branch with the ``master`` branch (this is useful when the ``master`` branch is being modified while you are working on ``<branch_name>``), you can use

.. code-block:: sh

   git pull upstream master

and fix any conflict that may occur.

Submit a Pull Request
^^^^^^^^^^^^^^^^^^^^^

A Pull Request (PR) is the way to efficiently visualize the changes you made and to propose your new feature/improvement/fix to the WarpX project.
Right after you push changes, a banner should appear on the Github page of your fork, with your ``<branch_name>``.

- Click on the ``compare & pull request`` button to prepare your PR.
- It is time to communicate your changes: write a title and a description for your PR.
  People who review your PR are happy to know

  * what feature/fix you propose, and why
  * how you made it (added new/edited files, created a new class than inherits from...)
  * how you tested it and what was the output you got
  * and anything else relevant to your PR (attach images and scripts, link papers, *etc.*)
- Press ``Create pull request``.
  Now you can navigate through your PR, which highlights the changes you made.

Please DO NOT write large Pull Requests, as they are very difficult and time-consuming to review.
As much as possible, split them into small targeted PRs.
For example, if find typos in the documentation open a pull request that only fixes typos.
If you want to fix a bug, make a small pull request that only fixes a bug.

If you want to implement a feature and are not too sure how to split it, just open an issue about your plans and ping other WarpX developers on it to chime in.
Generally, write helper functionality first, test it and then write implementation code.
Submit tests, documentation changes and implementation of a feature together for pull request review.

Even before your work is ready to merge, it can be convenient to create a PR (so you can use Github tools to visualize your changes).
In this case, please put the ``[WIP]`` tag (for Work-In-Progress) at the beginning of the PR title.
Another tag you may want to use is ``[mini]``, if your changes are very few lines and quick to review.
You can also use the GitHub project tab in your fork to organize the work into separate tasks/PRs and share it with the WarpX community to get feedback.

Include a test to your PR
"""""""""""""""""""""""""

A new feature is great, a **working** new feature is even better!
Please test your code and add your test to the automated test suite.
It's the way to protect your work from adventurous developers.
Instructions are given in the :ref:`testing section <developers-testing>` of our `developer's documentation <https://warpx.readthedocs.io/en/latest/developers/testing.html>`_.

Include documentation about your PR
"""""""""""""""""""""""""""""""""""

Now, let users know about your new feature by describing its usage in the `WarpX documentation <https://warpx.readthedocs.io>`_.
Our documentation uses `Sphinx <http://www.sphinx-doc.org/en/master/usage/quickstart.html>`_, and it is located in ``Docs/``.
For instance, if you introduce a new runtime parameter in the input file, you can add it to :ref:`Docs/source/running_cpp/parameters.rst <running-cpp-parameters>`.
If Sphinx is installed on your computer, you should be able to generate the html documentation with

.. code-block:: sh

   make html

in ``Docs/``. Then open ``Docs/build/html/index.html`` with your favorite web browser and look
for your changes.

Once your code is ready with documentation and automated test, congratulations!
You can create the PR (or remove the ``[WIP]`` tag if you already created it).
Reviewers will interact with you if they have comments/questions.

Style and conventions
---------------------

- For indentation, WarpX uses four spaces (no tabs)
  - Some text editors automatically modify the files you open.
    We recommend to turn on to remove trailing spaces and replace Tabs with 4 spaces.
- The number of characters per line should be <100

  - Exception: in documentation files (``.rst``/``.md``) use one sentence per line independent of its number of characters, which will allow easier edits.
- Space before and after assignment operator (``=``)
- To define a function , for e.g., ``myfunction()`` use a space between the name of the function and the paranthesis - ``myfunction ()``.
  To call the function, the space is not required, i.e., just use ``myfunction()``.

  - The reason this is beneficial is that when we do a ``git grep`` to search for ``myfunction ()``, we can clearly see the locations where ``myfunction ()`` is defined and where ``myfunction()`` is called.
  - Also, using ``git grep "myfunction ()"`` searches for files only in the git repo, which is more efficient compared to the ``grep "myfunction ()"`` command that searches through all the files in a directory, including plotfiles for example.
- It is recommended that style changes are not included in the PR where new code is added.
  This is to avoid any errors that may be introduced in a PR just to do style change.
- WarpX uses ``CamelCase`` convention for file names and class names, rather than ``snake_case``.
- The names of all member variables should be prefixed with ``m_``.
  This is particularly useful to avoid capturing member variables by value in a lambda function, which causes the whole object to be copied to GPU when running on a GPU-accelerated architecture.
  This convention should be used for all new piece of code, and it should be applied progressively to old code.
- ``#include`` directives in C++ have a distinct order to avoid bugs, see :ref:`the WarpX repo structure <developers-repo-structure>` for details
