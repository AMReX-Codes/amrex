# Contribute to WarpX

We welcome new contributors! Here is how to participate to the WarpX 
development.

## Git workflow

The WarpX project uses [git](https://git-scm.com) for version control. If you 
are new to git, you can follow one of these tutorials:
- [Learn git with bitbucket](https://www.atlassian.com/git/tutorials/learn-git-with-bitbucket-cloud)
- [git - the simple guide](http://rogerdudler.github.io/git-guide/)

### Create your branch from the `dev` branch

Follow instructions from the 
[WarpX documentation](https://ecp-warpx.github.io/doc_versions/dev/building/building.html) 
to install the code. Make sure you are on WarpX `dev` branch: run
```
git checkout dev
```
in the WarpX directory.

Create a branch `<branch_name>` (the branch name should reflect the piece 
of code you want to add, like `fix_spectral_solver`) with
```
git checkout -b <branch_name>
```
and do the hacks you want. Add the files you modified/added to the git staging 
area with 
```
git add file_I_created and_file_I_modified
```

### Commit & push your changes

Periodically commit your changes with
```
git commit -m "Write 50 chars description to explain your work"
```
The commit message (between quotation marks) is super important to follow the 
developments and identify bugs.

For the moment, commits are on your local repo only. You can push them to 
the WarpX Github repo with
```
git push -u origin <branch_name>
```

If you want to synchronize your branch with the `dev` branch (this is useful) 
when the `dev` branch is modified while you are working on `<branch_name>`, 
you can use
```
git pull --rebase origin dev
```

### Check that you did not break the code

Once your new feature is ready, you can check that you did not break anything. 
WarpX has automated tests running at each `git push`. For easier debugging, 
it is convenient to run the tests on your local machine. The code can be 
tested by running
```
./run_tests.sh
```
from the root folder of WarpX (after downloading the sources of `amrex` and 
`picsar`, as explained in the documentation).

The tests can be influenced by environment variables:
- `export WARPX_TEST_DIM=3` or `export WARPX_TEST_DIM=2` in order to select 
only the tests that correspond to this dimension
- `export WARPX_TEST_ARCH=CPU` or `export WARPX_TEST_ARCH=GPU` in order to 
run the tests on CPU or GPU respectively.
- `export WARPX_TEST_COMMIT=...` in order to test a specific commit.

### Submit a Pull Request

#### Include a test to your PR

#### Include documentation to your PR

## Style and conventions