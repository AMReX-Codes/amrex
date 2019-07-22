# Contribute to WarpX

We welcome new contributors! Here is how to participate to the WarpX 
development.

## Git workflow

The WarpX project uses [git](https://git-scm.com) for version control. If you 
are new to git, you can follow one of these tutorials:
- [Learn git with bitbucket](https://www.atlassian.com/git/tutorials/learn-git-with-bitbucket-cloud)
- [git - the simple guide](http://rogerdudler.github.io/git-guide/)

### Make your own fork and create a branch on it

The basic WarpX workflow is:
- fork the main repo (or update it if you already created it);
- Implement your changes and push them on a new branch `<branch_name>` on 
your fork;
- Create a Pull Request from branch `<branch_name>` on your fork to branch 
`dev` on the main WarpX repo.

First, let us setup your local git repo. Make your own fork of the main 
(`upstream`) WarpX repo: 
on the [WarpX Github page](https://github.com/ECP-WarpX/WarpX), press the 
fork button. Then, you can execute:
```
# These 4 first lines are the same as for a standard WarpX install
mkdir warpx_directory
cd warpx_directory
git clone --branch master https://bitbucket.org/berkeleylab/picsar.git
git clone --branch development https://github.com/AMReX-Codes/amrex.git

# Clone your fork on your local computer. You can get this address on your fork's Github page.
git clone --branch dev https://github.com/<myGithubUsername>/ECP-WarpX/WarpX.git
cd warpx
# Keep track of the main WarpX repo, to remain up-to-date.
git remote add upstream https://github.com/ECP-WarpX/WarpX.git
```
Now you are free to play with your fork (for additional information, you can visit the 
[Github fork help page](https://help.github.com/en/articles/fork-a-repo)).

> Note: you do not have to re-do the setup above every time. 
> Instead, in the future, you need to update the `dev` branch
> on your fork with
> ```
> git checkout dev
> git pull upstream dev
> ```

Make sure you are on WarpX `dev` branch with
```
git checkout dev
```
in the WarpX directory.

Create a branch `<branch_name>` (the branch name should reflect the piece 
of code you want to add, like `fix_spectral_solver`) with
```
git checkout -b <branch_name>
```
and do the coding you want. It is probably a good time to look at the 
[AMReX documentation](https://amrex-codes.github.io/amrex/docs_html/) and 
at the [AMReX Doxygen reference](https://ccse.lbl.gov/pub/AMReX_Docs/index.html). 
Add the files you work on to the git staging area with 
```
git add <file_I_created> <and_file_I_modified>
```
### Commit & push your changes

Periodically commit your changes with
```
git commit -m "This is a 50-char description to explain my work"
```

The commit message (between quotation marks) is super important in order to 
follow the developments and identify bugs.

For the moment, commits are on your local repo only. You can push them to 
your fork with
```
git push -u origin <branch_name>
```

If you want to synchronize your branch with the `dev` branch (this is useful 
when the `dev` branch is being modified while you are working on 
`<branch_name>`), you can use
```
git pull upstream dev
```
and fix any conflict that may occur.

### Check that you did not break the code

Once your new feature is ready, you can check that you did not break anything. 
WarpX has automated tests running for each Pull Request. For easier debugging, 
it can be convenient to run the tests on your local machine with
```
./run_tests.sh
```
from WarpX root folder. The tests can be influenced by environment variables:
- `export WARPX_TEST_DIM=3` or `export WARPX_TEST_DIM=2` in order to select 
only the tests that correspond to this dimension
- `export WARPX_TEST_ARCH=CPU` or `export WARPX_TEST_ARCH=GPU` in order to 
run the tests on CPU or GPU respectively.
- `export WARPX_TEST_COMMIT=...` in order to test a specific commit.

### Submit a Pull Request

A Pull Request (PR) is the way to efficiently visualize the changes you made 
and to propose your new feature/improvement/fix to the WarpX project. 
Right after you push changes, a banner should appear on the Github page of 
your fork, with your `<branch_name>`. 
- Click on the `compare & pull request` button to prepare your PR. 
- Change the PR destination from `master` to `dev` (make sure that the PR is 
from `<yourFork>/<branch_name>` to `ECP-WarpX/WarpX/dev`). 
- It is time to communicate your changes: write a title and a description for 
your PR. People who review your PR are happy to know
  * what feature/fix you propose, and why
  * how you made it (created a new class than inherits from...)
  * and anything relevant to your PR (performance tests, images, *etc.*)
- Press `Create pull request`. Now you can navigate through your PR, which 
highlights the changes you made.

Pull Requests DO NOT have to be large: it is much easier to review small 
targeted PRs than a huge chunk of code, so feel free to split your work 
in small pieces.

Even before your work is ready to merge, it can be convenient to create a PR 
(so you can use Github tools to visualize your changes). In this case, please 
put the `[WIP]` tag (for Work In Progress) at the beginning of the PR title.

#### Include a test to your PR

A new feature is great, a **working** new feature is even better! Please test 
your code and add your test to the automated test suite. It's the way to 
protect your work from adventurous developers. There are three steps to follow 
to add a new automated test (illustrated here for PML boundary conditions):
- An input file for your test, in folder `Example/Tests/...`. For the PML 
test, the input file is at 
[Examples/Tests/PML/inputs2d](./Examples/Tests/PML/inputs2d). You can also 
re-use an existing input file (even better!) and pass specific parameters at 
runtime (see below).
- A Python script that reads simulation output and tests correctness versus 
theory or calibrated results. For the PML test, see
[Examples/Tests/PML/analysis_pml.py](/Examples/Tests/PML/analysis_pml.py). 
It typically ends with Python statement `assert( error<0.01 )`.
- Add an entry to [Regression/WarpX-tests.ini](./Regression/WarpX-tests.ini), 
so that a WarpX simulation runs your test in the continuous integration 
process on [Travis CI](https://docs.travis-ci.com/user/tutorial/), and the 
Python script is executed to assess the correctness. For the PML test, the 
entry is
```
[pml_x_yee]
buildDir = .
inputFile = Examples/Tests/PML/inputs2d
runtime_params = warpx.do_dynamic_scheduling=0 algo.maxwell_fdtd_solver=yee
dim = 2
addToCompileString =
restartTest = 0
useMPI = 1
numprocs = 2
useOMP = 1
numthreads = 2
compileTest = 0
doVis = 0
analysisRoutine = Examples/Tests/PML/analysis_pml_yee.py
```
If you re-use an existing input file, you can add arguments to 
`runtime_params`, like 
`runtime_params = amr.max_level=1 amr.n_cell=32 512 max_step=100 plasma_e.zmin=-200.e-6`
.

#### Include documentation to your PR

Now, let users know about your new feature by adding it to the 
[WarpX documentation](https://ecp-warpx.github.io). Our documentation uses 
[Sphinx](http://www.sphinx-doc.org/en/master/usage/quickstart.html), and it is 
located in `Docs/`. For instance, if you introduce a new runtime parameter in 
the input file, you can add it to 
[Docs/source/running_cpp/parameters.rst](Docs/source/running_cpp/parameters.rst).
If Sphinx is installed on your computer, you should be able to generate the 
html documentation with
```
make html
```
in `Docs/`. Then open `html/index.html` with your favorite web browser and look 
for your changes.

Once your code is ready with documentation and automated test, 
congratulations! you can create the PR (or remove the [WIP] tag if you already 
created it). Reviewers will interact with you if they have comments/questions.

## Style and conventions
- For indentation, WarpX uses four spaces (no tabs)
- The number of characters per line should be <80
- To define a function , for e.g., myfunction() use a space between the name of the function and the paranthesis - myfunction (). To call the function, the space is not required, i.e., just use myfunction(). The reason this is beneficial is that when we do a 'git grep ' to search for myfunction (), we can clearly see the locations where myfunction () is defined and where myfunction() is called. 
- Also, using 'git grep "myfunction ()"' searches for files only in the git repo, which is more efficient compared to the 'grep "myfunction ()"' command that searches through all the files in a directory, including plotfiles for example. 
- It is recommended that style changes are not included in the PR where new code is added. Some text editors may do this automatically and it is suggested that any automatic style changes in the text editor are removed. This is to avoid any errors that may be introduced in a PR just to do style change. 

