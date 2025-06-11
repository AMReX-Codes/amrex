## Development Model

Development generally follows the following ideas:

  * New features are merged into to the `development` branch using
    Pull Requests (PRs).

    Nightly regression testing is used to ensure that no answers
    change (or if they do, that the changes were expected).

  * Bug fixes, questions and contributions of new features are welcome!

       * Bugs should be reported through GitHub Issues.
       * We suggest asking questions through GitHub Discussions.
       * All contributions should be done via pull requests.
         A pull request should be generated from your fork of
         amrex and target the `development` branch. See below for
         details on how this process works.

         In general we squash commits upon merge to have a clean history.
         *Please ensure that your PR title and description are descriptive,
         since these will be used for a squashed commit message.*

         Please note the following:
            If you choose to make contributions to the code
            then you hereby grant a non-exclusive, royalty-free perpetual license
            to install, use, modify, prepare derivative works,
            incorporate into other computer software,
            distribute, and sublicense such enhancements or derivative works
            thereof, in binary and source code form.

  * On the first workday of each month, we make a tagged release.  The merge window into
    `development` is closed a few days before the release day.  While the merge window is closed,
    only bug fixes should be merged into `development`.  Once the release is done, the merge window
    reopens.

## Git workflow

AMReX uses [git](https://git-scm.com) for version control. If you
are new to git, you can follow one of these tutorials:
- [Learn git with bitbucket](https://www.atlassian.com/git/tutorials/learn-git-with-bitbucket-cloud)
- [git - the simple guide](http://rogerdudler.github.io/git-guide/)

### Make your own fork and create a branch on it

The basic workflow is:
- Fork the main repo (or update it if you already created it).
- Implement your changes and push them on a new branch `<branch_name>` on
your fork.
- Create a Pull Request from branch `<branch_name>` on your fork to branch
`development` on the main AMReX repository.

First, let us setup your local git repo. To make your own fork of the main
repository, press the fork button on the [AMReX Github page](https://github.com/AMReX-Codes/amrex).

> Note: If you already had a fork of AMReX prior to 4/17/2020, we recommend deleting it and re-forking.
> This is due to a history re-write on the main repository. Note that you will lose any branches
> on your fork that haven't been merged into main development yet.

Then, clone amrex on your local computer. If you plan on doing a lot of amrex development,
we recommend configuring your clone to use ssh access so you won't have to enter your Github
password every time, which you can do using these commands:

```
git clone git@github.com:AMReX-Codes/amrex.git
cd amrex

# Add your own fork.
# Here <Jane> is the name you give to your fork. It does not need to be your github name.
# <myGithubUsername> is your GitHub name.
git remote add <Jane> git@github.com:<myGithubUsername>/amrex.git
git fetch <Jane>

# Don't push to the main repo. Instead pushes go to your fork.
git remote set-url --push origin git@github.com:<myGithubUsername>/amrex.git
```

For instructions on setting up SSH access to your Github account on a new
machine, see
[here.](https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh)

If you instead prefer to use HTTPS authentication, configure your local clone as follows:

```
git clone https://github.com/AMReX-Codes/amrex.git
cd amrex

# Add your own fork.
# Here <Jane> is the name you give to your fork. It does not need to be your github name.
# <myGithubUsername> is your GitHub name.
git remote add <Jane> https://github.com/<myGithubUsername>/amrex.git
git fetch <Jane>

# Don't push to the main repo. Instead pushes go to your fork.
git remote set-url --push origin https://github.com/<myGithubUsername>/amrex.git
```

Now you are free to play with your fork (for additional information, you can visit the
[Github fork help page](https://help.github.com/en/articles/fork-a-repo)).

> Note: you do not have to re-do the setup above every time.

Make sure you are on the `development` branch with
```
git checkout development
git pull
```
in the AMReX directory.

Create a branch `<branch_name>` (the branch name should reflect the piece
of code you want to add, like `high_order_interpolation`) with
```
git checkout -b <branch_name>
```
and do the coding you want.
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
git push -u <Jane> <branch_name>
```

If you want to synchronize your branch with the `development` branch (this is useful
when `development` is being modified while you are working on
`<branch_name>`), you can use
```
# merge amrex main repo's development into current branch
git pull origin development
```
and fix any conflicts that may occur.

Do not merge your branch for PR into your local `development` branch,
because it will make your local `development` branch diverge from the
matching branch in the main repository after your PR is merged.

### Submit a Pull Request

A Pull Request is the way to efficiently visualize the changes you made
and to propose your new feature/improvement/fix to the AMReX project.
Right after you push changes, a banner should appear on the Github page of
your fork, with your `<branch_name>`.
- Click on the `compare & pull request` button to prepare your PR.
- It is time to communicate your changes: write a title and a description for
your PR. People who review your PR are happy to know
  * what feature/fix you propose, and why
  * how you made it (created a new class than inherits from...)
  * and anything relevant to your PR (performance tests, images, *etc.*)
- Press `Create pull request`. Now you can navigate through your PR, which
highlights the changes you made.

Please DO NOT write large Pull Requests, as they are very difficult and
time-consuming to review. As much as possible, split them into small,
targeted PRs.
For example, if find typos in the documentation open a pull request that only fixes typos.
If you want to fix a bug, make a small pull request that only fixes a bug.
If you want to implement a large feature, write helper functionality first, test it and submit those as a first pull request.
If you want to implement a feature and are not too sure how to split it,
just open a discussion about your plans and ping other AMReX developers on it to chime in.

Even before your work is ready to merge, it can be convenient to create a PR
(so you can use Github tools to visualize your changes). In this case, please
make a "draft" PR using the drop-down menu next to the "Create pull request" button.

Once your pull request is made, we will review and potentially merge it.
We recommend always creating a new branch for each pull request, as per the above instructions.
Once your pull request is merged, you can delete your local PR branch with
```
git branch -D <branch_name>
```

and you can delete the remote one on your fork with
```
git push <Jane> --delete <branch_name>
```

Generally speaking, you want to follow the following rules.

  * Do not merge your branch for PR into your local `development` branch that tracks AMReX
    `development` branch.  Otherwise your local `development` branch will diverge from AMReX
    `development` branch.

  * Do not commit in your `development` branch that tracks AMReX `development` branch.

  * Always create a new branch based off the latest `development` branch for
    each pull request, unless you are going to use git to fix it later.

If you have accidentally committed in `development` branch, you can fix it as follows,
```
git checkout -b new_branch # save your changes in a branch
git checkout development
git fetch origin
git reset --hard origin/development
```
After this, the local `development` should be in sync with AMReX `development` and your recent
commits have been saved in `new_branch` branch.

## AMReX Coding Style Guide

### Code Guidelines

AMReX developers should adhere to the following coding guidelines:
  * Indentations should use 4 spaces, not tabs.
  * Use curly braces for single statement blocks. For example:
```cpp
       for (int n=0; n<10; ++n) {
           Print() << "Like this!";
       }
```
  or
```cpp
       for (int n=0; n<10; ++n) { Print() << "Like this!"; }
```
  but not
```cpp

       for (int n=0; n<10; ++n) Print() << "Not like this.";
```
  or
```cpp
       for (int n=0; n<10; ++n)
          Print() << "Not like this.";
```
  * When declaring and defining a function, add a space after the function name and before the
parenthesis of the parameter list (but not when simply calling the function). For example:
```cpp
        void CorrectFunctionDec (int input)
```
  Not
```cpp
        void IncorrectFunctionDec(int input)
```
  This makes it easy to find where functions are defined with grep.
  * Member variables should be prefixed with `m_`. For example:
```cpp
       amrex::Real m_variable;
```
These guidelines should be adhered to in new contributions to AMReX, but
please refrain from making stylistic changes to unrelated sections of code in your PRs.

### API Documentation Using Doxygen

The Doxygen documentation is designed for advanced user-developers. It aims
to maximize the efficiency of a search-and-find style of locating information.
Doxygen style comment blocks should proceed the namespace, class, function, etc.
to be documented where appropriate. For example:
```cpp
   /**
    * \brief A one line description.
    *
    * \param[in] variable A short description of the variable.
    * \param[inout] data The value of data is read and changed.
    *
    * A longer description can be included here.
    */

    void MyFunction (int variable, MultiFab& data){
    ...
```
Additional information regarding Doxygen comment formatting can be found
in the [Doxygen Manual](https://www.doxygen.nl/manual/).
