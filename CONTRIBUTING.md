## Development Model

Development generally follows the following ideas:

  * New features are merged into to the `development` branch using
    Pull Requests (PRs).

    Nightly regression testing is used to ensure that no answers
    change (or if they do, that the changes were expected).

    If a change is critical, we can cherry-pick the commit from
    `development` to `master`.

  * Bug fixes, questions and contributions of new features are welcome!

       * Bugs should be reported through GitHub issues
       * We suggest asking questions through GitHub issues as well
       * *Any contributions of new features that have the potential
         to change answers should be done via pull requests.*
         A pull request should be generated from your fork of
         amrex and target the `development` branch. See below for
         details on how this process works.

         In general we squash commits upon merge to have a clean history.
         *Please ensure that your PR title and first post are descriptive,
         since these will be used for a squashed commit message.*

         Please note the following:
            If you choose to make contributions to the code 
            then you hereby grant a non-exclusive, royalty-free perpetual license 
            to install, use, modify, prepare derivative works, 
            incorporate into other computer software,
            distribute, and sublicense such enhancements or derivative works
            thereof, in binary and source code form.

  * On the first workday of each month, we perform a merge of
    `development` into `master`.  For this merge to take place, we
    need to be passing the regression tests.

    To accommodate this need, we close the merge window into
    `development` a few days before the merge day.  While the merge
    window is closed, only bug fixes should be pushed into
    `development`.  Once the merge from `development` -> `master` is
    done, the merge window reopens.

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

First, let us setup your local git repo. Make your own fork of the main
(`upstream`) repository:
on the [AMReX Github page](https://github.com/AMReX-Codes/amrex), press the
fork button. 

If you already had a fork of AMReX prior to 4/17/2020, we recommend deleting it and re-forking.
This is due to a history re-write on the main repository. Note that you will lose any branches
on your fork that haven't been merged into main development yet.

Then, you can execute:
```
# Clone your fork on your local computer. You can get this address on your fork's Github page.
git clone --branch development https://github.com/<myGithubUsername>/amrex.git

# Navigate into your repo, add a new remote for the main AMReX repo, and fetch it
cd amrex
git remote add upstream https://github.com/AMReX-Codes/amrex
git remote set-url --push upstream https://github.com/<myGithubUsername>/amrex.git
git fetch upstream

# We recommend setting your development and master branches to track the upstream ones instead of your fork:
git branch -u upstream/development
git checkout -t -b master upstream/master
```
Now you are free to play with your fork (for additional information, you can visit the
[Github fork help page](https://help.github.com/en/articles/fork-a-repo)).

> Note: you do not have to re-do the setup above every time.
> Instead, in the future, you need to update the `development` branch
> on your fork with
> ```
> git checkout development
> git pull development
> ```

Make sure you are on the `development` branch with
```
git checkout development
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
git push -u origin <branch_name>
```

If you want to synchronize your branch with the `development` branch (this is useful
when `development` is being modified while you are working on
`<branch_name>`), you can use
```
git pull upstream development
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
If you want to implement a feature and are not too sure how to split it, just open an issue about your plans and ping other AMReX developers on it to chime in.

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
git push origin --delete <branch_name>
```

## Core Developers

People who make a number of substantive contributions will be named
"core developers" of AMReX.  The criteria for becoming a core
developer are flexible, but generally involve one of the following:

  * 100 non-trivial commits to `amrex/Src/`  *and/or*

  * addition of a new algorithm / module  *and/or*

  * substantial input into the code design process or testing

If a core developer is inactive for multiple years, we may reassess their
status as a core developer.

The current list of core developers is: Ann Almgren (LBNL), Vince Beckner, John Bell (LBNL), Johannes Blaschke (LBNL), Cy Chan (LBNL), Marcus Day (LBNL), Brian Friesen (NERSC), Kevin Gott (NERSC), Daniel Graves (LBNL), Max Katz (NVIDIA), Andrew Myers (LBNL), Tan Nguyen (LBNL), Andrew Nonaka (LBNL), Michele Rosso (LBNL), Sam Williams (LBNL), Weiqun Zhang (LBNL), Michael Zingale (Stony Brook University).
