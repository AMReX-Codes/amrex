import os
import shutil
import test_util

class Repo(object):
    """ a simple class to manage our git operations """
    def __init__(self, suite, directory, name,
                 branch_wanted=None, hash_wanted=None,
                 build=0, comp_string=None):

        self.suite = suite
        self.dir = directory
        self.name = name
        self.branch_wanted = branch_wanted
        self.hash_wanted = hash_wanted

        self.build = build   # does this repo contain build directories?
        self.comp_string = comp_string   # environment vars needed to build

        # for storage
        self.branch_orig = None
        self.hash_current = None

        self.update = True
        if hash_wanted:
            self.update = False

    def git_update(self):
        """ Do a git update of the repository.  If githash is not empty, then
            we will check out that version instead of git-pulling. """

        os.chdir(self.dir)

        # find out current branch so that we can go back later if we need.
        stdout0, stderr0, rc = test_util.run("git rev-parse --abbrev-ref HEAD")
        self.branch_orig = stdout0.rstrip('\n')

        if self.branch_orig != self.branch_wanted:
            self.suite.log.log("git checkout {} in {}".format(self.branch_wanted, self.dir))
            stdout, stderr, rc = test_util.run("git checkout {}".format(self.branch_wanted),
                                               stdin=True)

            if not rc == 0:
                self.suite.log.fail("ERROR: git checkout was unsuccessful")

        else:
            self.branch_wanted = self.branch_orig

        if self.hash_wanted == "" or self.hash_wanted is None:
            self.suite.log.log("'git pull' in {}".format(self.dir))

            # we need to be tricky here to make sure that the stdin is
            # presented to the user to get the password.
            stdout, stderr, rc = test_util.run("git pull", stdin=True,
                                               outfile="git.{}.out".format(self.name))

        else:
            stdout, stderr, rc = test_util.run("git checkout {}".format(self.hash_wanted),
                                               outfile="git.{}.out".format(self.name))

        if not rc == 0:
            self.suite.log.fail("ERROR: git update was unsuccessful")

        shutil.copy("git.{}.out".format(self.name), self.suite.full_web_dir)

    def save_head(self):

        os.chdir(self.dir)

        self.suite.log.log("saving git HEAD for {}/".format(self.name))

        stdout, stderr, rc = test_util.run("git rev-parse HEAD",
                                           outfile="git.{}.HEAD".format(self.name) )

        self.hash_current = stdout
        shutil.copy("git.{}.HEAD".format(self.name), self.suite.full_web_dir)

    def make_changelog(self):
        """ generate a ChangeLog git repository, and copy it to the
            web directory"""

        os.chdir(self.dir)

        self.suite.log.log("generating ChangeLog for {}/".format(self.name))

        test_util.run("git log --name-only",
                      outfile="ChangeLog.{}".format(self.name), outfile_mode="w")
        shutil.copy("ChangeLog.{}".format(self.name), self.suite.full_web_dir)

    def git_back(self):
        """ switch the repo back to its original branch """

        os.chdir(self.dir)
        self.suite.log.log("git checkout {} in {}".format(self.branch_orig, self.dir))

        stdout, stderr, rc = test_util.run("git checkout {}".format(self.branch_orig),
                                           stdin=True,
                                           outfile="git.{}.out".format(self.name))

        if not rc == 0:
            self.suite.log.fail("ERROR: git checkout was unsuccessful")
