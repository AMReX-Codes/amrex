import io
import os
import subprocess
import sys

def run(command, stdin=False, outfile=None):
    """ run a command in the unix shell """

    sin = None
    if stdin: sin = subprocess.PIPE
    p0 = subprocess.Popen(command, stdin=sin, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, shell=True)

    stdout0 = p0.communicate()
    if stdin: p0.stdin.close()
    rc = p0.returncode
    p0.stdout.close()

    if outfile is not None:
        try: cf = io.open(outfile, "w", encoding="latin-1")
        except IOError:
            sys.exit("ERROR: unable to open file for writing: {}".format(outfile))
        else:
            for line in stdout0:
                if line is not None:
                    cf.write(line.decode('latin-1'))
            cf.close()

    return stdout0, rc


class Preprocessor(object):
    """ hold the information about preprocessing """

    def __init__(self, temp_dir=None, cpp_cmd=None,
                 defines=None, f90_preprocess=None):

        self.temp_dir = temp_dir
        self.cpp_cmd = cpp_cmd
        self.defines = defines
        self.f90_preprocess = f90_preprocess

    def preprocess(self, sf, add_name="F90PP"):
        """ preprocess the file described by a SourceFile object sf """

        # we want to do:
        # $(FORT_CPP) $(CPPFLAGS) $< | $(F90PREP) > $(f77TempDir)/$*.f90
        # we output to the temporary directory

        processed_file = "{}/{}-{}".format(self.temp_dir, add_name,
                                           os.path.basename(sf.name))

        if self.f90_preprocess != "" and self.f90_preprocess is not None:
            command = "{} {} {} | {}".format(self.cpp_cmd, self.defines,
                                             sf.name, self.f90_preprocess)
        else:
            command = "{} {} {}".format(self.cpp_cmd, self.defines,
                                        sf.name)

        stdout, rc = run(command, outfile=processed_file)

        #if rc == 0:
        sf.cpp_name = processed_file
        #else:
        #    raise ValueError("cpp process failed for {}".format(sf.name))

        return command

