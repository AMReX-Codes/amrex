#!/usr/bin/env python

"""
Search the C++ source code for Fortran calls marked with
#pragma gpu and maintain a list of these.

Then copy the C++ headers for Fortran files (typically *_F.H) into a
temp directory, modifying any subroutines marked with #pragma gpu to
have both a device and host signature.

Finally, modify the C++ source files to replace the #pragma gpu
with a CUDA kernel launch.
"""

from __future__ import print_function

import sys

if sys.version_info < (2, 7):
    sys.exit("ERROR: need python 2.7 or later for dep.py")

if sys.version[0] == "2":
    reload(sys)
    sys.setdefaultencoding('utf8')

import os
import re
import argparse
import shutil
import find_files_vpath as ffv
import preprocess
from multiprocessing import Pool


TEMPLATE = """
__global__ static void cuda_{}
{{
{}
{}
   int blo[3];
   int bhi[3];
   for (int k = lo[2] + blockIdx.z * blockDim.z + threadIdx.z; k <= hi[2]; k += blockDim.z * gridDim.z) {{
     blo[2] = k;
     bhi[2] = k;
     for (int j = lo[1] + blockIdx.y * blockDim.y + threadIdx.y; j <= hi[1]; j += blockDim.y * gridDim.y) {{
       blo[1] = j;
       bhi[1] = j;
       for (int i = lo[0] + blockIdx.x * blockDim.x + threadIdx.x; i <= hi[0]; i += blockDim.x * gridDim.x) {{
         blo[0] = i;
         bhi[0] = i;
         {};
       }}
     }}
   }}
}}
"""

# for finding just the variable definitions in the function signature (between the ())
decls_re = re.compile("(.*?)(\\()(.*)(\\))", re.IGNORECASE|re.DOTALL)

# for finding a header entry for a function binding to a fortran subroutine
fortran_binding_re = re.compile("(void)(\\s+)((?:[a-z0-9][a-z_0-9]+))",
                                re.IGNORECASE|re.DOTALL)

class HeaderFile(object):
    """ hold information about one of the headers """

    def __init__(self, filename):

        self.name = filename

        # when we preprocess, the output has a different name
        self.cpp_name = None


class CppFile(object):
    """ hold information about one of the C++ files """

    def __init__(self, filename):

        self.name = filename

        self.cpp_name = None



def find_targets_from_pragmas(outdir, cxx_files, macro_list, cpp):
    """read through the C++ files and look for the functions marked with
    #pragma gpu -- these are the routines we intend to offload (the targets),
    so make a list of them"""

    targets = dict()

    for c in cxx_files:
        cxx = "/".join([c[1], c[0]])

        # Make a temporary copy of the original C++ file,
        # and do a bit of manual preprocessing. We want to
        # preprocess out all of the -D defines, but keep
        # things like BL_TO_FORTRAN_ANYD. A trick to do this
        # is to manually go through and remove all of the #include
        # statements, so that the file doesn't know how to replace
        # those macros, and then run it through cpp to remove
        # everything else not defined with -D.

        temp_cxx = "/".join([outdir, "temp-" + c[0]])
        shutil.copyfile(cxx, temp_cxx)

        noinclude_cxx = "/".join([outdir, "temp-noinclude-" + c[0]])

        try:
            hin = open(temp_cxx, "r")
        except IOError:
            sys.exit("Cannot open temporary C++ file {}".format(temp_cxx))

        try:
            hout = open(noinclude_cxx, "w")
        except IOError:
            sys.exit("Cannot open temporary C++ file {}".format(noinclude_cxx))

        line = hin.readline()
        while line:
            if not ("#" in line and "include" in line):
                hout.write(line)
            line = hin.readline()

        hin.close()
        hout.close()
        os.remove(temp_cxx)

        # Now run the noinclude file through cpp.

        cppf = CppFile(noinclude_cxx)

        # preprocess -- this will create a new file in our temp_dir that
        # was run through cpp and has the name CPP-filename
        cpp.preprocess(cppf, add_name="CPP")
        preprocessed_cxx = "/".join([outdir, "CPP-temp-noinclude-" + c[0]])

        os.remove(noinclude_cxx)

        # Now re-open the preprocessed C++ file as input
        try:
            hin = open(preprocessed_cxx, "r")
        except IOError:
            sys.exit("Cannot open temporary C++ file {}".format(preprocessed_cxx))

        # look for the appropriate pragma, and once found, capture the
        # function call following it
        line = hin.readline()
        while line:

            # if the line starts with "#pragma gpu", then we need
            # to take action
            if line.startswith("#pragma gpu"):
                # we don't need to reproduce the pragma line in the
                # output, but we need to capture the whole function
                # call that follows
                func_call = ""
                line = hin.readline()
                while not line.strip().endswith(";"):
                    func_call += line
                    line = hin.readline()
                # last line -- remove the semi-colon
                func_call += line.rstrip()[:-1]

                # now split it into the function name and the
                # arguments
                dd = decls_re.search(func_call)
                func_name = dd.group(1).strip().replace(" ", "")

                # Now, for each argument in the function, record
                # if it has a special macro in it. Save the position
                # in the argument list corresponding to this macro.
                # We'll use this later when updating the signature
                # in the header file.

                targets[func_name] = []
                args = dd.group().strip().split('(', 1)[1].rsplit(')', 1)[0].split(',')

                # Some of the arguments may have commas inside them, e.g.
                # BL_TO_FORTRAN_N(fab, comp). We want to make sure these only
                # look like one argument. So we'll check if there are parentheses
                # in the argument, and join the arguments together if so. This
                # is only sufficient to capture one comma, it needs to be expanded
                # if the argument could have multiple commas, e.g. a macro with
                # three arguments.

                for i, arg in enumerate(args):
                    if 'BL_TO_FORTRAN_N_ANYD' in arg:
                        args[i:i+2] = [''.join(args[i:i+2])]

                for j, macro in enumerate(macro_list):
                    targets[func_name].append([])

                    # Keep a running counter of which argument index we are at.
                    # This has to account for macros which may expand into
                    # multiple arguments. At this time we only know about
                    # the specific macros in the macro_list.

                    i = 0
                    for arg in args:
                        if macro in arg:
                            targets[func_name][j].append(i)

                        if 'BL_TO_FORTRAN_ANYD' in arg:
                            i += 3
                        elif 'BL_TO_FORTRAN_N_ANYD' in arg:
                            i += 3
                        elif 'BL_TO_FORTRAN_FAB' in arg:
                            i += 4
                        elif 'BL_TO_FORTRAN_BOX' in arg:
                            i += 2
                        else:
                            i += 1

            line = hin.readline()

        os.remove(preprocessed_cxx)

    return targets


def convert_headers(inputs):
    """rewrite the C++ headers that contain the Fortran routines"""

    header_file = inputs[0]
    outdir      = inputs[1]
    targets     = inputs[2]
    macro_list  = inputs[3]
    cpp         = inputs[4]

    print('looking for targets: {}'.format(list(targets)))
    print('looking in header file: {}'.format(header_file))

    # first preprocess the header and store it in a temporary
    # location.  The preprocessed header will only be used for the
    # search for the signature, not as the basis for writing the final
    # CUDA header

    hdr = "/".join([header_file[1], header_file[0]])
    hf = HeaderFile(hdr)

    # preprocess -- this will create a new file in our temp_dir that
    # was run through cpp and has the name CPP-filename
    cpp.preprocess(hf, add_name="CPP")

    # now scan the preprocessed headers and find any of our function
    # signatures and output to a new unpreprocessed header

    # open the preprocessed header file -- this is what we'll scan
    try:
        hin = open(hf.cpp_name, "r")
    except IOError:
        sys.exit("Cannot open header {}".format(hf.cpp_name))

    # we'll keep track of the signatures that we need to mangle
    signatures = {}

    line = hin.readline()
    while line:

        # if the line does not start a function signature that
        # matches one of our targets, then we ignore it.
        # Otherwise, we need to capture the function signature
        found = None

        # strip comments
        idx = line.find("//")
        tline = line[:idx]

        for target in list(targets):
            target_match = fortran_binding_re.search(tline)
            if target_match:
                if target == target_match.group(3):
                    found = target
                    print('found target {} in header {}'.format(target, hf.cpp_name))
                    break

        # we found a target function, so capture the entire
        # signature, which may span multiple lines
        if found is not None:
            launch_sig = ""
            sig_end = False
            while not line.strip().endswith(";"):
                launch_sig += line
                line = hin.readline()
            launch_sig += line

            signatures[found] = [launch_sig, targets[found]]

        line = hin.readline()

    hin.close()

    # we've now finished going through the header. Note: there may
    # be more signatures here than we really need, because some may
    # have come in via '#includes' in the preprocessing.


    # Now we'll go back to the original file, parse it, making note
    # of any of the signatures we find, but using the preprocessed
    # version in the final output.

    # open the CUDA header for output
    _, tail = os.path.split(hf.name)
    ofile = os.path.join(outdir, tail)
    try:
        hout = open(ofile, "w")
    except IOError:
        sys.exit("Cannot open output file {}".format(ofile))

    # and back to the original file (not preprocessed) for the input
    try:
        hin = open(hf.name, "r")
    except IOError:
        sys.exit("Cannot open output file {}".format(ofile))

    signatures_needed = {}

    line = hin.readline()
    while line:

        # if the line does not start a function signature that we
        # need, then we ignore it
        found = None

        # strip comments
        idx = line.find("//")
        tline = line[:idx]

        # if the line is not a function signature that we already
        # captured then we just write it out
        for target in list(signatures):

            target_match = fortran_binding_re.search(tline)
            if target_match:
                if target == target_match.group(3):
                    found = target
                    signatures_needed[found] = signatures[found]

                    print('found target {} in unprocessed header {}'.format(target, hf.name))
                    break

        if found is not None:
            hout.write(line)
            launch_sig = "" + line
            sig_end = False
            while not sig_end:
                line = hin.readline()
                hout.write(line)
                launch_sig += line
                if line.strip().endswith(";"):
                    sig_end = True

        else:
            # this was not one of our device headers
            hout.write(line)

        line = hin.readline()

    # we are done with the pass through the header and we know all
    # of the signatures that need to be CUDAed

    # remove any dupes in the signatures needed
    signatures_needed = list(set(signatures_needed))

    # now do the CUDA signatures
    hout.write("\n")
    hout.write("#include <AMReX_ArrayLim.H>\n")
    hout.write("#include <AMReX_BLFort.H>\n")
    hout.write("#include <AMReX_CudaDevice.H>\n")
    hout.write("\n")

    hdrmh = os.path.basename(hf.name).strip(".H")

    # Add an include guard -- do we still need this?
    hout.write("#ifndef _cuda_" + hdrmh + "_\n")
    hout.write("#define _cuda_" + hdrmh + "_\n\n")

    # Wrap the device declarations in extern "C"
    hout.write("#ifdef AMREX_USE_GPU_PRAGMA\n")
    hout.write("extern \"C\" {\n\n")

    for name in list(signatures_needed):

        func_sig = signatures[name][0]

        # First write out the device signature
        device_sig = "__device__ {};\n\n".format(func_sig)

        idx = func_sig.find(name)

        # here's the case-sensitive name
        case_name = func_sig[idx:idx+len(name)]

        # Add _device to the function name.

        device_sig = device_sig.replace(case_name, case_name + "_device")

        # Now write out the global signature. This involves
        # getting rid of the data type definitions and also
        # replacing the lo and hi (which must be in the function
        # definition) with blo and bhi.
        dd = decls_re.search(func_sig)
        vars = []

        has_lo = False
        has_hi = False

        intvect_vars = []
        real_vars = []

        for n, v in enumerate(dd.group(3).split(",")):

            # we will assume that our function signatures _always_ include
            # the name of the variable
            _tmp = v.split()
            var = _tmp[-1].replace("*", "").replace("&", "").strip()

            # Replace AMReX Fortran macros
            var = var.replace("BL_FORT_FAB_ARG_3D", "BL_FORT_FAB_VAL_3D")
            var = var.replace("BL_FORT_IFAB_ARG_3D", "BL_FORT_FAB_VAL_3D")

            # Get the list of all arguments which contain each macro.

            args = signatures[name][0].split('(', 1)[1].rsplit(')', 1)[0].split(',')

            for i, arg_positions in enumerate(signatures[name][1]):

                if arg_positions != []:

                    if macro_list[i] == 'AMREX_INT_ANYD':

                        # Replace AMREX_INT_ANYD with the necessary machinery, a set
                        # of three constant ints which will be passed by value.
                        # We only want to do this replacement once, otherwise it will
                        # replace every time for each argument position, so we will
                        # semi-arbitrarily do it in the loop index corresponding to
                        # the actual argument position.

                        for arg_position in arg_positions:

                            if n == arg_position:
                                arg = args[arg_position]
                                v = arg.split()[-1]
                                func_sig = func_sig.replace(arg, "const int {}_1, const int {}_2, const int {}_3".format(v, v, v))
                                device_sig = device_sig.replace(arg, "const int* {}".format(v))
                                intvect_vars.append(v)

                    elif macro_list[i] == 'AMREX_REAL_ANYD':

                        # Same as the above, but with reals.

                        for arg_position in arg_positions:

                            if n == arg_position:
                                arg = args[arg_position]
                                v = arg.split()[-1]
                                func_sig = func_sig.replace(arg, "const amrex::Real {}_1, const amrex::Real {}_2, const amrex::Real {}_3".format(v, v, v))
                                device_sig = device_sig.replace(arg, "const amrex::Real* {}".format(v))
                                real_vars.append(v)

                    elif macro_list[i] == 'BL_TO_FORTRAN_ANYD' or macro_list[i] == 'BL_TO_FORTRAN_N_ANYD' or macro_list[i] == 'BL_TO_FORTRAN_FAB':

                        # Treat this as a real* followed by two copies of AMREX_INT_ANYD,
                        # corresponding to the lo and hi bounds of the box. BL_TO_FORTRAN_FAB
                        # has a fourth argument corresponding to the number of components,
                        # which can be ignored for this logic.

                        for arg_position in arg_positions:
                            if n == arg_position:
                                for pos in [1, 2]:
                                    arg = args[arg_position+pos]
                                    v = arg.split()[-1]
                                    func_sig = func_sig.replace(arg, "const int {}_1, const int {}_2, const int {}_3".format(v, v, v))
                                    device_sig = device_sig.replace(arg, "const int* {}".format(v))
                                    intvect_vars.append(v)

                    elif macro_list[i] == 'BL_TO_FORTRAN_BOX':

                        # This is two copies of AMREX_INT_ANYD.

                        for arg_position in arg_positions:
                            if n == arg_position:
                                for pos in [0, 1]:
                                    arg = args[arg_position+pos]
                                    v = arg.split()[-1]
                                    func_sig = func_sig.replace(arg, "const int {}_1, const int {}_2, const int {}_3".format(v, v, v))
                                    device_sig = device_sig.replace(arg, "const int* {}".format(v))
                                    intvect_vars.append(v)

            if var == "lo":
                var = "blo"
                has_lo = True

            elif var == "hi":
                var = "bhi"
                has_hi = True

            vars.append(var)

        if not has_lo or not has_hi:
            sys.exit("ERROR: function signature must have variables lo and hi defined:\n--- function name:\n {} \n--- function signature:\n {}\n---".format(name, func_sig))

        # reassemble the function sig
        all_vars = ", ".join(vars)
        new_call = "{}({})".format(case_name + "_device", all_vars)

        # Collate all the IntVects that we are going to make
        # local copies of.

        intvects = ""

        if len(intvect_vars) > 0:
            for intvect in intvect_vars:
                intvects += "   int {}[3] = {{{}_1, {}_2, {}_3}};\n".format(intvect, intvect, intvect, intvect)

        # Same for reals.

        reals = ""

        if len(real_vars) > 0:
            for real in real_vars:
                reals += "   amrex::Real {}[3] = {{{}_1, {}_2, {}_3}};\n".format(real, real, real, real)


        hout.write(device_sig)
        hout.write(TEMPLATE.format(func_sig[idx:].replace(';',''), intvects, reals, new_call))
        hout.write("\n")


    # Close out the extern "C" region
    hout.write("\n}\n")
    hout.write("#endif\n")

    # Close out the include guard
    hout.write("\n")
    hout.write("#endif\n")

    hin.close()
    hout.close()


def convert_cxx(inputs):
    """look through the C++ files for "#pragma gpu" and switch it
    to the appropriate CUDA launch macro"""

    cxx_file = inputs[0]
    outdir   = inputs[1]
    cpp      = inputs[2]
    defines  = inputs[3]

    print('looking in C++ file: {}'.format(cxx_file))

    cxx = "/".join([cxx_file[1], cxx_file[0]])

    # open the original C++ file
    try:
        hin = open(cxx, "r")
    except IOError:
        sys.exit("Cannot open header {}".format(cxx))

    # open the C++ file for output
    _, tail = os.path.split(cxx)
    ofile = os.path.join(outdir, tail)
    try:
        hout = open(ofile, "w")
    except IOError:
        sys.exit("Cannot open output file {}".format(ofile))

    print("Working on file {}".format(ofile))

    # look for the appropriate pragma, and once found, capture the
    # function call following it
    line = hin.readline()
    while line:

        # if the line starts with "#pragma gpu", then we need
        # to take action
        if line.startswith("#pragma gpu"):
            # capture pragma options
            # they should be in the format
            # #pragma gpu box(bx) smem(bytes)
            # where
            # box is the box to loop over
            # smem is the number of bytes of dynamic shared memory to reserve

            split_line = line.split()

            box = None
            for entry in split_line:
                if "box(" in entry:
                    box = entry[len("box("):-1]

            smem = 0
            for entry in split_line:
                if "smem(" in entry:
                    smem = entry[len("smem("):-1]

            # we don't need to reproduce the pragma line in the
            # output, but we need to capture the whole function
            # call that follows
            func_call = ""
            line = hin.readline()
            while not line.strip().endswith(";"):
                func_call += line
                line = hin.readline()
            # last line -- remove the semi-colon
            func_call += line.rstrip()[:-1]

            # now split it into the function name and the
            # arguments
            dd = decls_re.search(func_call)
            func_name = dd.group(1).strip().replace(" ", "")
            args = dd.group(3)

            # Convert BL_TO_FORTRAN* to a form that we know will copy in
            # the box indices by value. This is necessary to be in sync
            # with what we did above to the headers.
            args = args.replace('BL_TO_FORTRAN', 'BL_TO_FORTRAN_GPU')

            # Finally output the code in the form we want, with
            # the device launch. We'll also have an option to
            # drop back to a host launch; this is primarily for
            # debugging at this time, but could be used later
            # for dividing the work between the host and device.

            hout.write("#if defined(__CUDA_ARCH__)\n")

            # For the device launch, we need to replace certain macros.
            host_args = args
            host_args = host_args.replace("AMREX_INT_ANYD", "AMREX_ARLIM_3D")
            host_args = host_args.replace("AMREX_REAL_ANYD", "AMREX_ZFILL")
            host_args = host_args.replace("BL_TO_FORTRAN_GPU", "BL_TO_FORTRAN")

            hout.write("{}_device\n ({});\n".format(func_name, host_args))

            hout.write("#else\n")

            hout.write("if (amrex::Gpu::inLaunchRegion()) {\n")
            hout.write("    dim3 {}numBlocks, {}numThreads;\n".format(func_name, func_name))
            if box:
                hout.write("    amrex::Cuda::Device::box_threads_and_blocks({}, {}numBlocks, {}numThreads);\n".format(box, func_name, func_name))
            else:
                hout.write("    amrex::Cuda::Device::grid_stride_threads_and_blocks({}numBlocks, {}numThreads);\n".format(func_name, func_name))
            hout.write("#if ((__CUDACC_VER_MAJOR__ > 9) || (__CUDACC_VER_MAJOR__ == 9 && __CUDACC_VER_MINOR__ >= 1))\n" \
                       "    AMREX_GPU_SAFE_CALL(cudaFuncSetAttribute(&cuda_{}, cudaFuncAttributePreferredSharedMemoryCarveout, 0));\n" \
                       "#endif\n".format(func_name))
            hout.write("    cuda_{}<<<{}numBlocks, {}numThreads, {}, amrex::Cuda::Device::cudaStream()>>>\n    ({});\n".format(func_name, func_name, func_name, smem, args))

            # Catch errors in the launch configuration.

            hout.write("    AMREX_GPU_SAFE_CALL(cudaGetLastError());\n")

            if 'AMREX_DEBUG' in defines:
                hout.write("AMREX_GPU_SAFE_CALL(cudaDeviceSynchronize());\n")

            # For the host launch, we need to replace certain macros.
            host_args = args
            host_args = host_args.replace("AMREX_INT_ANYD", "AMREX_ARLIM_3D")
            host_args = host_args.replace("AMREX_REAL_ANYD", "AMREX_ZFILL")
            host_args = host_args.replace("BL_TO_FORTRAN_GPU", "BL_TO_FORTRAN")

            hout.write("} else {\n")
            hout.write("    {}\n ({});\n".format(func_name, host_args))
            hout.write("}\n")

            hout.write("#endif\n")


        else:
            # we didn't find a pragma
            hout.write(line)

        line = hin.readline()

    hout.close()

    cpp.preprocess(CppFile(ofile), add_name="CPP")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--vpath",
                        help="the VPATH to search for files")
    parser.add_argument("--headers",
                        help="the names of the header files to convert")
    parser.add_argument("--cxx",
                        help="the names of the C++ files to process pragmas")
    parser.add_argument("--output_dir",
                        help="where to write the new header files",
                        default="")
    parser.add_argument("--cpp",
                        help="command to run C preprocessor.  If omitted, then no preprocessing is done",
                        default="")
    parser.add_argument("--defines",
                        help="defines to send to preprocess the files",
                        default="")
    parser.add_argument("--exclude_defines",
                        help="space separated string of directives to remove from defines",
                        default="")
    parser.add_argument("--num_workers",
                        help="number of parallel workers",
                        default="1")


    args = parser.parse_args()

    defines = args.defines

    if args.exclude_defines != "":
        excludes = args.exclude_defines.split()
        for ex in excludes:
            defines = defines.replace(ex, "")

    print("defines: ", defines)

    if args.cpp != "":
        cpp_pass = preprocess.Preprocessor(temp_dir=args.output_dir, cpp_cmd=args.cpp,
                                           defines=defines)
    else:
        cpp_pass = None

    headers, _ = ffv.find_files(args.vpath, args.headers)
    cxx, _ = ffv.find_files(args.vpath, args.cxx)

    num_workers = int(args.num_workers)
    pool = Pool(num_workers)

    # part I: we need to find the names of the Fortran routines that
    # are called from C++ so we can modify the header in the
    # corresponding *_F.H file.

    # A list of specific macros that we want to look for in each target.

    macro_list = ['AMREX_INT_ANYD', 'AMREX_REAL_ANYD', 'BL_TO_FORTRAN_ANYD', 'BL_TO_FORTRAN_N_ANYD', 'BL_TO_FORTRAN_BOX', 'BL_TO_FORTRAN_FAB']

    # look through the C++ for routines launched with #pragma gpu
    targets = find_targets_from_pragmas(args.output_dir, cxx, macro_list, cpp_pass)

    # copy the headers to the output directory, replacing the
    # signatures of the target Fortran routines with the CUDA pair
    inputs = [[header, args.output_dir, targets, macro_list, cpp_pass] for header in headers]
    pool.map(convert_headers, inputs)

    # part II: for each C++ file, we need to expand the `#pragma gpu`
    inputs = [[cxx_file, args.output_dir, cpp_pass, defines] for cxx_file in cxx]
    pool.map(convert_cxx, inputs)

    pool.close()
    pool.join()
