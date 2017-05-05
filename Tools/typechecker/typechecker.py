#!/usr/bin/env python

import os
import sys
import argparse
from pycparser import parse_file, c_parser, c_ast, c_generator

# supported C types: char, short, int, long, float, double, and their pointer types, and void*

def typechecker(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir",
                        required=True)
    script_args = parser.parse_args()

    ret_type, arg_type = getFortranArg('mysub', script_args.workdir)
    print ret_type, arg_type

    ret_type, arg_type = getFortranArg('myfunc', script_args.workdir)
    print ret_type, arg_type

    ret_type, arg_type = getFortranArg('myfunc2', script_args.workdir)
    print ret_type, arg_type

def getFortranArg(funcname, workdir):
    """ given a function name and a directory for searching, return a tuple
        of function return type and list of arguments.  For example,
        (INTEGER 4, [('INTERGER 4', 'value'), ('REAL 8', 'pointer')])."""

    # gfortran's internal parse tree look like this
    # symtree: 'fort_fab_copy'|| symbol: 'fort_fab_copy' 
    #   type spec : (UNKNOWN 0)
    #   attributes: (PROCEDURE MODULE-PROC  BIND(C) SUBROUTINE)
    #   Formal arglist: lo hi dst dlo dhi src slo shi sblo ncomp
    # procedure name = fort_fab_copy
    #   symtree: 'dhi'         || symbol: 'dhi'          
    #     type spec : (INTEGER 4)
    #     attributes: (VARIABLE  DIMENSION DUMMY(IN))
    #     Array spec:(1 [0] AS_EXPLICIT 1 3 )
    #   symtree: 'dlo'         || symbol: 'dlo'          
    #     type spec : (INTEGER 4)
    #     attributes: (VARIABLE  DIMENSION DUMMY(IN))
    #     Array spec:(1 [0] AS_EXPLICIT 1 3 )
    #   ...
    #   symtree: 'fort_fab_copy'|| symbol: 'fort_fab_copy' from namespace 'basefab_nd_module'
    #   symtree: 'hi'          || symbol: 'hi'           
    #     type spec : (INTEGER 4)
    #     attributes: (VARIABLE  DIMENSION DUMMY(IN))
    #     Array spec:(1 [0] AS_EXPLICIT 1 3 )
    #   ...
    #   code:
    #   ...

    node_tok = "symtree: '" + funcname + "'"
    args_tok = "Formal arglist:"
    proc_tok = "procedure name = " + funcname + "\n"
    return_type = ''
    arguments_type = []
    for fname in os.listdir(workdir):
        if fname.endswith(('F90.orig', 'f90.orig')):
            f = open(os.path.join(workdir,fname), 'r')
            func_found = False;
            line = ''
            while True:  # try to find the function
                line = f.readline()
                if not line:
                    break
                if node_tok in line:
                    func_found = True
                    line = f.readline() # like, type spec : (UNKNOWN 0)
                                        # or,   type spec : (INTEGER 4)
                    if not 'type spec :' in line:
                        print("Why is this line not type spec? {0}".format(line))
                        sys.exit(1)
                    if "UNKNOWN" in line:
                        return_type = 'void'
                    else:
                        return_type = line.split('(')[1].split(')')[0]
                    break
            if func_found:
                num_white_spaces = len(line) - len(line.lstrip())
                # we now need to find the function's formal arglist
                # we may not find it because the funciton may not have any argument.
                # In that case, we stop searching at the end of this block (indicated
                # by the number of leading white space).
                while True:
                    line = f.readline()
                    current_num_white_spaces = len(line) - len(line.lstrip())
                    if current_num_white_spaces < num_white_spaces:
                        arglist = []
                    if args_tok in line:
                        arglist = line.replace(args_tok, "").split()
                        break
                if len(arglist) > 0:
                    while True: # need to find the procedure
                        line = f.readline()
                        if not line:
                            print("Cannot find procedure" + funcname)
                            sys.exit(1)
                        if proc_tok in line:
                            break
                    func_args = {}
                    while True: # need to build a dictionary of argument and their types
                        line = f.readline()
                        if "code:\n" in line:  # the code section starts after the argument section
                            break
                        if "symtree: '" in line:
                            # the line should look like,    symtree: 'dhi'         || symbol: 'dhi'
                            argname = line.split("'")[1]
                            if argname in arglist:
                                line = f.readline() # like,  type spec : (INTEGER 4)
                                argtype = line.split('(')[1].split(')')[0]
                                line = f.readline() # like,  attributes: (VARIABLE  VALUE DUMMY(IN))
                                if "VALUE" in line:
                                    argattrib = 'value'
                                else:
                                    argattrib = 'pointer'
                                func_args[argname] = (argtype, argattrib)
                        if len(func_args) == len(arglist):
                            break
                    for a in arglist:
                        if a in func_args.keys():
                            arguments_type.append(func_args[a])
            f.close()
    return return_type, arguments_type

if __name__ == "__main__":
    typechecker(sys.argv)
