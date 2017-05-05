#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import argparse
from pycparser import parse_file, c_ast

def typechecker(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--debug", action='store_true')
    script_args = parser.parse_args()

    class FuncDeclVisitor(c_ast.NodeVisitor):
        def visit_FuncDecl(self, node):
            check_doit(node, script_args.workdir, script_args.debug)

    for f in os.listdir(script_args.workdir):
        if f.endswith('-cppd.h'):
            fname = os.path.join(script_args.workdir,f)
            if script_args.debug:
                print("\nChecking "+fname+"...")
            ast = parse_file(fname)
            v = FuncDeclVisitor()
            v.visit(ast)


def check_doit(node, workdir, debug):
    f_ret_type, f_arg_type = getFortranArg(node.type.declname, workdir)
    if f_ret_type:
        # found a fortran function with that name
        c_ret_type = c_ast_get_type(node.type)
        c_arg_type = []
        for param in node.args.params:
            typ = c_ast_get_type(param.type).split(' ')
            if len(typ) == 1:
                c_arg_type.append([typ[0], 'value'])
            elif len(typ) == 2:
                c_arg_type.append(typ)
            else:
                print("check_doit: how did this happen? ", typ)
                sys.exit(1)
        if debug:
            print("\nFunction "+node.type.declname+": C vs. Fortran")
            if c_to_f_type(c_ret_type) == f_ret_type:
                print("    C return type {0} matches F {1}."
                      .format(c_ret_type,f_ret_type))
            else:
                print("    C return type {0} does NOT F match {1}."
                      .format(c_ret_type,f_ret_type))
            if len(c_arg_type) == len(f_arg_type):
                print("    number of arguments {0} matches {1}."
                      .format(len(c_arg_type),len(f_arg_type)))
            else:
                print("    number of arguments {0} does NOT matche {1}."
                      .format(len(c_arg_type),len(f_arg_type)))
            for ct,ft in zip(c_arg_type, f_arg_type):
                if c_to_f_type(ct[0]) == ft[0] and ct[1] == ft[1]:
                    print("    C arg type {0} matches F arg type {1}."
                          .format(ct, ft))
                else:
                    print("    C arg type {0} does NOT match F arg type {1}."
                          .format(ct, ft))


def c_to_f_type(ctyp):
    # supported C types: char, short, int, long, float, double, void
    if ctyp == 'char':
        return 'INTEGER 1'
    elif ctyp == 'short':
        return 'INTEGER 2'
    elif ctyp == 'int':
        return 'INTEGER 4'
    elif ctyp == 'long':
        return 'INTEGER 8'  # yes, this is our assumption
    elif ctyp == 'float':
        return 'REAL 4'
    elif ctyp == 'double':
        return 'REAL 8'
    elif ctyp == 'amrex_real':
        return 'REAL 8'
    elif ctyp == 'void':
        return 'void'


def c_ast_get_type(decl):
    typ = type(decl)
    if typ == c_ast.IdentifierType:
        return ' '.join(decl.names)
    elif typ == c_ast.PtrDecl:
        return c_ast_get_type(decl.type) + ' pointer'
    else:
        return c_ast_get_type(decl.type)
    
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
