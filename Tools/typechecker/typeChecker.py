import os
import sys
from pycparser import parse_file, c_parser, c_ast, c_generator

#These types are commonly used in AMReX. We will add more types into the dictionaries later if needed
c_primitive_types= ["char", "short", "int", "long", "float", "double"]
fortran_primitive_types= ["CHARACTER", "BYTE", "INTEGER 1", "INTEGER 2", "INTEGER 4", "INTEGER 8", "LOGICAL 1", "LOGICAL 2", "LOGICAL 4", "LOGICAL 8", "REAL 4", "REAL 8", "REAL 16", "COMPLEX 8", "COMPLEX 16", "DOUBLE PRECISION"]

def may_type_safe(ctype, cdim, ftype, fdim):
  if cdim==0: 
    return False #we can't pass by value 
  if not ctype in c_primitive_types: return True
  if not ftype in fortran_primitive_types: return True 
  #now cdim always >=1
  if cdim!=fdim: 
    if cdim==1 and fdim==0: return True #passing by reference
    else: return cdim < fdim #in C BoxLib users often employ flattened arrays
  #now cdim ==fdim
  if ctype=="char": return ftype=="CHARACTER" or ftype=="BYTE" or ftype=="INTEGER 1" or ftype=="LOGICAL 1"
  elif ctype=="short": return ftype=="INTEGER 2" or ftype=="LOGICAL 2"
  elif ctype=="int": return ftype=="INTEGER 4" or ftype=="LOGICAL 4"
  elif ctype=="long": return ftype=="INTEGER 4" or ftype=="LOGICAL 4"
  elif ctype=="float": return ftype=="REAL 4" or ftype=="COMPLEX 4"
  elif ctype=="double": return ftype=="DOUBLE PRECISION" or ftype=="REAL 8" or ftype=="COMPLEX 8"

def checkFortranArguments(funcName, functionNode):
  str= "symtree: '"+funcName+"'"
  for file in os.listdir("."):
    if file.endswith(('F90.orig', 'f90.orig')):
      f = open(file, 'r')
      a_line = f.readline()
      arglist = ""
      state=0
      while a_line:
        if str in a_line:
          if(state==0): #first occurence
            f.readline()
            f.readline()
            a_line = f.readline()
            if not "arglist" in a_line:
              a_line = f.readline()
            arglist= a_line.replace("Formal arglist:", "")
            state = state+1
            str= "procedure name = "+funcName
          else: #state==1 --- we can get the type or arguments
            argarray= arglist.split()
            f_argcount= len(argarray)
            if f_argcount != len(functionNode.args.params):
              print "Argument mismatch: %s has %d parameters but received %d arguments" %(funcName, len(functionNode.args.params), f_argcount)
              return False
            a_line = f.readline()
            while a_line:
              if "procedure name" in a_line: break
              if("symtree" in a_line):
                argname = a_line.split()[1] 
                argname = argname[1:]
                argname = argname[:-1]
                if(argname in argarray): #this local variable is an argument
                  #find the type of this argument
                  a_line = f.readline()
                  word_array = a_line.split()
                  f_type = word_array[3]
                  f_type = f_type[1:]
                  f_type_size= word_array[4]
                  f_type_size= f_type_size[:-1]
                  f_type=  f_type+" "+f_type_size
                  f.readline() #skip 1 line
                  f_dimension =0
		  a_line = f.readline()
                  if "Array spec" in a_line: 
		    word_array= a_line.split()[1]
		    f_dimension= word_array[-1:] 
                  c_type =functionNode.args.params[argarray.index(argname)].type
                  if isinstance(c_type, c_ast.PtrDecl) !=True: #XML input
                    c_dimension=functionNode.args.params[argarray.index(argname)].dim
                    if may_type_safe(c_type, c_dimension, f_type, int(f_dimension)) == False:
	              print "Possible type mismatch when calling routine %s with argument number %d" %(funcName, argarray.index(argname))
                  else: #c header input
                    c_dimension=0
                    while isinstance(c_type,c_ast.PtrDecl)==True: 
                      c_dimension = c_dimension+1
                      c_type= c_type.type
                    if may_type_safe(c_type.type.names[0], c_dimension, f_type, int(f_dimension)) == False:
	              print "Possible type mismatch when calling routine %s with argument number %d" %(funcName, argarray.index(argname))
                  continue
              a_line = f.readline()
            state= state+1
        if state==2: break
        a_line = f.readline()
  return True


#Move to the directory which holds parsed trees of C and Fortran headers
os.chdir(sys.argv[1])



#############################################################
class argument():
    def __init__(self, name, type, dim):
        self.name = name
        self.type = type
        self.dim = dim

class ParamList():
    def __init__(self, params):
        self.params = params

class FuncDecl():
    def __init__(self, args, type):
        self.args = args
        self.type = type

#If the C input is the parsed C header files using gccXML (.xml)
for file in os.listdir("."):
  if file.endswith(".xml"):
    wholeFileStr= open(file, 'r').read()
    f = open(file, 'r')
    line = f.readline()
    while line:
      if "Function id=" in line:
        funcName= line.split()[2]
        funcName= funcName[6:-1]
        if(funcName[:10] != "__builtin_"): 
	  #print (funcName)
	  #Construct an abstract syntax tree of a function declaration, given the input has been parsed into an XML form
          argList = []
          #look for types of all arguments
          line = f.readline()
          while line.split()[0] != "</Function>":
	    for word in line.split():
              if word[:5] == "type=": 
                dimension =0
                typeID = word[6:-1]
                while True:
                  keyword = "id=\""+typeID+"\""
	  	  before, keyword, after = wholeFileStr.partition(keyword)
                  tag = before.split()[-1]
                  if tag == "<PointerType" or tag == "<ReferenceType":
                    typeID = after.split()[0]
                    typeID = typeID[6:-1]
		    dimension= dimension+1
                  elif tag == "<Typedef":
                    typeID = after.split()[1]
                    typeID = typeID[6:-1]
                  elif tag=="<CvQualifiedType":
                    typeID = after.split()[0]
                    typeID = typeID[6:-1]
                  elif before.split()[-1] == "<FundamentalType":
	            typeName = after.split()[0]
                    typeName = typeName[6:-1]
                    argList.append(argument("", typeName, dimension))
                    break
                  else:
		    print "Unknown type found"
                    break
            line = f.readline()
          node= FuncDecl(ParamList(argList), "void") 
          checkFortranArguments(funcName, node)
      line = f.readline()
    f.close()


#############################################################
#If the C input is C header files (.h)

class FuncDeclVisitor(c_ast.NodeVisitor):
  def visit_FuncDecl(self, node):
    checkFortranArguments(node.type.declname, node)
    #node.show()

def preprocessCHeaders():
  for file in os.listdir("."):
    if file.endswith(".h"):
      f = open(file, 'rw')
      if os.path.exists(file[:-1]+"i"):
        print("Error")
        break
      else:
        f1 = open(file[:-1]+"i", 'w')
      line = f.readline()
      while line:
        line= line.replace("extern \"C\"","\nextern \"C\"\n")
        line= line.replace("{","\n{\n")
        line= line.replace("}","\n}\n")
        f1.write(line)
        line = f.readline()
  
#create a new file without namespace, extern "C" {} so pycparser can process
preprocessCHeaders()
for file in os.listdir("."):
  if file.endswith(".i"):
    f = open(file, 'r')
    if os.path.exists(file[:-1]+"j"):
      print("Error")
      break
    else:
      f1 = open(file[:-1]+"j", 'w')
    a_line = f.readline()
    state_stack = []
    while a_line:
      if "extern \"C\"" in a_line or "namespace" in a_line:
        state_stack.append(1) #remove extern "C" and namespace
      elif "{" in a_line:
        if(len(state_stack) >0): state= state_stack.pop()
        else: state=0
        if state==1: #remove this bracket
          state_stack.append(1) 
          state_stack.append(2) 
        else: 
          state_stack.append(0) 
          f1.write(a_line)
      elif "}" in a_line:
        state= state_stack.pop()
        if state==2: 
          state_stack.pop()
        else: 
          f1.write(a_line)
      else: 
        if(a_line.strip() != ''): 
          if(len(state_stack) >0): state= state_stack.pop()
          else: state=-1
          if state==1: 
	    pass  #extern C isn't followed by {
          elif state!=-1: state_stack.append(state)
          replaced_line= a_line.replace('&','')  #remove all & occurences
          f1.write(replaced_line)
      a_line = f.readline()
    f.close()
    f1.close()
    os.remove(file)

for file in os.listdir("."):
  if file.endswith(".j"):
    f = open(file, 'r')
    #print ("Processing file "+file)
    ast = parse_file(file, use_cpp=True)
    v = FuncDeclVisitor()
    v.visit(ast)
    f.close()
    os.remove(file)
  


