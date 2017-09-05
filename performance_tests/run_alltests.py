import os, sys, shutil
import argparse, re

# Read command-line arguments
# ---------------------------

# Create parser and read arguments
parser = argparse.ArgumentParser(
    description='Run performance tests and write results in files')
parser.add_argument('--recompile', dest='recompile', action='store_true')
parser.add_argument('--no-recompile', dest='recompile', action='store_false')
parser.add_argument( '--compiler', choices=['gcc', 'intel'], default='gcc',
    help='which compiler to use')
parser.add_argument( '--architecture', choices=['cpu', 'knl'], default='cpu',
    help='which architecture to cross-compile for NERSC machines')
args = parser.parse_args()

# Define environment variables
res_dir_base = os.environ['SCRATCH'] + '/performance_warpx/'
bin_dir = '../Bin/'
bin_name = 'main3d.' + args.compiler + '.TPROF.MPI.OMP.ex'

# Initialize tests
# ----------------

# Test list
test_list = ['uniform_plasma']
n_tests   = len(test_list)

# Set default options for compilation and run
os.system('module unload darshan')
os.system('module load craype-hugepages4M')
os.system('module load python/2.7-anaconda')
os.system('module load h5py-parallel')
if args.architecture == 'knl':
    if args.compiler == 'intel':
        os.system('module unload PrgEnv-gnu')
        os.system('module load PrgEnv-intel')
    elif args.compiler == 'gcc':
        os.system('module unload PrgEnv-intel')
        os.system('module load PrgEnv-gnu')
    os.system('module unload craype-haswell')
    os.system('module load craype-mic-knl')
elif args.architecture == 'haswell':
    if args.compiler == 'intel':
        os.system('module unload PrgEnv-gnu')
        os.system('module load PrgEnv-intel')
    elif args.compiler == 'gcc':
        os.system('module unload PrgEnv-intel')
        os.system('module load PrgEnv-gnu')
    os.system('module unload craype-mic-knl')
    os.system('module load craype-haswell')

# Recompile if requested
if args.recompile == True:
    os.chdir(build_dir)
    with open('GNUmakefile') as makefile_handler:
        makefile_text = makefile_handler.read()
    makefile_text = re.sub('\nCOMP.*', '\nCOMP=%s' %args.compiler, makefile_text)
    os.system("make")
    
# Create main result directory if does not exist
if not os.path.exists(result_dir):
    os.mkdir(result_dir)
        
# Define function to run a test
# -----------------------------
def runcase(run_name, res_dir, n_node=1, n_mpi=1, n_omp=1):
    # Run a simulation inside an allocation
    # Clean res_dir
    if os.path.exists(res_dir):
        shutil.rmtree(res_dir)
    os.makedirs(res_dir)
    # Copy files to res_dir
    shutil.copyfile(bin_dir + bin_name, res_dir)
    shutil.copyfile(run_name, res_dir)    
    os.chdir(res_dir)
    os.environ('OMP_NUM_THREADS') = n_omp
    # number of logical cores per MPI process
    cflag_value = (68/n_mpi) * 4 
    exec_command = 'srun --cpu_bind=cores ' +  
                    ' -n ' + str(n_nodes*n_mpi +
                    ' -c ' + str(cflag_value)  +
                    ' ./'  + bin_name + ' ' + run_name + 
                    ' > my_output.txt'
    print(os.getcwd())
    print(os.listdir())
    print(exec_command)
    # TO KEEP FOR LATER
    # os.system(exec_command)    

# Loop over the tests and return run time + details
# -------------------------------------------------

for count, run_name in enumerate(test_list):
    runcase(run_name, res_dir + run_name, n_node=1, n_mpi=1, n_omp=1)
    
# Write data in a file

