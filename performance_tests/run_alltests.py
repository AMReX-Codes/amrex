import os, sys, shutil
import argparse, re, time

# Read command-line arguments
# ---------------------------

# Create parser and read arguments
parser = argparse.ArgumentParser(
    description='Run performance tests and write results in files')
parser.add_argument('--recompile', dest='recompile', action='store_true')
parser.add_argument('--no-recompile', dest='recompile', action='store_false')
parser.add_argument( '--compiler', choices=['gnu', 'intel'], default='gnu',
    help='which compiler to use')
parser.add_argument( '--architecture', choices=['cpu', 'knl'], default='cpu',
    help='which architecture to cross-compile for NERSC machines')
args = parser.parse_args()

# Dictionary for compiler names. Use for WarpX executable name
compiler_name = {'intel': 'intel', 'gnu': 'gcc'}

# Define environment variables
res_dir_base = os.environ['SCRATCH'] + '/performance_warpx/'
bin_dir = '../Bin/'
bin_name = 'main3d.' + args.compiler + '.TPROF.MPI.OMP.ex'
log_file = 'performance_log.txt'
log_dir  = os.environ['HOME'] + 'cori/warpx/performance_tests/'

# Initialize tests
# ----------------

# Test list
test_list = ['uniform_plasma']
nsteps_list = [100]
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
    makefile_text = re.sub('\nCOMP.*', '\nCOMP=%s' %compiler_name[args.compiler], makefile_text)
    os.system("make")
    
# Create main result directory if does not exist
if not os.path.exists(res_dir_base):
    os.mkdir(res_dir_base)
    
# Create log_file for performance tests if does not exist
if not os.path.isfile(res_dir_base + log_file):
    log_line = 'run_name year month day time_initialization time_one_iteration\n'
    f_log = open(res_dir_base + log_file, 'a')
    f_log.write(log_line)
    f_log.close()

# Create log_folder for performance tests if does not exist
if not os.path.exists(log_dir):
    os.mkdir(log_dir)
        
# Define function to run a test
# -----------------------------
def runcase(run_name, res_dir, n_node=1, n_mpi=1, n_omp=1):
    # Run a simulation inside an allocation
    print('Running test ' + run_name + ' ...')
    # Clean res_dir
    if os.path.exists(res_dir):
        shutil.rmtree(res_dir)
    os.makedirs(res_dir)
    # Copy files to res_dir
    shutil.copyfile(bin_dir + '/' + bin_name, res_dir + '/' + bin_name)
    shutil.copyfile(run_name, res_dir + '/inputs')
    os.chdir(res_dir)
    os.environ['OMP_NUM_THREADS'] = str(n_omp)
    # number of logical cores per MPI process
    cflag_value = (68/n_mpi) * 4 
    #[LOCAL]
    exec_command = 'srun --cpu_bind=cores '     + \ 
                    ' -n ' + str(n_node*n_mpi) + \
                    ' -c ' + str(cflag_value)   + \
                    ' ./'  + bin_name + ' inputs > my_output.txt'
    os.system('chmod 700 ' + bin_name)
    os.system(exec_command)    

# Loop over the tests and return run time + details
# -------------------------------------------------

for count, run_name in enumerate(test_list):
    res_dir = res_dir_base + run_name
    runcase(run_name, res_dir, n_node=1, n_mpi=1, n_omp=1)
    partition_limit = 'NCalls  Incl. Min  Incl. Avg  Incl. Max   Max %'
    with open('my_output.txt') as file_handler:
        makefile_text = file_handler.read()
    # Get total simulation time
    line_match_totaltime = re.search('TinyProfiler total time across processes.*', makefile_text)
    total_time = float(line_match_totaltime.group(0).split()[8])
    search_area = makefile_text.partition(partition_limit)[2]
    line_match_looptime = re.search('\nWarpX::Evolve().*', search_area)
    time_wo_initialization = float(line_match_looptime.group(0).split()[3])
    time_one_iteration = time_wo_initialization/nsteps_list[count]
    time_initialization = total_time - time_wo_initialization
    # Write data in a file
    day = time.strftime('%d') + ' '
    month = time.strftime('%m') + ' '
    year = time.strftime('%Y') + ' '
    log_line = run_name + ' ' + year + ' ' + month + ' ' + day + ' ' + \
       str(time_initialization) + ' ' + str(time_one_iteration) + '\n'

    f_log = open(res_dir_base + log_file, 'a')
    f_log.write(log_line)
    f_log.close()

shutil.copyfile(res_dir_base + log_file, log_dir + year + '_' + month + '_' + day + '_' + log_file)