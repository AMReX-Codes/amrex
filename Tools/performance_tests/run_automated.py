#!/usr/common/software/python/2.7-anaconda-4.4/bin/python

import os, sys, shutil, datetime
import argparse, re, time
import pandas as pd
from functions_perftest import *

# typical use: python run_automated.py --n_node_list='1,8,16,32' --automated

parser = argparse.ArgumentParser( description='Run performance tests and write results in files' )
parser.add_argument('--recompile'   , dest='recompile'  , action='store_true' , default=False)
parser.add_argument('--no-recompile', dest='recompile'  , action='store_false', default=False)
parser.add_argument('--commit'      , dest='commit'     , action='store_true' , default=False)
parser.add_argument('--automated'   , dest='automated'  , action='store_true' , default=False  , help='Use to run the automated test list')
parser.add_argument('--n_node_list' , dest='n_node_list',                       default=[]     , help='list ofnumber of nodes for the runs', type=str)
parser.add_argument('--start_date'  , dest='start_date' )
parser.add_argument('--compiler'    , choices=['gnu', 'intel'],                 default='intel', help='which compiler to use')
parser.add_argument('--architecture', choices=['cpu', 'knl']  ,                 default='knl'  , help='which architecture to cross-compile for NERSC machines')
parser.add_argument('--mode'        , choices=['run', 'read'] ,                 default='run'  , help='whether to run perftests or read their perf output. run calls read')

args = parser.parse_args()
do_commit = args.commit
n_node_list_string   = args.n_node_list.split(',')
n_node_list = [int(i) for i in n_node_list_string]
start_date = args.start_date

if args.automated == True:
    do_commit = False
    run_name = 'automated_tests'
    perf_database_file = 'automated_tests_database.h5'
    rename_archive = False
    csv_type = 'alldata' # 'alldata' or 'splitdata'
    write_csv = True
    store_full_input = False
    move_csv = False
    commit_csv = False
    pull_3_repos = False

test_list = []
n_repeat = 2
# run_name, n_nodes, MPI per node, OMP_NUM_THREADS
# n_nodes from this list is not used.
test_list.extend([['automated_test_1_uniform_rest_32ppc', 0, 8, 8]]*n_repeat)
test_list.extend([['automated_test_2_uniform_rest_1ppc',  0, 8, 8]]*n_repeat)
test_list.extend([['automated_test_3_uniform_drift_4ppc', 0, 8, 8]]*n_repeat)
test_list.extend([['automated_test_4_labdiags_2ppc',      0, 8, 8]]*n_repeat)
test_list.extend([['automated_test_5_loadimbalance',      0, 8, 8]]*n_repeat)
test_list.extend([['automated_test_6_output_2ppc',        0, 8, 8]]*n_repeat)     
n_tests = len(test_list)
ncell_dict = {'automated_test_1_uniform_rest_32ppc': [128, 128, 128],
              'automated_test_2_uniform_rest_1ppc' : [256, 256, 512],
              'automated_test_3_uniform_drift_4ppc': [128, 128, 128],
              'automated_test_4_labdiags_2ppc'     : [ 64,  64, 128],
              'automated_test_5_loadimbalance'     : [128, 128, 128],
              'automated_test_6_output_2ppc'       : [128, 256, 256]}
nstep_dict = {'automated_test_1_uniform_rest_32ppc': 10,
              'automated_test_2_uniform_rest_1ppc' : 10,
              'automated_test_3_uniform_drift_4ppc': 10,
              'automated_test_4_labdiags_2ppc'     : 50,
              'automated_test_5_loadimbalance'     : 10,
              'automated_test_6_output_2ppc'       : 1}

# Dictionaries
compiler_name = {'intel': 'intel', 'gnu': 'gcc'}
module_name = {'cpu': 'haswell', 'knl': 'mic-knl'}
module_Cname = {'cpu': 'haswell', 'knl': 'knl,quad,cache'}
cwd = os.getcwd() + '/'
res_dir_base = os.environ['SCRATCH'] + '/performance_warpx/'
bin_dir = cwd + 'Bin/'
bin_name = 'perf_tests3d.' + args.compiler + '.' + module_name[args.architecture] + '.TPROF.MPI.OMP.ex'
log_dir  = cwd
perf_database_file = cwd + perf_database_file

day = time.strftime('%d')
month = time.strftime('%m')
year = time.strftime('%Y')

# Initialize tests
# ----------------
if args.mode == 'run':
    start_date = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    # Set default options for compilation and execution
    config_command = ''
    config_command += 'module unload darshan;' 
    config_command += 'module load craype-hugepages4M;'
    if args.architecture == 'knl':
        if args.compiler == 'intel':
            config_command += 'module unload PrgEnv-gnu;'
            config_command += 'module load PrgEnv-intel;'
        elif args.compiler == 'gnu':
            config_command += 'module unload PrgEnv-intel;'
            config_command += 'module load PrgEnv-gnu;'
        config_command += 'module unload craype-haswell;'
        config_command += 'module load craype-mic-knl;'
    elif args.architecture == 'cpu':
        if args.compiler == 'intel':
            config_command += 'module unload PrgEnv-gnu;'
            config_command += 'module load PrgEnv-intel;'
        elif args.compiler == 'gnu':
            config_command += 'module unload PrgEnv-intel;'
            config_command += 'module load PrgEnv-gnu;'
        config_command += 'module unload craype-mic-knl;'
        config_command += 'module load craype-haswell;'
    # Create main result directory if does not exist
    if not os.path.exists(res_dir_base):
        os.mkdir(res_dir_base)    

# Recompile if requested
if args.recompile == True:
    with open(cwd + 'GNUmakefile_perftest') as makefile_handler:
        makefile_text = makefile_handler.read()
    makefile_text = re.sub('\nCOMP.*', '\nCOMP=%s' %compiler_name[args.compiler], makefile_text)
    with open(cwd + 'GNUmakefile_perftest', 'w') as makefile_handler:
        makefile_handler.write( makefile_text )
    os.system(config_command + " make -f GNUmakefile_perftest realclean ; " + " rm -r tmp_build_dir *.mod; make -j 8 -f GNUmakefile_perftest")
    if os.path.exists( cwd + 'store_git_hashes.txt' ):
        os.remove( cwd + 'store_git_hashes.txt' )
    store_git_hash(repo_path=cwd + '../../../picsar/', filename=cwd + 'store_git_hashes.txt', name='picsar')
    store_git_hash(repo_path=cwd + '../../../amrex/' , filename=cwd + 'store_git_hashes.txt', name='amrex' )
    store_git_hash(repo_path=cwd + '../../'          , filename=cwd + 'store_git_hashes.txt', name='warpx' )

# This function runs a batch script with dependencies to perform the analysis 
# when performance runs are done.
def process_analysis():
    dependencies = ''
    f_log = open(cwd + 'log_jobids_tmp.txt' ,'r')
    line = f_log.readline()
    dependencies += line.split()[3] + ':'
    batch_string = ''
    batch_string += '#!/bin/bash\n'
    batch_string += '#SBATCH --job-name=warpx_1node_read\n'
    batch_string += '#SBATCH --time=00:07:00\n'
#     batch_string += '#SBATCH -C haswell\n'
    batch_string += '#SBATCH -C knl\n'
    batch_string += '#SBATCH -N 1\n'
    batch_string += '#SBATCH -S 4\n'
    batch_string += '#SBATCH -q regular\n'
    batch_string += '#SBATCH -e read_error.txt\n'
    batch_string += '#SBATCH -o read_output.txt\n'
    batch_string += '#SBATCH --mail-type=end\n'
    batch_string += '#SBATCH --account=m2852\n'
    batch_string += 'python ' + __file__ + ' --no-recompile --compiler=' + \
                    args.compiler + ' --architecture=' + args.architecture + \
                    ' --mode=read' + \
                ' --n_node_list=' + '"' + args.n_node_list + '"' + \
                ' --start_date=' + start_date
    if do_commit == True:
        batch_string += ' --commit'
    if args.automated == True:
        batch_string += ' --automated'
    batch_string += '\n'
    batch_file = 'slurm_perfread'
    f_exe = open(batch_file,'w')
    f_exe.write(batch_string)
    f_exe.close()
    os.system('chmod 700 ' + batch_file)
    os.system('sbatch  --dependency afterok:' + dependencies[0:-1] + ' ' + batch_file)
    return 0

# Loop over the tests and return run time + details
# -------------------------------------------------
if args.mode == 'run':
    if os.path.exists( 'log_jobids_tmp.txt' ):
        os.remove( 'log_jobids_tmp.txt' )
    for n_node in n_node_list:
        res_dir = res_dir_base
        res_dir += '_'.join([run_name, args.compiler, args.architecture, str(n_node)]) + '/'
        runtime_param_list = []
        for count, current_run in enumerate( test_list ):
            ncell_scaled_list = scale_ncell(ncell_dict[current_run[0]], n_node)
            max_step = str(nstep_dict[ current_run[0] ])
            runtime_param_string  = ' amr.n_cell=' + ' '.join(str(i) for i in ncell_scaled_list)
            runtime_param_string += ' max_step=' + max_step
            runtime_param_list.append( runtime_param_string )
        # Run the simulations.
        run_batch_nnode(test_list, res_dir, bin_name, config_command,\
                        architecture=args.architecture, Cname=module_Cname[args.architecture], \
                        n_node=n_node, runtime_param_list=runtime_param_list)
    os.chdir(cwd)
    process_analysis()
for n_node in n_node_list:
    print(n_node)
    if args.mode == 'read':
        for count, current_run in enumerate(test_list):
            print('read ' + str(current_run))
            input_file = current_run[0]
            n_mpi    = current_run[2]
            n_omp    = current_run[3]
            n_steps  = nstep_dict[ input_file ]# get_nsteps(cwd  + input_file)
            res_dir = res_dir_base
            res_dir += '_'.join([run_name, args.compiler,\
                                 args.architecture, str(n_node)]) + '/'
            # Read performance data from the output file
            output_filename = 'out_' + '_'.join([input_file, str(n_node), str(n_mpi), str(n_omp), str(count)]) + '.txt'
            # Read data for all test to put in hdf5 a database
            # ------------------------------------------------
            # This is an hdf5 file containing ALL the simulation parameters and results. Might be too large for a repo
            df_newline = extract_dataframe(res_dir + output_filename, n_steps)
            # Add all simulation parameters to the dataframe
            df_newline['git_hashes'] = get_file_content(filename=cwd+'store_git_hashes.txt')
            df_newline['start_date'] = start_date
            df_newline['run_name'] = run_name
            df_newline['input_file'] = input_file
            df_newline['n_node'] = n_node
            df_newline['n_mpi_per_node'] = n_mpi
            df_newline['n_omp'] = n_omp
            df_newline['n_steps'] = n_steps
            df_newline['rep'] = count%n_repeat
            df_newline['date'] = datetime.datetime.now()
            if store_full_input:
                df_newline['inputs_content'] = get_file_content( filename=cwd+input_file )
            if os.path.exists(perf_database_file):
                df_base = pd.read_hdf(perf_database_file, 'all_data')
                updated_df = df_base.append(df_newline, ignore_index=True)
            else:
                updated_df = df_newline
            updated_df.to_hdf(perf_database_file, key='all_data', mode='w')

        # Extract small data from data frame and write them to 
        if write_csv == True:
            # First, generate csv files
            df = pd.read_hdf( perf_database_file )
            # One large file
            if csv_type == 'alldata':
                df.loc[:,'step_time'] = pd.Series(df['time_running']/df['n_steps'], index=df.index)
                df.loc[:, ['date', 'input_file', 'git_hashes', 'n_node', 'n_mpi_per_node', 'n_omp', 'rep', 'start_date', 'time_initialization', 'step_time'] ].to_csv( 'cori_knl.csv' )      
            # Split full data frame into one dataframe per nnode/test, and then write to csv.
            # Errors may occur depending on the version of pandas. I had errors with v0.21.0 solved with 0.23.0
            if csv_type == 'splitdata':
                # one csv file per nnode
                for n_node in pd.unique(df['n_node']):
                    df_small = df[ df['n_node'] == n_node ]
                    df_small.loc[:,'step_time'] = pd.Series(df_small['time_running']/df_small['n_steps'], index=df_small.index)
                    df_small = df_small.loc[:, ['date', 'input_file', 'git_hashes', 'n_node', 'n_mpi', 'rep', 'n_omp', 'start_date', 'time_initialization', 'step_time'] ]
                    agg_df_small = df_small.groupby(['input_file', 'n_node', 'start_date']).aggregate('mean')
                    metadata = df_small[df_small['rep']==0][['date', 'git_hashes', 'input_file', 'n_node', 'start_date']].copy()
                    df_to_write = pd.merge(metadata, agg_df_small, how='left', on=['input_file', 'n_node', 'start_date']).drop('rep', axis=1)
                    df_to_write.to_csv( 'cori_knl_alltests_nnode_' + str(n_node) + '.csv' )

                # one csv file per test
                for test_name in pd.unique(df['input_file']):
                    df_small = df[ df['input_file'] == test_name ]
                    df_small.loc[:,'step_time'] = pd.Series(df_small['time_running']/df_small['n_steps'], index=df_small.index)
                    df_small = df_small.loc[:, ['date', 'input_file', 'git_hashes', 'n_node', 'n_mpi', 'rep', 'n_omp', 'start_date', 'time_initialization', 'step_time'] ]
                    agg_df_small = df_small.groupby(['input_file', 'n_node', 'start_date']).aggregate('mean')
                    metadata = df_small[df_small['rep']==0][['date', 'git_hashes', 'input_file', 'n_node', 'start_date']].copy()
                    df_to_write = pd.merge(metadata, agg_df_small, how='left', on=['input_file', 'n_node', 'start_date']).drop('rep', axis=1)
                    df_to_write.to_csv( 'cori_knl_allnnodes_' + test_name + '.csv' )
            # Second, move files to perf_logs repo
            if move_csv:
                if csv_type == 'alldata':
                    shutil.move( 'cori_knl.csv', os.environ('PERF_LOGS_REPO') + '/logs_csv/cori_knl.csv' )
        # Rename directory with precise date for archive purpose
        if rename_archive == True:
            loc_counter = 0
            res_dir_arch = res_dir_base
            res_dir_arch += '_'.join([year, month, day, run_name, args.compiler,\
                                      args.architecture, str(n_node), str(loc_counter)]) + '/'
            while os.path.exists( res_dir_arch ):
                loc_counter += 1
                res_dir_arch = res_dir_base
                res_dir_arch += '_'.join([year, month, day, run_name, args.compiler,\
                                          args.architecture, str(n_node), str(loc_counter)]) + '/'
            os.rename( res_dir, res_dir_arch )
