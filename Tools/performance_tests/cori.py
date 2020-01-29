# Copyright 2019 Axel Huebl, Luca Fedeli, Maxence Thevenet
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import os, copy

from functions_perftest import test_element

module_name = {'cpu': 'haswell.', 'knl': 'mic-knl.', 'gpu':'.'}

def executable_name(compiler, architecture):
    return 'perf_tests3d.' + compiler + \
        '.' + module_name[architecture] + 'TPROF.MPI.OMP.ex'

def get_config_command(compiler, architecture):
    config_command = ''
    config_command += 'module unload darshan;'
    if architecture == 'knl':
        if compiler == 'intel':
            config_command += 'module unload PrgEnv-gnu;'
            config_command += 'module load PrgEnv-intel;'
        elif compiler == 'gnu':
            config_command += 'module unload PrgEnv-intel;'
            config_command += 'module load PrgEnv-gnu;'
        config_command += 'module unload craype-haswell;'
        config_command += 'module load craype-mic-knl;'
    elif architecture == 'cpu':
        if compiler == 'intel':
            config_command += 'module unload PrgEnv-gnu;'
            config_command += 'module load PrgEnv-intel;'
        elif compiler == 'gnu':
            config_command += 'module unload PrgEnv-intel;'
            config_command += 'module load PrgEnv-gnu;'
        config_command += 'module unload craype-mic-knl;'
        config_command += 'module load craype-haswell;'
    return config_command

# This function runs a batch script with
# dependencies to perform the analysis
# after all performance tests are done.
def process_analysis(automated, cwd, compiler, architecture, n_node_list, start_date, path_source, path_results):
    dependencies = ''
    f_log = open(cwd + 'log_jobids_tmp.txt' ,'r')
    for line in f_log.readlines():
        dependencies += line.split()[3] + ':'

    batch_string = '''#!/bin/bash
#SBATCH --job-name=warpx_1node_read
#SBATCH --time=00:07:00
#SBATCH -C knl
#SBATCH -N 1
#SBATCH -S 4
#SBATCH -q regular
#SBATCH -e read_error.txt
#SBATCH -o read_output.txt
#SBATCH --mail-type=end
#SBATCH --account=m2852
module load h5py-parallel
'''
    batch_string += 'python run_automated.py --compiler=' + \
        compiler + ' --architecture=' + architecture + \
        ' --mode=read' + \
        ' --n_node_list=' + '"' + n_node_list + '"' + \
        ' --start_date=' + start_date + \
        ' --path_source=' + path_source + \
        ' --path_results=' + path_results
    if automated == True:
        batch_string += ' --automated'
    batch_string += '\n'
    batch_file = 'slurm_perfread'
    f_exe = open(batch_file,'w')
    f_exe.write(batch_string)
    f_exe.close()
    os.system('chmod 700 ' + batch_file)
    print( 'process_analysis line:  ' + 'sbatch  --dependency afterok:' + dependencies[0:-1] + ' ' + batch_file)
    os.system('sbatch  --dependency afterok:' + dependencies[0:-1] + ' ' + batch_file)

# Calculate simulation time. Take 5 min + 5 min / simulation
def time_min(nb_simulations):
    return 5. + nb_simulations*5.

def get_submit_job_command():
    return ' sbatch '

def get_batch_string(test_list, job_time_min, Cname, n_node):

    job_time_str = str(int(job_time_min/60)) + ':' + str(int(job_time_min%60)) + ':00'

    batch_string = ''
    batch_string += '#!/bin/bash\n'
    batch_string += '#SBATCH --job-name=' + test_list[0].input_file + '\n'
    batch_string += '#SBATCH --time=' + job_time_str + '\n'
    batch_string += '#SBATCH -C ' + Cname + '\n'
    batch_string += '#SBATCH -N ' + str(n_node) + '\n'
    batch_string += '#SBATCH -q regular\n'
    batch_string += '#SBATCH -e error.txt\n'
    batch_string += '#SBATCH --account=m2852\n'
    return batch_string

def get_run_string(current_test, architecture, n_node, count, bin_name, runtime_param_string):
    srun_string = ''
    srun_string += 'export OMP_NUM_THREADS=' + str(current_test.n_omp) + '\n'
    # number of logical cores per MPI process
    if architecture == 'cpu':
        cflag_value = max(1, int(32/current_test.n_mpi_per_node) * 2) # Follow NERSC directives
    elif architecture == 'knl':
        cflag_value = max(1, int(64/current_test.n_mpi_per_node) * 4) # Follow NERSC directives
    output_filename = 'out_' + '_'.join([current_test.input_file, str(n_node), str(current_test.n_mpi_per_node), str(current_test.n_omp), str(count)]) + '.txt'
    srun_string += 'srun --cpu_bind=cores '+ \
        ' -n ' + str(n_node*current_test.n_mpi_per_node) + \
        ' -c ' + str(cflag_value)   + \
        ' ./'  + bin_name + \
        ' ' + current_test.input_file + \
        runtime_param_string + \
        ' > ' + output_filename + '\n'
    return srun_string

def get_test_list(n_repeat):
    test_list_unq = []
    # n_node is kept to None and passed in functions as an external argument
    # That way, several test_element_instance run with the same n_node on the same batch job
    test_list_unq.append( test_element(input_file='automated_test_1_uniform_rest_32ppc',
                                       n_mpi_per_node=8,
                                       n_omp=8,
                                       n_cell=[128, 128, 128],
                                       max_grid_size=64,
                                       blocking_factor=32,
                                       n_step=10) )
    test_list_unq.append( test_element(input_file='automated_test_2_uniform_rest_1ppc',
                                       n_mpi_per_node=8,
                                       n_omp=8,
                                       n_cell=[256, 256, 512],
                                       max_grid_size=64,
                                       blocking_factor=32,
                                       n_step=10) )
    test_list_unq.append( test_element(input_file='automated_test_3_uniform_drift_4ppc',
                                       n_mpi_per_node=8,
                                       n_omp=8,
                                       n_cell=[128, 128, 128],
                                       max_grid_size=64,
                                       blocking_factor=32,
                                       n_step=10) )
    test_list_unq.append( test_element(input_file='automated_test_4_labdiags_2ppc',
                                       n_mpi_per_node=8,
                                       n_omp=8,
                                       n_cell=[64, 64, 128],
                                       max_grid_size=64,
                                       blocking_factor=32,
                                       n_step=50) )
    test_list_unq.append( test_element(input_file='automated_test_5_loadimbalance',
                                       n_mpi_per_node=8,
                                       n_omp=8,
                                       n_cell=[128, 128, 128],
                                       max_grid_size=64,
                                       blocking_factor=32,
                                       n_step=10) )
    test_list_unq.append( test_element(input_file='automated_test_6_output_2ppc',
                                       n_mpi_per_node=8,
                                       n_omp=8,
                                       n_cell=[128, 256, 256],
                                       max_grid_size=64,
                                       blocking_factor=32,
                                       n_step=0) )
    test_list = [copy.deepcopy(item) for item in test_list_unq for _ in range(n_repeat) ]
    return test_list
