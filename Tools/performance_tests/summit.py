# Copyright 2019 Axel Huebl, Luca Fedeli, Maxence Thevenet
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# requirements:
# - module load python/3.7.0-anaconda3-5.3.0

import os, copy
from functions_perftest import test_element

def executable_name(compiler,architecture):
    return 'perf_tests3d.' + compiler + '.TPROF.MPI.CUDA.ex'

def get_config_command(compiler, architecture):
    config_command = ''
    config_command += 'module load gcc;'
    config_command += 'module load cuda;'
    return config_command

# This function runs a batch script with
# dependencies to perform the analysis
# after all performance tests are done.
def process_analysis(automated, cwd, compiler, architecture, n_node_list, start_date, path_source, path_results):

    batch_string = '''#!/bin/bash
#BSUB -P APH114
#BSUB -W 00:10
#BSUB -nnodes 1
#BSUB -J perf_test
#BSUB -o read_output.txt
#BSUB -e read_error.txt
'''
    f_log = open(cwd + 'log_jobids_tmp.txt' ,'r')
    for line in f_log.readlines():
        dependency = line.split()[1][1:-1]
        batch_string += '#BSUB -w ended(' + dependency + ')\n'

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
    batch_file = 'bsub_perfread'
    f_exe = open(batch_file,'w')
    f_exe.write(batch_string)
    f_exe.close()
    os.system('chmod 700 ' + batch_file)
    print( 'process_analysis line:  ' + 'bsub ' + batch_file)
    os.system('bsub ' + batch_file)

# Calculate simulation time. Take 2 min + 2 min / simulation
def time_min(nb_simulations):
    return 2. + nb_simulations*2.

def get_submit_job_command():
    return ' bsub '

def get_batch_string(test_list, job_time_min, Cname, n_node):

    job_time_str = str(int(job_time_min/60)) + ':' + str(int(job_time_min%60))

    batch_string = ''
    batch_string += '#!/bin/bash\n'
    batch_string += '#BSUB -P APH114\n'
    batch_string += '#BSUB -W ' + job_time_str + '\n'
    batch_string += '#BSUB -nnodes ' + str(n_node) + '\n'
    batch_string += '#BSUB -J ' + test_list[0].input_file + '\n'
    batch_string += '#BSUB -e error.txt\n'
    batch_string += 'module load gcc\n'
    batch_string += 'module load cuda\n'
    return batch_string

def get_run_string(current_test, architecture, n_node, count, bin_name, runtime_param_string):

    output_filename = 'out_' + '_'.join([current_test.input_file, str(n_node), str(current_test.n_mpi_per_node), str(current_test.n_omp), str(count)]) + '.txt'

    ngpu = str(current_test.n_mpi_per_node)
    srun_string = ''
    srun_string += 'jsrun '
    srun_string += ' -n ' + str(n_node)
    srun_string += ' -a ' + ngpu + ' -g ' + ngpu + ' -c ' + ngpu + ' --bind=packed:1 '
    srun_string += ' ./' + bin_name + ' '
    srun_string += current_test.input_file + ' '
    srun_string += runtime_param_string
    srun_string += ' > ' + output_filename + '\n'
    return srun_string

def get_test_list(n_repeat):
    test_list_unq = []
    # n_node is kept to None and passed in functions as an external argument
    # That way, several test_element_instance run with the same n_node on the same batch job
    test_list_unq.append( test_element(input_file='automated_test_1_uniform_rest_32ppc',
                                       n_mpi_per_node=6,
                                       n_omp=1,
                                       n_cell=[128, 128, 192],
                                       max_grid_size=256,
                                       blocking_factor=32,
                                       n_step=10) )
    test_list_unq.append( test_element(input_file='automated_test_2_uniform_rest_1ppc',
                                       n_mpi_per_node=6,
                                       n_omp=1,
                                       n_cell=[256, 512, 768],
                                       max_grid_size=512,
                                       blocking_factor=256,
                                       n_step=10) )
    test_list_unq.append( test_element(input_file='automated_test_3_uniform_drift_4ppc',
                                       n_mpi_per_node=6,
                                       n_omp=1,
                                       n_cell=[128, 128, 384],
                                       max_grid_size=256,
                                       blocking_factor=64,
                                       n_step=10) )
    test_list_unq.append( test_element(input_file='automated_test_4_labdiags_2ppc',
                                       n_mpi_per_node=6,
                                       n_omp=1,
                                       n_cell=[384, 256, 512],
                                       max_grid_size=256,
                                       blocking_factor=128,
                                       n_step=50) )
    test_list_unq.append( test_element(input_file='automated_test_5_loadimbalance',
                                       n_mpi_per_node=6,
                                       n_omp=1,
                                       n_cell=[64, 64, 192],
                                       max_grid_size=64,
                                       blocking_factor=32,
                                       n_step=10) )
    test_list_unq.append( test_element(input_file='automated_test_6_output_2ppc',
                                       n_mpi_per_node=6,
                                       n_omp=1,
                                       n_cell=[384, 256, 512],
                                       max_grid_size=256,
                                       blocking_factor=64,
                                       n_step=0) )
    test_list = [copy.deepcopy(item) for item in test_list_unq for _ in range(n_repeat) ]
    return test_list
