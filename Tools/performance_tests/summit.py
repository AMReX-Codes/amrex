import os

def executable_name(compiler,architecture):
    return 'perf_tests3d.' + compiler + '.TPROF.MPI.ACC.CUDA.ex'

def get_config_command(compiler, architecture):
    config_command = ''
    config_command += 'module load pgi;'
    config_command += 'module load cuda;'
    return config_command

# This function runs a batch script with 
# dependencies to perform the analysis 
# after all performance tests are done.
def process_analysis(automated, cwd, compiler, architecture, n_node_list, start_date):

    batch_string = ''
    batch_string += '#!/bin/bash\n'
    batch_string += '#BSUB -P APH114\n'
    batch_string += '#BSUB -W 00:10\n'
    batch_string += '#BSUB -nnodes 1\n'
    batch_string += '#BSUB -J perf_test\n'
    batch_string += '#BSUB -o read_output.txt\n'
    batch_string += '#BSUB -e read_error.txt\n'
    f_log = open(cwd + 'log_jobids_tmp.txt' ,'r')
    for line in f_log.readlines():
        dependency = line.split()[1][1:-1]
        batch_string += '#BSUB -w ended(' + dependency + ')\n'

    batch_string += 'python run_automated.py --compiler=' + \
        compiler + ' --architecture=' + architecture + \
        ' --mode=read' + \
        ' --n_node_list=' + '"' + n_node_list + '"' + \
        ' --start_date=' + start_date
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

def get_batch_string(test_list, job_time_min, Cname, n_node):

    job_time_str = str(int(job_time_min/60)) + ':' + str(int(job_time_min%60))

    batch_string = ''
    batch_string += '#!/bin/bash\n'
    batch_string += '#BSUB -P APH114\n'
    batch_string += '#BSUB -W ' + job_time_str + '\n'
    batch_string += '#BSUB -nnodes ' + str(n_node) + '\n'
    batch_string += '#BSUB -J ' + test_list[0].input_file + '\n'
    batch_string += '#BSUB -e error.txt\n'
    batch_string += 'module load pgi\n' 
    batch_string += 'module load cuda\n' 
    return batch_string

def get_run_string(current_test, architecture, n_node, count, bin_name, runtime_param_list):

    output_filename = 'out_' + '_'.join([current_test.input_file, str(n_node), str(current_test.n_mpi_per_node), str(current_test.n_omp), str(count)]) + '.txt'

    srun_string = ''
    srun_string += 'jsrun '
    srun_string += ' -n ' + str(n_node)
    srun_string += ' -a 6 -g 6 -c 6 --bind=packed:1 '
    srun_string += ' ./' + bin_name + ' '
    srun_string += current_test.input_file + ' '
    srun_string += runtime_param_list[ count ]
    srun_string += ' > ' + output_filename + '\n'
    return srun_string
