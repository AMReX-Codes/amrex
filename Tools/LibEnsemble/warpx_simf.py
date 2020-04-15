import os
import time
import numpy as np

from libensemble.executors.executor import Executor
from libensemble.message_numbers import WORKER_DONE, TASK_FAILED
from read_sim_output import read_sim_output
from write_sim_input import write_sim_input

"""
This file is part of the suite of scripts to use LibEnsemble on top of WarpX
simulations. It defines a sim_f function that takes LibEnsemble history and
input parameters, run a WarpX simulation and returns 'f'.
"""


def run_warpx(H, persis_info, sim_specs, libE_info):
    """
    This function runs a WarpX simulation and returns quantity 'f' as well as
    other physical quantities measured in the run for convenience. Status check
    is done periodically on the simulation, provided by LibEnsemble.
    """

    # Setting up variables needed for input and output
    # keys              = variable names
    # x                 = variable values
    # libE_output       = what will be returned to libE

    calc_status = 0  # Returns to worker

    input_file = sim_specs['user']['input_filename']
    time_limit = sim_specs['user']['sim_kill_minutes'] * 60.0
    machine_specs = sim_specs['user']['machine_specs']

    exctr = Executor.executor  # Get Executor

    # Modify WarpX input file with input parameters calculated by gen_f
    # and passed to this sim_f.
    write_sim_input(input_file, H['x'])

    # Passed to command line in addition to the executable.
    # Here, only input file
    app_args = input_file
    os.environ["OMP_NUM_THREADS"] = machine_specs['OMP_NUM_THREADS']

    # Launch the executor to actually run the WarpX simulation
    if machine_specs['name'] == 'summit':
        task = exctr.submit(calc_type='sim',
                            extra_args=machine_specs['extra_args'],
                            app_args=app_args,
                            stdout='out.txt',
                            stderr='err.txt',
                            wait_on_run=True)
    else:
        task = exctr.submit(calc_type='sim',
                            num_procs=machine_specs['cores'],
                            app_args=app_args,
                            stdout='out.txt',
                            stderr='err.txt',
                            wait_on_run=True)

    # Periodically check the status of the simulation
    poll_interval = 1  # secs
    while(not task.finished):
        time.sleep(poll_interval)
        task.poll()
        if task.runtime > time_limit:
            task.kill()  # Timeout

    # Set calc_status with optional prints.
    if task.finished:
        if task.state == 'FINISHED':
            calc_status = WORKER_DONE
        elif task.state == 'FAILED':
            print("Warning: Task {} failed: Error code {}"
                  .format(task.name, task.errcode))
            calc_status = TASK_FAILED
        elif task.state == 'USER_KILLED':
            print("Warning: Task {} has been killed"
                  .format(task.name))
        else:
            print("Warning: Task {} in unknown state {}. Error code {}"
                  .format(task.name, task.state, task.errcode))

    # Safety
    time.sleep(0.2)

    try:
        # Get output from a run and delete output files
        warpx_out = read_sim_output(task.workdir)
    except Exception:
        warpx_out = np.nan
        print('Warning - output is Nan')

    # Pass the sim output values to LibEnsemble.
    # When optimization is ON, 'f' is then passed to the generating function
    # gen_f to generate new inputs for next runs.
    # All other parameters are here just for convenience.
    libE_output = np.zeros(1, dtype=sim_specs['out'])
    libE_output['f'] = warpx_out[0]
    libE_output['energy_std'] = warpx_out[1]
    libE_output['energy_avg'] = warpx_out[2]
    libE_output['charge'] = warpx_out[3]
    libE_output['emittance'] = warpx_out[4]
    libE_output['ramp_down_1'] = H['x'][0][0]
    libE_output['ramp_down_2'] = H['x'][0][1]
    libE_output['zlens_1'] = H['x'][0][2]
    libE_output['adjust_factor'] = H['x'][0][3]

    return libE_output, persis_info, calc_status
