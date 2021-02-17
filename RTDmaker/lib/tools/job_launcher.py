import os
import sys
import time
import multiprocessing
from subprocess import Popen
from multiprocessing.pool import ThreadPool


def work(command, logfile, job_id, tot):

    line_info = f'\n{time.asctime()} Starting Job {job_id} (out of {tot})'
    line_command = f'\n{time.asctime()} Job {job_id} command: {command}\n'

    with open(logfile, "a+") as fh:
        fh.write(line_info)
        fh.write(line_command)

    try:
        process = Popen(command, shell=True)
        process.wait()
        # If the output of the process needs further processing/parsing, it can be done here
        # Source: https://stackoverflow.com/questions/26774781/
        # python-multiple-subprocess-with-a-pool-queue-recover-output-as-soon-as-one-finis
    except Exception as e:
        line_error = f'\n{time.asctime()} Error while executing Job {job_id}:\n'
        print(line_error)
        with open(logfile, "a+") as fh:
            fh.write(line_error)
            fh.write(str(e) + "\n")

    line_end = f'\n{time.asctime()} Job {job_id} completed!\n'
    print(line_end)

    with open(logfile, "a+") as fh:
        fh.write(line_end)


def launch_jobs(commands_list, logfile=None, n_jobs=None, core_proportion=(1, 3), max_cores=8, log_dir=None):

    if not log_dir:
        log_dir = os.getcwd()

    # Create a log file to track the completed jobs
    if not logfile:
        logfile = os.path.join(log_dir, f"{time.asctime().replace(' ', '_')}_jobs_logfile.txt")

    # Attach ID to the jobs to track their execution
    indexed_commands = [(i + 1, command) for i, command in enumerate(commands_list)]
    tot = len(indexed_commands)

    if n_jobs:
        if n_jobs > multiprocessing.cpu_count():
            sys.exit(f"The system does not posses that many cores. It must be {multiprocessing.cpu_count()} or less.")
        else:
            n_cores = n_jobs
    else:
        # Use a predetermined proportion of available total cores (ex: one third (1/3) of the available cores)
        n_cores = int((multiprocessing.cpu_count() / core_proportion[1]) * core_proportion[0])

    # Limit the maximum number of cores to use at a time
    if n_cores > max_cores:
        n_cores = max_cores

    line_start = f'{time.asctime()} Launching {tot} jobs, using {n_cores} cores\n'
    print(line_start)

    with open(logfile, "a+") as fh:
        fh.write(line_start)

    # Launch "n" number of jobs at a time; whenever a job is finish, launch a new one
    # The number of jobs is determined by the number of available/selected cores to use
    tp = ThreadPool(n_cores)
    for (job_id, command) in indexed_commands:
        tp.apply_async(work, (command, logfile, job_id, tot, ))
    tp.close()
    tp.join()

    line_end = f'\n{time.asctime()} All Jobs completed! ({tot})\n'
    print(line_end)

    with open(logfile, "a+") as fh:
        fh.write(line_end)
