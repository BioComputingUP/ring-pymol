import os
import subprocess
import sys
import time

from pymol import cmd
from tqdm import tqdm


def run_ring_local(ring_pth, file_pth, obj_name, run_config, tmp_dir, log_f, progress_f, verbose=False):
    progress_f(1)

    p = subprocess.Popen(
        [ring_pth, "-i", file_pth, "--out_dir", tmp_dir,
         "-g", run_config["-g"],
         "-o", run_config["-o"],
         "-s", run_config["-s"],
         "-k", run_config["-k"],
         "-a", run_config["-a"],
         "-b", run_config["-b"],
         "-w", run_config["-w"],
         run_config["water"],
         run_config["add_h"],
         run_config["edges"], "--all_models", "--md", "-v", "--no_colors"],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, universal_newlines=True)

    log_f("Local RING generation started")

    if verbose:
        log_f("Executing : " + ' '.join(p.args))

    log_file = open(tmp_dir + '/ring.log', 'w')

    for line in iter(p.stdout.readline, ''):
        if line.startswith('{') and '}' in line:
            progress = float(line.split('%')[0].split('{')[1].strip())
            progress_f(progress)
            if verbose and '(' in line and ')' in line:
                times = line.split('(')[1].split(')')[0].split('<')
                log_f('Elapsed : ' + times[0].strip() + ' | Remaining : ' + times[1].strip())
        if 'warning:' in line:
            warn = line.split('warning:')[1].strip()
            log_file.write(warn + '\n')
            if verbose:
                log_f(warn, warning=True)

    p.wait()
    log_file.close()

    # Check if the log file is empty
    if os.stat(tmp_dir + '/ring.log').st_size != 0:
        log_f(f"Local RING generation finished with warnings, check log file {tmp_dir}/ring.log for details", warning=True)

    if p.returncode != 0:
        raise Exception("Local RING generation failed")
