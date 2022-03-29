import subprocess

from pymol import cmd


def run_ring_local(ring_pth, file_pth, obj_name, run_config, log_f, progress_f):
    progress_f(1)

    p = subprocess.Popen(
            [ring_pth, "-i", file_pth, "--out_dir", "/tmp/ring/",
             "-g", run_config["-g"],
             "-o", run_config["-o"],
             "-s", run_config["-s"],
             "-k", run_config["-k"],
             "-a", run_config["-a"],
             "-b", run_config["-b"],
             "-w", run_config["-w"],
             "--all_chains", run_config["edges"], "--all_models", "--md"], stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE, universal_newlines=True)

    log_f("Local RING generation started")

    n_states = cmd.count_states(obj_name)
    prev_state = 0
    while p.poll() is None:
        line = p.stderr.readline()
        if line != "":
            if "model" in line:
                current_state = int(line.split("model ")[1].strip())
                if current_state > prev_state:
                    progress_f((current_state / n_states) * 100)
                    prev_state = current_state
