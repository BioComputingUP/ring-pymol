import os
import shutil
import threading
import time
from threading import Thread

import requests

drmaatic_url = "https://drmaatic.biocomputingup.it"
list_job_url = "{}/job/".format(drmaatic_url)

_log_f = lambda x: x


class Status:
    statusMap = {
        "job has been rejected from the ws": "failed",
        "job has been received from the ws": "pending",
        "job has been created and sent to the DRM": "pending",
        "process status cannot be determined": "pending",
        "job is queued and active": "running",
        "job is queued and in system hold": "running",
        "job is queued and in user hold": "running",
        "job is queued and in user and system hold": "running",
        "job is running": "running",
        "job is system suspended": "pending",
        "job is user suspended": "pending",
        "job finished normally": "success",
        "job finished, but failed": "failed",
        "job has been deleted": "deleted"
    }

    def __init__(self, status):
        self.status = self.decode_status(status)

    def __repr__(self):
        return self.status

    def __eq__(self, other):
        return self.status == other

    def decode_status(self, status_long):
        return self.statusMap[status_long]


class Job:
    _status: [Status, None] = None
    _uuid: [str, None] = None

    def __init__(self, uuid=None, status=None):
        self.uuid = uuid
        self.status = status

    def __repr__(self) -> str:
        return "{} - {}".format(self.uuid, self.status)

    @property
    def status(self):
        return self._status

    @status.setter
    def status(self, status):
        self._status = Status(status) if status is not None else status

    @property
    def uuid(self):
        return self._uuid

    @uuid.setter
    def uuid(self, uuid):
        self._uuid = uuid

    def is_finished(self) -> bool:
        return self._status == "failed" or self._status == "deleted" or self._status == "success"


def check_for_job(job):
    try:
        job_url = "{}/{}".format(list_job_url, job.uuid)

        while not job.is_finished():
            response = requests.get(job_url, timeout=500)
            response.raise_for_status()
            job.status = response.json()["status"]
            if not job.is_finished():
                time.sleep(3)

    except requests.exceptions.RequestException as err:
        return err


def post_job(job, file_pth, params):
    try:
        files = {'input_file': open(file_pth, 'rb')}

        response = requests.post(list_job_url, files=files, data=params, timeout=10000)
        response.raise_for_status()
        job.uuid = response.json()["uuid"]
        job.status = response.json()["status"]

    except requests.exceptions.RequestException as err:
        return err


def download_results(job, extract_pth):
    if job.status == "failed":
        return
    try:
        output_url = "{}/{}/{}".format(list_job_url, job.uuid, "download")

        response = requests.get(output_url, timeout=5)
        response.raise_for_status()
        file_name = response.headers["content-disposition"].split("filename=")[1]
        file_pth = "{}/{}".format(extract_pth, file_name)

        with open(file_pth, "wb") as f:
            f.write(response.content)

        shutil.unpack_archive(file_pth, extract_pth)
        os.remove(file_pth)

    except requests.exceptions.RequestException as err:
        return err


def config_to_parameters(config):
    convert = {"-g": "seq_sep",
               "-o": "len_salt",
               "-s": "len_ss",
               "-k": "len_pipi",
               "-a": "len_pica",
               "-b": "len_hbond",
               "-w": "len_vdw"}

    new_config = {}

    for key, value in config.items():
        # replace key of config with value
        if key in convert:
            new_config[convert[key]] = value.strip("--")

    new_config[config["edges"].strip("--")] = True

    return new_config


def run_ring_api(file_pth, run_config, tmp_dir, log_f, progress_f):
    global _log_f

    _log_f = log_f

    job: Job = Job()

    file_name = os.path.basename(file_pth)
    _log_f(file_pth, file_name)

    parameters = {"task": "ring-plugin-api",
                  "original_name": file_name
                  }

    parameters.update(config_to_parameters(run_config))

    _log_f("Remote RING generation started")
    _log_f("Sending job to remote server")
    t_post_job = Thread(target=post_job, args=(job, file_pth, parameters))
    t_post_job.start()

    prev_progress = 0
    while t_post_job.is_alive():
        progress_f(min([prev_progress, 15]))
        prev_progress = (prev_progress + 0.01)

    t_check_job = Thread(target=check_for_job, args=(job,))
    t_check_job.start()

    prev_progress = 15
    timer = time.time() - 5
    while t_check_job.is_alive():
        if time.time() - timer > 5:
            timer = time.time()
            _log_f("Running job {}".format(job))

        progress_f(min([prev_progress, 85]))
        prev_progress = (prev_progress + 0.00001)

    if job.status == "success":
        _log_f("Computation terminated, downloading results")
    else:
        _log_f("Error in the execution of RING, please retry later or launch locally", error=True)

    t_download_results = Thread(target=download_results, args=(job, tmp_dir))
    t_download_results.start()

    prev_progress = 85

    _log_f("Downloading results...")
    while t_download_results.is_alive():
        progress_f(min([prev_progress, 100]))
        prev_progress = (prev_progress + 0.01)

    _log_f("Download completed!")
    progress_f(100)
