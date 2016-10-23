#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Description

Details
"""
import datetime

__author__ = "Denys Berkovskyy"
__email__ = "berkovskyy@gmail.com"
__status__ = "Development"

import errno
import os
import shutil
import subprocess
import time
import string
import json
import csv

from functionparameters import FunctionParameters
from cluster_parser import ClusterQstatInfo


class GlobalSettings(object):
    """Storage for global Job Manager settings: projects list and matlab version"""
    RECENT_PROJECTS_COUNT = 10

    SETTINGS_FILENAME = 'matlab_job_manager.conf'

    def __init__(self):
        self._projects_list = list()
        self._matlab_version = None

    @classmethod
    def read_from_file(cls, filename: str=None):
        """Read settings from file"""
        if filename is None:
            filename = cls.SETTINGS_FILENAME
        with open(filename, 'r') as fin:
            settings_dict = json.load(fin)
        result = cls()
        result._projects_list = settings_dict['projects']
        result._matlab_version = settings_dict['matlab version']
        return result

    def write_to_file(self, filename: str=None) -> None:
        """Write settings to file"""
        if filename is None:
            filename = self.SETTINGS_FILENAME
        settings_dict = dict()
        settings_dict['projects'] = self._projects_list
        settings_dict['matlab version'] = self._matlab_version
        with open(filename, 'w') as fout:
            json.dump(settings_dict, fout)

    def get_projects(self) -> list:
        """Get list of recent projects as list

        The most recent project is last in list"""
        return self._projects_list

    def add_project(self, project_location: str) -> None:
        """Add project to settings

        Global settings store list of recent project. Maximum number of projects stored is
        defined in RECENT_PROJECTS_COUNT variable. If new project already exists in list it is
        moved to end of the list, otherwise it is appended to the list. First projects in list
        are removed until length of list is less than RECENT_PROJECTS_COUNT"""
        for i in range(self._projects_list.count(project_location)):
            self._projects_list.remove(project_location)
        self._projects_list.append(project_location)
        if len(self._projects_list) > self.RECENT_PROJECTS_COUNT:
            self._projects_list = self._projects_list[-self.RECENT_PROJECTS_COUNT:-1]

    def get_recent_project(self) -> str:
        """Get the most recent project path or None"""
        if len(self._projects_list):
            return self._projects_list[-1]
        else:
            return None

    def get_matlab_version(self) -> str:
        """Get version of matlab"""
        return self._matlab_version

    def set_matlab_version(self, version: str) -> None:
        """Set version of matlab"""
        self._matlab_version = version

    @staticmethod
    def check_headers_for_matlab_version(matlab_version):
        return os.path.exists(os.path.join('.', 'MatlabHeaders', '{}-Log.txt'.format(matlab_version))) and\
                os.path.exists(os.path.join('.', 'MatlabHeaders', '{}-Output.txt'.format(matlab_version)))


class ProjectSettings(object):
    CONF_FILENAME = 'project.conf'

    def __init__(self, settings: GlobalSettings, location: str):
        self._settings = settings
        self._location = location
        self._function_name = None
        self._function_parameters = None
        self._wall_seconds = 100 * 3600
        self._ppn = 1

        self._jobs = None

    def get_settings(self) -> GlobalSettings:
        return self._settings

    @classmethod
    def read_from_file(cls, settings: GlobalSettings, location: str, filename: str=None):
        if filename is None:
            filename = cls.CONF_FILENAME
        with open(os.path.join(location, filename), 'r') as fin:
            settings_dict = json.load(fin)
        result = cls(settings, location)
        result._function_name = settings_dict['function name']
        result._function_parameters = settings_dict['function parameters']
        result._wall_seconds = settings_dict['wall time seconds']
        result._ppn = settings_dict['ppn']
        result.reread_jobs()
        return result

    def reread_jobs(self):
        try:
            self._jobs = ProjectJobs.read_from_file(self)
        except IOError:
            self._jobs = None

    def write_to_file(self, filename: str=None) -> None:
        if filename is None:
            filename = self.CONF_FILENAME
        settings_dict = dict()
        settings_dict['function name'] = self._function_name
        settings_dict['function parameters'] = self._function_parameters
        settings_dict['wall time seconds'] = self._wall_seconds
        settings_dict['ppn'] = self._ppn
        with open(os.path.join(self._location, filename), 'w') as fout:
            json.dump(settings_dict, fout)

    def create_project_folder(self) -> None:
        folders_list = ['', 'Code', 'Data', 'Output', 'Log', 'Temp']
        for folder_name in folders_list:
            try:
                os.makedirs(os.path.join(self._location, folder_name))
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

    def check_output_folders_empty(self) -> bool:
        return os.listdir(os.path.join(self.get_location(), 'Log')) or \
               os.listdir(os.path.join(self.get_location(), 'Output')) or \
               os.listdir(os.path.join(self.get_location(), 'Temp'))

    def clear_output_folders(self) -> None:
        shutil.rmtree(os.path.join(self.get_location(), 'Log'))
        shutil.rmtree(os.path.join(self.get_location(), 'Output'))
        shutil.rmtree(os.path.join(self.get_location(), 'Temp'))
        self.create_project_folder()

    def set_function_name(self, value: str) -> None:
        self._function_name = value

    def set_function_parameters(self, value: str) -> None:
        self._function_parameters = value

    def get_location(self) -> str:
        return self._location

    def get_function_name(self) -> str:
        return self._function_name

    def get_function_parameters(self) -> str:
        return self._function_parameters

    def get_wall_time_str(self) -> str:
        if self._wall_seconds is None:
            return None
        else:
            s = self._wall_seconds % 60
            m = (self._wall_seconds // 60) % 60
            h = self._wall_seconds // 3600
            return '{}:{:0>2}:{:0>2}'.format(h, m, s)

    def get_wall_hours(self) -> int:
        if self._wall_seconds is None:
            return None
        else:
            return self._wall_seconds // 3600

    def set_wall_hours(self, value: int) -> None:
        self._wall_seconds = value * 3600

    def set_wall_minutes(self, value: int) -> None:
        self._wall_seconds = value * 60

    def set_wall_seconds(self, value: int) -> None:
        self._wall_seconds = value

    def set_wall_time(self, s: int, m: int=0, h: int=0) -> None:
        self._wall_seconds = s + m * 60 + h * 3600

    def set_ppn(self, value: int) -> None:
        self._ppn = value

    def get_ppn(self) -> int:
        return self._ppn

    def generate_jobs(self) -> None:
        self._jobs = ProjectJobs(self)
        parameters_generator = FunctionParameters.parse(self.get_function_parameters()).string_generator()
        internal_id = 0
        for param in parameters_generator:
            tjob = JobInfo(self, internal_id, '{}({})'.format(self._function_name, param))
            self._jobs.add_job(tjob)
            internal_id += 1

    def set_jobs(self, jobs) -> None:
        self._jobs = jobs

    def get_jobs(self):
        return self._jobs


class JobInfo(object):
    # Load job template
    with open('jobfile.pytmplt', 'r') as fin:
        JOB_TEMPLATE = string.Template(fin.read())

    def __init__(self, project: ProjectSettings, internal_id: int, function_string: str, cluster_id=None, state=None,
                 time_submit=None, time_start=None, time_run=None, hostname=None, queue=None,
                 description=None, notes=None):
        self.project = project
        self.internal_id = internal_id
        self.function_string = function_string
        self.cluster_id = cluster_id
        self.state = state
        self.time_submit = time_submit
        self.time_start = time_start
        self.time_run = time_run
        self.hostname = hostname
        self.queue = queue
        self.description = description
        self.notes = notes

    def start_job(self):
        # Check function exists
        matlab_func_location = os.path.join(self.project.get_location(), 'Code',
                                            self.function_string.split('(')[0] + '.m')
        if not os.path.exists(matlab_func_location):
            raise Exception('Matlab function file \'{}\' is missing'.format(matlab_func_location))

        # Create job string by substituting values into template
        job_string = self.JOB_TEMPLATE.substitute(internal_id=self.internal_id,
                                                  project_path=self.project.get_location(),
                                                  function_string=self.function_string,
                                                  matlab_version=self.project.get_settings().get_matlab_version(),
                                                  wall_time=self.project.get_wall_time_str(),
                                                  ppn=self.project.get_ppn())

        # Start job by submitting it with qsub command
        job_obj = subprocess.Popen(
            ['qsub', '-N', str(self.internal_id), '-d', os.path.join(self.project.get_location(), 'Temp'),
             #'-o', os.path.join(self.project.get_location(),
             #                   'Temp', '{}_job_output.txt'.format(self.internal_id)),
             #'-e', os.path.join(self.project.get_location(),
             #                   'Temp', '{}_job_error.txt'.format(self.internal_id)),
             '-'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)

        # Send script to qsub and get output
        job_out, job_err = job_obj.communicate(job_string.encode('UTF8'))
        if job_obj.returncode != 0 or len(job_err) != 0:
            raise Exception(
                'Error submitting job. Return code {};\nstdout: {}\nerrorout: {}'.format(job_obj.returncode, job_out,
                                                                                         job_err))
        # Extract job id and update job object
        job_out = job_out.decode('utf8').strip()
        self.time_submit = datetime.datetime.utcnow()
        self.cluster_id = job_out
        return job_out


class ProjectJobs(object):
    def __init__(self, project: ProjectSettings):
        self._project = project
        self._jobs = list()
        self._cluster_info_present = False

    def add_job(self, job: JobInfo) -> None:
        self._jobs.append(job)

    @classmethod
    def read_from_file(cls, project: ProjectSettings):
        result = cls(project)
        with open(os.path.join(project.get_location(), 'Jobs.csv'), 'r', newline='') as fin:
            jobreader = csv.reader(fin)
            header_row = True
            for row in jobreader:
                if header_row:
                    header_row = False
                else:
                    internal_id, function_string, clusted_id, time_submit, time_start, time_run, hostname, queue,\
                        job_description, job_notes = [None if len(x) == 0 else x for x in row]
                    time_submit, time_start = [
                        datetime.datetime.strptime(x, '%Y-%m-%d %H:%M:%S.%f') if x is not None else None for x in
                        [time_submit, time_start]]
                    if time_run is not None:
                        time_run = float(time_run)
                    tjob = JobInfo(project, int(internal_id), function_string, clusted_id, None,
                                   time_submit, time_start, time_run, hostname, queue, job_description, job_notes)
                    result.add_job(tjob)
        return result

    def write_to_file(self) -> None:
        with open(os.path.join(self._project.get_location(), 'Jobs.csv'), 'w', newline='') as fout:
            jobwriter = csv.writer(fout)
            jobwriter.writerow(
                ['Internal ID', 'Function', 'Cluster ID',
                 'Submit time', 'Start time', 'Processing Time', 'Hostname', 'Queue', 'Description', 'Notes'])
            for job in self._jobs:
                time_list = [None if x is None else x.strftime('%Y-%m-%d %H:%M:%S.%f') for x in
                             [job.time_submit, job.time_start]]
                jobwriter.writerow(
                    [job.internal_id, job.function_string, job.cluster_id,
                     time_list[0], time_list[1], job.time_run, job.hostname, job.queue, job.description, job.notes])

    def start_jobs(self):
        first_job = True
        for job in self._jobs:
            if job.cluster_id is None:
                if first_job:
                    first_job = False  # No need to wait before first job
                else:
                    time.sleep(1)  # Wait before submitting job
                yield job.internal_id, job.start_job()

    def get_total_jobs_count(self) -> int:
        return len(self._jobs)

    def get_submitted_jobs_count(self) -> int:
        result = 0
        for job in self._jobs:
            if job.cluster_id is not None:
                result += 1
        return result

    def get_running_jobs_count(self) -> int:
        result = 0
        for job in self._jobs:
            if job.state == 'R':
                result += 1
        return result

    def get_finished_jobs_count(self) -> int:
        result = 0
        for job in self._jobs:
            if job.time_run is not None:
                result += 1
        return result

    def get_finished_jobs_time(self) -> float:
        result = 0
        for job in self._jobs:
            if job.time_run is not None:
                result += job.time_run
        return result

    def update_jobs_info(self):
        cluster_jobs_info = ClusterQstatInfo()
        self._cluster_info_present = True
        for job in self._jobs:
            if job.time_run is None:
                num_id = job.cluster_id.split('.')[0]
                tjob = cluster_jobs_info.get_job_info(num_id)
                if tjob is None:
                    job.state = 'C'
                else:
                    job.state = tjob.job_state
                if job.state == 'C':
                    logs_copied = False  # Flag which indicates that some logs need attention

                    pdir = self._project.get_location()  # Project dir
                    matlab_version = self._project.get_settings().get_matlab_version()

                    # Define filenames for files
                    filenames = {
                        'info_in': (pdir, 'Temp', '{}_{}_info.txt'.format(job.internal_id, num_id)),
                        'description_in': (pdir, 'Temp', '{}_{}_description.txt'.format(job.internal_id, num_id)),

                        'job_output_in': (pdir, 'Temp', '{}.o{}'.format(job.internal_id, num_id)),
                        'job_output_out': (pdir, 'Log', '{}_job_output.txt'.format(job.internal_id)),
                        'job_error_in': (pdir, 'Temp', '{}.e{}'.format(job.internal_id, num_id)),
                        'job_error_out': (pdir, 'Log', '{}_job_error.txt'.format(job.internal_id)),

                        'matlab_output_in': (pdir, 'Temp', '{}_{}_matlab_output.txt'.format(job.internal_id, num_id)),
                        'matlab_output_out': (pdir, 'Log', '{}_matlab_output.txt'.format(job.internal_id)),
                        'matlab_error_in': (pdir, 'Temp', '{}_{}_matlab_error.txt'.format(job.internal_id, num_id)),
                        'matlab_error_out': (pdir, 'Log', '{}_matlab_error.txt'.format(job.internal_id)),
                        'matlab_log_in': (pdir, 'Temp', '{}_{}_matlab_log.txt'.format(job.internal_id, num_id)),
                        'matlab_log_out': (pdir, 'Log', '{}_matlab_log.txt'.format(job.internal_id)),

                        'matlab_ref_log': ('.', 'MatlabHeaders', '{}-Log.txt'.format(matlab_version)),
                        'matlab_ref_output': ('.', 'MatlabHeaders', '{}-Output.txt'.format(matlab_version))
                    }

                    # Join path parts
                    for key in filenames:
                        filenames[key] = os.path.join(*filenames[key])

                    move_files_list = list()  # List of files which will be moved: (src, dst) tuples
                    delete_files_list = list()  # List of files which will be deleted

                    # List of files which should be empty as tuples (src, dst)
                    # File will be copied to dst if file is not empty
                    empty_files_list = [('job_output_in', 'job_output_out'),
                                        ('job_error_in', 'job_error_out'),
                                        ('matlab_error_in', 'matlab_error_out')]
                    for fileptr in empty_files_list:
                        is_empty = True
                        with open(filenames[fileptr[0]]) as fin:
                            with open(filenames[fileptr[1]], 'w') as fout:
                                for line in fin:
                                    if not (fileptr[0] == 'job_error_in' and is_empty and
                                            line.strip().startswith('Warning: Permanently added \'') and
                                            line.strip().endswith('\' (RSA) to the list of known hosts.')):
                                        is_empty = False
                                        fout.write(line)
                        delete_files_list.append(fileptr[0])
                        if is_empty:
                            delete_files_list.append(fileptr[1])
                        else:
                            logs_copied = True

                    # List of files which will be checked against reference as tuples (reference, src, dst)
                    # Difference will be copied to dst if file is not same as reference
                    matlab_check_filelist = [('matlab_ref_log', 'matlab_log_in', 'matlab_log_out'),
                                             ('matlab_ref_output', 'matlab_output_in', 'matlab_output_out')]
                    for fileptr in matlab_check_filelist:
                        if os.path.exists(filenames[fileptr[0]]):  # Do check only if reference file exists, otherwise move file
                            difference_found = False
                            with open(filenames[fileptr[0]]) as fref:
                                with open(filenames[fileptr[1]]) as fin:
                                    with open(filenames[fileptr[2]], 'w') as fout:
                                        for line in fin:
                                            if not difference_found:
                                                if line.strip() != fref.readline().strip():
                                                    difference_found = True
                                            if difference_found:
                                                fout.write(line)
                            if difference_found:
                                logs_copied = True
                            else:
                                delete_files_list.append(fileptr[2])
                            delete_files_list.append(fileptr[1])
                        else:
                            move_files_list.append((fileptr[1], fileptr[2]))

                    # Get and update job info
                    with open(filenames['info_in']) as fin:
                        cluster_id, time_start, time_work, hostname, queue = [x.strip() for x in list(fin)]
                    if job.cluster_id != cluster_id:
                        raise Exception('Cluster IDs do not match: \'{}\' != \'{}\''.format(job.cluster_id, cluster_id))
                    job.time_start = datetime.datetime.strptime(time_start, '%Y-%m-%d %H:%M:%S.%f')
                    job.time_run = float(time_work)
                    job.hostname = hostname
                    job.queue = queue
                    if logs_copied:
                        job.notes = 'Check Logs'
                    delete_files_list.append('info_in')

                    try:
                        with open(filenames['description_in']) as fin:
                            job.description = fin.read().strip()
                        delete_files_list.append('description_in')
                    except IOError:
                        pass

                    for fileptr in move_files_list:
                        shutil.move(filenames[fileptr[0]], filenames[fileptr[1]])

                    for fileptr in delete_files_list:
                        os.remove(filenames[fileptr])
            self.write_to_file()


def main():
    raise NotImplementedError()


if __name__ == '__main__':
    main()
