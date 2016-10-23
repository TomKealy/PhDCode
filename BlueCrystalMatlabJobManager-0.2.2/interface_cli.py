#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Command Line Interface for BlueCrystal Matlab Job Manager"""

__author__ = "Denys Berkovskyy"
__email__ = "berkovskyy@gmail.com"
__status__ = "Development"

import os
import datetime
import time

from functionparameters import FunctionParameters

from job_manager import GlobalSettings, ProjectSettings
import cluster_parser


INFO_MESSAGE = "Use this script to create and start and manage BlueCrystal matlab jobs. " \
               "This script creates 'Code', 'Data', 'Output', 'Log', 'Temp' folders. " \
               "Copy your matlab code to 'Code' folder and your data to 'Data' folder. " \
               "Your code should write results to 'Output' folder. This script will store " \
               "temporary files in 'Temp' folder and will copy any useful logs to 'Logs' folder. \n\n" \
               "This script will ask number of questions. For multiple selection question " \
               "you can make selection either using number or text. For some questions " \
               "script will show last used value in square brackets. You can use this value " \
               "by entering empty string (pressing enter).\n\n" \
               "You can specify function parameters either as ranges, " \
               "for which script will generate intermediate values or as list of string, which " \
               "script will use 'as is'.\n\n" \
               "!NB! matlab working dir will be project folder and not 'Code' folder."


def read_input_with_hint(prompt, hint):
    """Ask user a question optionally show a suggested value.

    Suggested value will be shown in square brackets and will
    be used if user enters empty string."""

    if hint is not None:
        hint_str = ' [\'{}\']'.format(hint)
    else:
        hint_str = ''
    result = None
    while not result:
        result = input(prompt + hint_str + ':\n')
        if not result:
            result = hint
    return result.strip()


class SelectionOption(object):
    """Class for representation of user selectable options.

    key is primary key (usually single character) which is returned as a result
    msg is description of the option
    aliases is text entries which will be accepted as user inputs"""

    def __init__(self, key: str, msg: str, aliases: list):
        self.key = key
        self.msg = msg
        self.aliases = aliases


def get_selection(prompt: str, options: list):
    """Ask user a question and get answer from defined options.

    The options will be shown as list under question"""

    good_inputs = list()

    full_prompt = list()
    full_prompt.append(prompt)

    for opt in options:
        good_inputs.append(opt.key)
        for als in opt.aliases:
            good_inputs.append(als)
        full_prompt.append(' {}: {}'.format(opt.key, opt.msg))
    full_prompt.append('')

    full_prompt = '\n'.join(full_prompt)

    result = None
    while True:
        user_input = input(full_prompt).lower()
        if user_input in good_inputs:
            for opt in options:
                if user_input == opt.key or user_input in opt.aliases:
                    result = opt.key
                    break
            break
        print('Unknown option: {}'.format(user_input))
    return result


def get_yes_no_input(prompt) -> bool:
    """Ask user a question and get a yes/no answer. Returns bool."""

    options = [SelectionOption('0', 'No', ['n', 'no']),
               SelectionOption('1', 'Yes', ['y', 'yes'])]
    user_selection = get_selection(prompt, options)
    if user_selection == '0':
        return False
    elif user_selection == '1':
        return True
    else:
        raise NotImplementedError(user_selection)


def restart_project(settings: GlobalSettings, project: ProjectSettings):
    """Start/Restart project

    Start new project. If previous project exists in the folder
    values from the project will be suggested as answers for questions"""

    if project.check_output_folders_empty():
        do_clear = get_yes_no_input('Do you want to clear working directories (Log, Output, Temp)?')
        if do_clear:
            project.clear_output_folders()

    project.set_function_name(read_input_with_hint('Please enter matlab function name', project.get_function_name()))
    print('Using \'{}\' as matlab function name'.format(project.get_function_name()))

    function_parameters = None
    while not FunctionParameters.verify(function_parameters):
        function_parameters = read_input_with_hint('Please enter matlab function parameters',
                                                   project.get_function_parameters())
    print('Using \'{}\' as matlab function parameters'.format(function_parameters))
    project.set_function_parameters(function_parameters)

    matlab_versions = cluster_parser.MatlabVersions().get_versions()
    if settings.get_matlab_version() is None and len(matlab_versions):
        settings.set_matlab_version(matlab_versions[-1])
    settings.set_matlab_version(read_input_with_hint('Please enter matlab version ({})'.format(' / '.join(matlab_versions)),
                                                     settings.get_matlab_version()))
    print('Using \'{}\' as matlab version'.format(settings.get_matlab_version()))
    if not GlobalSettings.check_headers_for_matlab_version(settings.get_matlab_version()):
        print('Warning: matlab log headers are for matlab version {} are missing, \
               log files will not be deleted automatically'.format(settings.get_matlab_version()))

    project.set_wall_hours(int(read_input_with_hint('Please enter maximum time per job in hours', str(project.get_wall_hours()))))
    print('Using \'{}\' as wall time'.format(project.get_wall_time_str()))

    project.set_ppn(int(read_input_with_hint('Please enter processors per node value', str(project.get_ppn()))))
    print('Using \'{}\' as ppn'.format(project.get_ppn()))

    settings.write_to_file()
    project.write_to_file()

    project.generate_jobs()
    project.get_jobs().write_to_file()

    do_submit = get_yes_no_input('Do you want to submit jobs?')
    if do_submit:
        continue_project(project)


def continue_project(project: ProjectSettings):
    """Continue project: submit jobs and show results"""
    total_jobs = project.get_jobs().get_total_jobs_count()
    try:
        first_job = True
        for internal_id, job_id in project.get_jobs().start_jobs():
            if first_job:
                print('Submitting jobs...')
                first_job = False
            time_left = datetime.timedelta(0, total_jobs - internal_id - 1)
            print('Submitted job {}/{} (cluster id is {}). Time left: {}'.format(internal_id + 1, total_jobs,
                                                                                 job_id, time_left))
            project.get_jobs().write_to_file()

        print('All jobs submitted, showing updated stats every 5 sec (You can stop at any time by pressing Ctrl-C)')
        while True:
            project.get_jobs().update_jobs_info()
            print_jobs_stat(project.get_jobs())
            if project.get_jobs().get_finished_jobs_count() == project.get_jobs().get_total_jobs_count():
                print('Done')
                break
            time.sleep(5)
    except KeyboardInterrupt:
        pass
    print_jobs_stat(project.get_jobs())


def print_jobs_stat(jobs):
    print('Jobs stats: {}/{}({})/{} (finished/submitted(running)/total). Total time of finished jobs is {}'
          .format(jobs.get_finished_jobs_count(), jobs.get_submitted_jobs_count(), jobs.get_running_jobs_count(),
                  jobs.get_total_jobs_count(), datetime.timedelta(0, jobs.get_finished_jobs_time())))


def interface():
    print(INFO_MESSAGE)

    # Load last settings
    try:
        settings = GlobalSettings.read_from_file()
    except (IOError, ValueError):
        settings = GlobalSettings()

    # Determine project location
    project_location = read_input_with_hint('Please enter project location relative to home folder',
                                            settings.get_recent_project())
    project_dir = os.path.join(os.path.expanduser('~'), project_location)
    print('Using \'{}\' as project location'.format(project_dir))
    settings.add_project(project_dir)
    settings.write_to_file()

    try:
        project = ProjectSettings.read_from_file(settings, project_dir)
    except (IOError, ValueError):
        project = ProjectSettings(settings, project_dir)
    project.create_project_folder()

    if project.get_jobs() is not None:
        options = [SelectionOption('1', 'Review parameters and restart project', ['restart']),
                   SelectionOption('2', 'Continue with current project', ['continue'])]
        user_selection = get_selection('Found active project in the folder. What do you want to do?', options)
        if user_selection == '1':
            restart_project(settings, project)
        elif user_selection == '2':
            continue_project(project)
        else:
            raise NotImplementedError(user_selection)
    else:
        restart_project(settings, project)


def main():
    interface()


if __name__ == '__main__':
    main()