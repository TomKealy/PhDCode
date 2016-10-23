#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Parser of BlueCrystal queue management commands output"""

__author__ = "Denys Berkovskyy"
__email__ = "berkovskyy@gmail.com"
__status__ = "Development"


import subprocess


class OutputParseException(Exception):
    """Parser exception"""
    def __init__(self, cmd: str, line: int, msg: str=None):
        super().__init__('Error parsing output of command \'{}\': {}. (Line {})'.format(cmd, msg, line))


class MatlabVersions(object):
    '''Access list of available matlab versions'''

    def __init__(self):
        self._versions = list()

        cmd_obj = subprocess.Popen(['module avail'], shell=True,
                                         stdout=None, stdin=None, stderr=subprocess.PIPE)
        outs, errs = cmd_obj.communicate()
        if cmd_obj.returncode != 0:
            raise Exception()

        lines = errs.decode().splitlines()
        for line in lines:
            line = line.strip()
            if line and not line.startswith('-'):
                if line.startswith('apps/matlab-'):
                    self._versions.append(line[len('apps/matlab-'):])

    def get_versions(self):
        return self._versions

    def check_version_available(self, version):
        return version in self._versions


class ClusterQstatInfo(object):
    """Access to qstat command results"""
    class ClusterJobInfo(object):
        """Individual job info"""
        def __init__(self, job_id, job_name, job_user, job_time, job_state, job_queue):
            self.job_id = job_id
            self.job_name = job_name
            self.job_user = job_user
            self.job_time = job_time
            self.job_state = job_state
            self.job_queue = job_queue

    def __init__(self):
        output = subprocess.check_output(['qstat']).decode('utf8')
        lines = output.splitlines()
        current_line = 0
        if lines[current_line].strip() != 'Job ID                    Name             User            Time Use S Queue':
            raise OutputParseException('qstat', current_line)
        current_line += 1
        if lines[current_line].strip() != '------------------------- ---------------- --------------- -------- - -----':
            raise OutputParseException('qstat', current_line)
        current_line += 1

        jobs_count = len(lines) - 2
        self._jobs_ids = dict()
        self._jobs_users = dict()
        for i in range(jobs_count):
            fields = [x.strip() for x in lines[current_line].split()]
            current_line += 1
            tjob = self.ClusterJobInfo(*fields)
            self._jobs_ids[tjob.job_id.split('.')[0]] = tjob
            if tjob.job_user not in self._jobs_users:
                self._jobs_users[tjob.job_user] = set()
            self._jobs_users[tjob.job_user].add(tjob)
        if current_line != len(lines):
            raise OutputParseException('qstat', current_line, 'Unknown lines after end')

    def get_job_info(self, job_id):
        return self._jobs_ids.get(job_id)

    def users_job_count_stat(self):
        result = dict()
        for user in self._jobs_users.keys():
            result[user] = len(self._jobs_users[user])
        return result


class ClusterShowqInfo(object):
    '''Access shoq command results'''

    @staticmethod
    def parse_table(lines, current_line, titles, headers):
        # Section title
        good_title = False
        for title in titles:
            if lines[current_line].strip().lower().startswith(title):
                good_title = True
                break
        if not good_title:
            raise OutputParseException('showq', current_line, 'Line {} does not start with {}'.format(lines[current_line], title))
        current_line += 1

        # Table header
        good_header = False
        for header in headers:
            if lines[current_line].lower().split() == header:
                good_header = True
                break
        if not good_header:
            raise OutputParseException('showq', current_line)
        current_line += 1

        # Empty line
        if lines[current_line].strip():
            raise OutputParseException('showq', current_line)
        current_line += 1

        body = list()
        # Table content
        while lines[current_line].strip():
            body.append(lines[current_line].strip())
            current_line += 1
        current_line += 1

        footer = list()
        while lines[current_line].strip():
            footer.append(lines[current_line].strip())
            current_line += 1

        current_line += 1

        return body, footer, current_line

    def __init__(self):
        output = subprocess.check_output('showq').decode('utf8')
        lines = output.splitlines()
        current_line = 0

        # Find first non empty line
        while not(lines[current_line].strip()):
            current_line += 1

        body, footer, current_line =\
            ClusterShowqInfo.parse_table(lines, current_line, ['active jobs'],
                                        [['jobid', 'username', 'state', 'procs', 'remaining', 'starttime'],
                                         ['jobname', 'username', 'state', 'proc', 'remaining', 'starttime']])
        for line in body:
            # 2266367              sw1496    Running    32    12:29:21  Thu Dec 19 06:07:49
            values = line.split()
            if len(values) != 9:
                raise OutputParseException('showq', current_line, 'Cannot parse job line, number of items in line is {}'.format(len(values)))
            job_id = values[0]
            job_owner = values[1]
            job_state = values[2]
            job_proc = int(values[3])
            job_remaining = values[4]
            job_start_time = values[5] + ' ' + values[6] + ' ' + values[7] + ' ' + values[8]

        body, footer, current_line =\
            ClusterShowqInfo.parse_table(lines, current_line, ['eligible', 'idle jobs'],
                                        [['jobid', 'username', 'state', 'procs', 'wclimit', 'queuetime'],
                                         ['jobname', 'username', 'state', 'proc', 'wclimit', 'queuetime']])
        for line in body:
            # 2242400             sl12655  BatchHold    32 15:00:00:00  Tue Nov 19 16:57:49
            values = line.split()
            if len(values) != 9:
                raise OutputParseException('showq', current_line, 'Cannot parse job line, number of items in line is {}'.format(len(values)))
            job_id = values[0]
            job_owner = values[1]
            job_state = values[2]
            job_proc = int(values[3])
            job_wclimit = values[4]
            job_queuetime = values[5] + ' ' + values[6] + ' ' + values[7] + ' ' + values[8]

        body, footer, current_line =\
            ClusterShowqInfo.parse_table(lines, current_line, ['blocked jobs'],
                                        [['jobid', 'username', 'state', 'procs', 'wclimit', 'queuetime'],
                                         ['jobname', 'username', 'state', 'proc', 'wclimit', 'queuetime']])
        for line in body:
            # 2242400             sl12655  BatchHold    32 15:00:00:00  Tue Nov 19 16:57:49
            values = line.split()
            if len(values) != 9:
                raise OutputParseException('showq', current_line, 'Cannot parse job line, number of items in line is {}'.format(len(values)))
            job_id = values[0]
            job_owner = values[1]
            job_state = values[2]
            job_proc = int(values[3])
            job_wclimit = values[4]
            job_queuetime = values[5] + ' ' + values[6] + ' ' + values[7] + ' ' + values[8]

        current_line += 1
        #if current_line

        #while current_line < len(lines) and not lines[current_line]:
        #    current_line += 1

        #if current_line != len(lines):
        #   raise OutputParseException('showq', current_line, 'Unknown lines after end')


def main():
    raise NotImplementedError()

if __name__ == "__main__":
    main()
