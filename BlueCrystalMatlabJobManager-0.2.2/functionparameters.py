#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Function parameters parser and generator"""

__author__ = "Denys Berkovskyy"
__email__ = "berkovskyy@gmail.com"
__status__ = "Development"


class ParametersParsingException(Exception):
    pass


class FunctionParameterRange():
    """Storage of parameter with range. Intermediate steps will be generated for such parameter"""
    def __init__(self):
        self.add_start_stop = False
        self.start = None
        self.stop = None
        self.step = None


class FunctionParameterList():
    """Storage of list parameters. Raw values will be used"""
    def __init__(self):
        self.values = list()


class FunctionParameters():
    """Storage of list of parameters. Each of parameter could be ranged or listed"""
    def __init__(self):
        self.entries = list()

    def string_generator(self, current_string=None):
        if current_string is None:
            yield from self.string_generator(list())
        else:
            if len(current_string) < len(self.entries):
                if isinstance(self.entries[len(current_string)], FunctionParameterRange):
                    entry = self.entries[len(current_string)]
                    if entry.stop is None:
                        yield from self.string_generator(current_string + [str(entry.start)])
                    else:
                        for value in range(entry.start, entry.stop, entry.step):
                            if entry.add_start_stop:
                                tend = value + entry.step - 1
                                if tend > entry.stop:
                                    tend = entry.stop
                                yield from self.string_generator(current_string + [str(value) + ',' + str(tend)])
                            else:
                                yield from self.string_generator(current_string + [str(value)])
                        if not entry.add_start_stop:
                            yield from self.string_generator(current_string + [str(entry.stop)])
                elif isinstance(self.entries[len(current_string)], FunctionParameterList):
                    for entry in self.entries[len(current_string)].values:
                        yield from self.string_generator(current_string + [str(entry)])
                else:
                    raise Exception('Unknown parameter type')
            else:
                yield ','.join(current_string)

    @classmethod
    def parse(cls, parameters_str):
        if parameters_str is None:
            raise ParametersParsingException('String is empty')
        result = cls()
        if parameters_str.startswith('(') and parameters_str.endswith(')'):
            parameters_str = parameters_str[1:-1]
        for group_str in parameters_str.split(';'):
            if ',' in group_str:
                add_list = True
                for item_str in group_str.split(','):
                    if len(item_str) > 0:
                        if add_list:
                            result.entries.append(FunctionParameterList())
                            add_list = False
                        result.entries[-1].values.append(item_str)
            elif '..' in group_str:
                result.entries.append(FunctionParameterRange())
                if group_str.count('..') != 1:
                    raise ParametersParsingException()
                if group_str[-1] == 'r':
                    result.entries[-1].add_start_stop = True
                    group_str = group_str[:-1]
                try:
                    result.entries[-1].start = int(group_str.split('..')[0])
                    group_str = group_str.split('..')[1]
                    if 's' in group_str:
                        result.entries[-1].stop = int(group_str.split('s')[0])
                        result.entries[-1].step = int(group_str.split('s')[1])
                    else:
                        result.entries[-1].stop = int(group_str)
                        result.entries[-1].step = 1
                except ValueError:
                    raise ParametersParsingException('Cannot convert value to integer')
            elif len(group_str) > 0:
                result.entries.append(FunctionParameterList())
                result.entries[-1].values.append(group_str)
        return result

    @classmethod
    def verify(cls, parameters):
        try:
            cls.parse(parameters)
            return True
        except ParametersParsingException:
            return False


def main():
    raise NotImplementedError()

if __name__ == '__main__':
    main()