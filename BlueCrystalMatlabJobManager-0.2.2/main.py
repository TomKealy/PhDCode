#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Startup script for BlueCrystal Matlab Job Manager"""

__author__ = "Denys Berkovskyy"
__email__ = "berkovskyy@gmail.com"
__status__ = "Development"

import sys
import interface_cli
#import interface_curses


VERSION_STRING = '0.2.2'


USAGE_STRING = '''\
Usage:
    main.py [help|cli|curses|gui]
'''


HELP_STRING = '''\
BlueCrystal Matlab Job Manager version {}
Detailed instructions are in README.md file.
'''.format(VERSION_STRING)


def main():
    if len(sys.argv) == 1:
        while True:
            options = [
                interface_cli.SelectionOption('0', 'Help', ['help']),
                interface_cli.SelectionOption('1', 'CLI', ['cli']),
                #interface_cli.SelectionOption('2', 'Curses', ['curses']),
                #interface_cli.SelectionOption('3', 'GUI', ['gui'])
            ]
            user_selection = interface_cli.get_selection('Select from following:', options)
            if user_selection == '0':
                print(USAGE_STRING)
                print(HELP_STRING)
            elif user_selection == '1':
                print('Starting CLI interface...')
                interface_cli.main()
                break
            elif user_selection == '2':
                # TODO: implement and remove
                print('Curses interface is not available')
                #print('Starting Curses interface...')
                #interface_curses.main()
                break
            elif user_selection == '3':
                print('GUI interface is not available')
            else:
                raise NotImplementedError(user_selection)
    elif len(sys.argv) == 2:
        if sys.argv[1] in ['0', 'help']:
            print(USAGE_STRING)
            print(HELP_STRING)
        elif sys.argv[1] in ['1', 'cli']:
            print('Starting CLI interface...')
            interface_cli.main()
        elif sys.argv[1] in ['2', 'curses']:
            # TODO: emplement and remove
            print('Curses interface is not available')
            #print('Starting Curses interface')
            #interface_curses.main()
        elif sys.argv[1] in ['3', 'gui']:
            print('GUI interface is not available')
        else:
            print('Unknown argument: {}'.format(sys.argv[1]))
            print(USAGE_STRING)
    else:
        print('Invalid number of arguments: {}'.format(len(sys.argv)))
        print(USAGE_STRING)


if __name__ == '__main__':
    main()
