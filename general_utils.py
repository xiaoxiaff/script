import re
import sys
import os
import subprocess as sub


def check_command_output(output, command_string):
    if re.search(r'\berror\b', output, re.I):
        print("\n\nError in executing command:")
        print(command_string)
        print("Exiting...\n\n")
        sys.exit()


def execute_command(command_string, verbose):
    command_args = command_string.split()
    print(command_args)
    if(verbose):
        print("\nExecuting command:")
        print(command_string + "\n\n")
        sub.check_call(command_args)
    else:
	    FNULL = open(os.devnull, 'w')
	    sub.check_call(command_args, stdout=FNULL, stderr=sub.STDOUT)

