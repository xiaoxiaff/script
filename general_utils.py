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
    if(verbose):
        print("\nExecuting command:")
        print(command_string + "\n\n")
        sub.check_call(command_args)
    else:
	    FNULL = open(os.devnull, 'w')
	    sub.check_call(command_args, stdout=FNULL, stderr=sub.STDOUT)


def simulate_reads(
    script_path, 
    number_of_transcripts, 
    readlen, 
    error_rate, 
    coverage, 
    output_dir):
    
    command = "Rscript --vanilla " \
        + script_path + " " \
        + str(number_of_transcripts) + " " \
        + str(readlen) + " " \
        + str(error_rate) + " " \
        + str(coverage) + " " \
        + str(output_dir) + " " \

    execute_command(command, True)
