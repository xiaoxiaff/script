import re
import sys
import os
import subprocess as sub
import numpy as np


def remove_file_if_exists(file_path):
    try:
        os.remove(file_path)
    except OSError:
        pass


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
    
    remove_file_if_exists(output_dir + '/transcript_names.txt')
    remove_file_if_exists(output_dir + '/num_of_reads.txt')

    command = "Rscript --vanilla " \
        + script_path + " " \
        + str(number_of_transcripts) + " " \
        + str(readlen) + " " \
        + str(error_rate) + " " \
        + str(coverage) + " " \
        + str(output_dir) + " " \

    execute_command(command, True)

    transcript_names = np.genfromtxt(
        output_dir + '/transcript_names.txt',
        names = None,
        dtype= None,
        usecols = (0))
    num_of_reads = np.genfromtxt(
        output_dir + '/num_of_reads.txt',
        names = None,
        dtype= None,
        usecols = (0))

    ground_truth_map = dict(zip(transcript_names, num_of_reads))

    return ground_truth_map
