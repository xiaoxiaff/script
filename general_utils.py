import re
import sys
import os
import subprocess as sub
import numpy as np
from sklearn.metrics import average_precision_score


def remove_file_if_exists(file_path):
    try:
        os.remove(file_path)
    except OSError:
        pass


def execute_command(command_string, verbose):
    command_args = command_string.split()
    try:
        FNULL = open(os.devnull, 'w')
        sub.check_call(command_args, stdout=FNULL, stderr=sub.STDOUT)
    except:
        print("\n*** Error while executing:\n" + command_string + "\n")
        raise


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


def get_average_percentage_error(ground_truth_map, quantificatoin_map):
    # print(ground_truth_map)
    # print(quantificatoin_map)

    errors = []

    for (key,ground_truth_value) in ground_truth_map.items():
        quantification_value = quantificatoin_map[key]
        error = abs(float(quantification_value) - float(ground_truth_value))/float(ground_truth_value)
        errors.append(error)

    return np.mean(errors)










