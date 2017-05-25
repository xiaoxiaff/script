import sys
import os
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from general_utils import execute_command
from general_utils import remove_file_if_exists
from general_utils import get_average_accuracy
from general_utils import cleanup_dir
import datetime
from timeit import timeit
from data_processing import plot_accuracy_for_tool
from data_processing import plot_runtime_for_tool
from data_processing import save_result_matrix_as_csv


################### Global Vars      ############################

project_dir = ""
log_dir = ""
simulation_script_path = ""
transcriptome_reference_file = ""
simulated_reads_dir = ""

number_of_transcripts = 10
verbose = True

log_file_name = ""

################### Default Settings ############################
default_readlen = 100
default_error_rate = 0.005
default_coverage = 20


def print_and_log(line):
    print(line)

    with open(log_dir + "/" + log_file_name, "a") as log:
        log.write(line + "\n")


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

    execute_command(command, verbose)

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


def get_index_dir_by_toolname(tool_name):
    return project_dir + "/" + tool_name + "/index"


def get_output_dir_by_toolname(tool_name):
    return project_dir + "/" + tool_name + "/output"


def convert_tool_name_to_module_name(tool_name):
    return tool_name + "_utils"


def run_with_k_for_tool(tool_name, k, ground_truth_map, transcriptome_reference_file, simulated_reads_dir):
    print_and_log("quant with {0}, k={1}...".format(tool_name,str(k)))
    tool = __import__(convert_tool_name_to_module_name(tool_name), fromlist=[''])
    tool.set_verbose(verbose)
    quantificatoin_map, runtime_ms = tool.run(k, transcriptome_reference_file, get_index_dir_by_toolname(tool_name), simulated_reads_dir, get_output_dir_by_toolname(tool_name))
    accuracy = get_average_accuracy(ground_truth_map, quantificatoin_map)
    print_and_log("\t{0}_accuracy={1:f}, runtime(ms)={2:f}".format(tool_name, accuracy, runtime_ms))
    return accuracy, runtime_ms


def run_with_simulation_parameters_for_tool(tool_name, k_range, number_of_transcripts, readlen, error_rate, coverage):
    global log_file_name
    log_file_name = "{0}_readlen{1}_error{2:.0f}_coverage{3}".format(tool_name,str(readlen),error_rate*1000.0,str(coverage))
    print("\n\n")
    print_and_log("Simulation settings:")
    print_and_log('{:>30}  {:>8}'.format('number_of_transctipts:', str(number_of_transcripts)))
    print_and_log('{:>30}  {:>8}'.format('readlen:', str(readlen)))
    print_and_log('{:>30}  {:>8}'.format('error_rate:', str(error_rate)))
    print_and_log('{:>30}  {:>8}'.format('coverage:', str(coverage)))
    ground_truth_map = simulate_reads(simulation_script_path, number_of_transcripts, readlen, error_rate, coverage, project_dir)
    print_and_log('{:>30}  {:>8}'.format('Total Number of Reads:', str(sum(ground_truth_map.values()))))
    print_and_log("")

    accuracies = []
    runtimes = []
    for k in k_range:
        if tool_name=="kallisto" and k>31:
            continue
        if tool_name=="sailfish" and k>31:
            continue
        accuracy, runtime_ms = run_with_k_for_tool(tool_name, k, ground_truth_map, transcriptome_reference_file, simulated_reads_dir)
        accuracies.append(accuracy) 
        runtimes.append(runtime_ms)

    return accuracies, runtimes


def get_np_data_file_name(tool_name, loop_type, file_type):
    return "{0}/{1}/{2}_{3}_{4}_matrix".format(project_dir,'np_data',tool_name,loop_type, file_type)


def run_loop_for_tool(tool_name, loop_type, loop_range, k_range):
    accuracy_matrix = []
    runtime_matrix = []
    current_readlen = default_readlen
    current_coverage = default_coverage
    current_error_rate = default_error_rate

    for loop_value in loop_range:
        if(loop_type=="coverage"):
            current_coverage = loop_value
        elif(loop_type=="error_rate"):
            current_error_rate = loop_value
        elif(loop_type=="readlen"):
            current_readlen = loop_value

        accuracies, runtimes = run_with_simulation_parameters_for_tool(tool_name, k_range, number_of_transcripts, current_readlen, current_error_rate, current_coverage)
        accuracy_matrix.append(accuracies)
        runtime_matrix.append(runtimes)

    accuracy_np_data_file_name = get_np_data_file_name(tool_name, loop_type, 'accuracy')
    runtime_np_data_file_name = get_np_data_file_name(tool_name, loop_type, 'runtime')

    np.save(accuracy_np_data_file_name, accuracy_matrix)
    np.save(runtime_np_data_file_name, runtime_matrix)


def run_coverage_for_tool(tool_name, coverage_range, k_range):
    accuracy_matrix = []
    runtime_matrix = []

    for coverage in coverage_range:
        accuracies, runtimes = run_with_simulation_parameters_for_tool(tool_name, k_range, number_of_transcripts, default_readlen, default_error_rate, coverage)
        accuracy_matrix.append(accuracies)
        runtime_matrix.append(runtimes)

    accuracy_np_data_file = get_np_data_file_name(tool_name, 'coverage', 'accuracy')
    runtime_np_data_file = get_np_data_file_name(tool_name, 'coverage', 'runtime')

    np.save(accuracy_np_data_file,accuracy_matrix)
    np.save(runtime_np_data_file,runtime_matrix)


def run_error_rate_for_tool(tool_name, error_rate_range, k_range):
    accuracy_matrix = []
    runtime_matrix = []
    
    for error_rate in error_rate_range:
        accuracies, runtimes = run_with_simulation_parameters_for_tool(tool_name, k_range, number_of_transcripts, default_readlen, error_rate, default_coverage)
        accuracy_matrix.append(accuracies)
        runtime_matrix.append(runtimes)

    accuracy_np_data_file = get_np_data_file_name(tool_name, 'error_rate', 'accuracy')
    runtime_np_data_file = get_np_data_file_name(tool_name, 'error_rate', 'runtime')

    np.save(accuracy_np_data_file,accuracy_matrix)
    np.save(runtime_np_data_file,runtime_matrix)


def run_readlen_for_tool(tool_name, readlen_range, k_range):
    accuracy_matrix = []
    runtime_matrix = []
    
    for readlen in readlen_range:
        accuracies, runtimes = run_with_simulation_parameters_for_tool(tool_name, k_range, number_of_transcripts, readlen, default_error_rate, default_coverage)
        accuracy_matrix.append(accuracies)
        runtime_matrix.append(runtimes)

    accuracy_np_data_file = get_np_data_file_name(tool_name, 'readlen', 'accuracy')
    runtime_np_data_file = get_np_data_file_name(tool_name, 'readlen', 'runtime')

    np.save(accuracy_np_data_file,accuracy_matrix)
    np.save(runtime_np_data_file,runtime_matrix)


def init():
    parser = OptionParser()
    parser.add_option("-t", "--transcript", dest="number_of_transcripts", default=10,
        action="store", type="int",
        help="int, number of transcripts to use in the simulation")

    parser.add_option("-q", "--quiet",
        action="store_false", dest="verbose", default=True,
        help="don't print execution outputs")

    (options, args) = parser.parse_args()


    global number_of_transcripts
    global project_dir
    global simulation_script_path
    global transcriptome_reference_file
    global simulated_reads_dir
    global verbose
    global log_dir

    number_of_transcripts = options.number_of_transcripts
    verbose = options.verbose

    OS = sys.platform
    # mac
    if(OS=="darwin"):
        project_dir = "/Users/liyuanqi/Google_Drive/UCLA_MSCS/Quarter3/CS229S/Project"
        simulation_script_path = project_dir + "/simulation_script.R"
    #ubuntu
    else:
        project_dir = "/home/ubuntu/cs229"
        simulation_script_path = project_dir + "/CS229S_Project/simulation_script.R"

    transcriptome_reference_file = project_dir + "/chr22_small.fa"
    simulated_reads_dir = project_dir + "/simulated_reads"

    start_time = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
    log_dir = project_dir + "/logs/" +start_time
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    print("\nRunning on platform: " + OS + "\n")
    print("Script Settings:")
    print('{:>30}  {:<50}'.format('Verbose:', str(verbose)))
    print('{:>30}  {:<50}'.format('Project Directory:', str(project_dir)))
    print('{:>30}  {:<50}'.format('Simulation Script:', str(simulation_script_path)))
    print('{:>30}  {:<50}'.format('Transcript Reference:', str(transcriptome_reference_file)))
    print('{:>30}  {:<50}'.format('Simulated Reads Directory:', str(simulated_reads_dir)))
    print("")


def save_result(tool_name, loop_type, k_range, loop_range):
    accuracy_np_data_file_name = get_np_data_file_name(tool_name,loop_type,"accuracy")
    runtime_np_data_file = get_np_data_file_name(tool_name,loop_type,"runtime")
    accuracy_matrix = np.load(accuracy_np_data_file_name+'.npy')
    runtime_matrix = np.load(runtime_np_data_file+'.npy')
    while(len(k_range)!=len(accuracy_matrix[0])):
        k_range = np.delete(k_range,len(k_range)-1)

    save_result_matrix_as_csv(tool_name,"accuracy",loop_type,k_range,loop_range,accuracy_matrix)
    save_result_matrix_as_csv(tool_name,"runtime",loop_type,k_range,loop_range,runtime_matrix)
    plot_accuracy_for_tool(tool_name, loop_type, loop_range, k_range, accuracy_matrix)
    plot_runtime_for_tool(tool_name, loop_type, loop_range, k_range, runtime_matrix)


def main():
    init()
    ## test ranges
    k_range = np.arange(31,35,2)
    coverage_range = np.arange(20,40,10)
    error_rate_range = np.arange(0.005,0.02,0.005)
    readlen_range = np.arange(80,110,10)

    ## real ranges
    # k_range = np.arange(21,40,2)
    # coverage_range = np.arange(10,50,10)
    # error_rate_range = np.arange(0.0,0.08,0.01)
    # readlen_range = np.arange(70,130,10)

    ranges_dict = {}
    ranges_dict["coverage"] = coverage_range
    ranges_dict["error_rate"] = error_rate_range
    ranges_dict["readlen"] = readlen_range

    tools = ["rnaskim"]
    for tool_name in tools:
        for loop_type in ranges_dict.keys():
            run_loop_for_tool(tool_name, loop_type, ranges_dict[loop_type], k_range)
            save_result(tool_name,loop_type,k_range,ranges_dict[loop_type])


if __name__ == "__main__":
    main()
