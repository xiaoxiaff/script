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


# def run_with_k_for_all_for_all(k, ground_truth_map, transcriptome_reference_file, simulated_reads_dir):
#     print_and_log("quant with k=" + str(k) + "...")

#     # salmon
#     print_and_log("run salmon...")
#     salmon_quantificatoin_map = salmon.run_salmon(k, transcriptome_reference_file, get_index_dir_by_toolname("salmon"), simulated_reads_dir, get_output_dir_by_toolname("salmon"))
#     salmon_accuracy = get_average_accuracy(ground_truth_map, salmon_quantificatoin_map)
   
#     # sailfish
#     print_and_log("run sailfish...")
#     sailfish_quantificatoin_map = sailfish.run_sailfish(k, transcriptome_reference_file, get_index_dir_by_toolname("sailfish"), simulated_reads_dir, get_output_dir_by_toolname("sailfish"))
#     sailfish_accuracy = get_average_accuracy(ground_truth_map, sailfish_quantificatoin_map)
       
#     # kallisto, max allowed k=31
#     print_and_log("run kallisto...")
#     kallisto_quantificatoin_map = kallisto.run_kallisto(k, transcriptome_reference_file, get_index_dir_by_toolname("kallisto"), simulated_reads_dir, get_output_dir_by_toolname("kallisto"))
#     kallisto_accuracy = get_average_accuracy(ground_truth_map, kallisto_quantificatoin_map)
    
#     # rnaskim
#     print_and_log("run rnaskim...")
#     rnaskim_quantificatoin_map = rnaskim.run_RNASkim(k, transcriptome_reference_file, get_index_dir_by_toolname("rnaskim"), simulated_reads_dir, get_output_dir_by_toolname("rnaskim"), 4)
#     rnaskim_accuracy = get_average_accuracy(ground_truth_map, rnaskim_quantificatoin_map)
  
#     print_and_log("\tsalmon_accuracy=\t" + str(salmon_accuracy))
#     print_and_log("\tsailfish_accuracy=\t" + str(sailfish_accuracy))
#     print_and_log("\tkallisto_accuracy=\t" + str(kallisto_accuracy))
#     print_and_log("\trnaskim_accuracy=\t" + str(rnaskim_accuracy))

#     return salmon_accuracy, sailfish_accuracy, kallisto_accuracy, rnaskim_accuracy


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


# def plot_result_all(readlen, error_rate, coverage, k_range, salmon_accuracies, sailfish_accuracies, kallisto_accuracies, rnaskim_accuracies):
#     plt.figure()
#     n_groups = len(k_range)
#     index = np.arange(n_groups)

#     axes = plt.gca()
#     # axes.set_xlim([xmin,xmax])
#     # axes.set_ylim([0,1])
#     salmon_label = "Salmon,bestk={0:d},accuracy={1:.4f}".format(k_range[np.argmax(salmon_accuracies)],max(salmon_accuracies))
#     sailfish_label = "Sailfish,bestk={0:d},accuracy={1:.4f}".format(k_range[np.argmax(sailfish_accuracies)],max(sailfish_accuracies))
#     kallisto_label = "Kallisto,bestk={0:d},accuracy={1:.4f}".format(k_range[np.argmax(kallisto_accuracies)],max(kallisto_accuracies))
#     rnaskim_label = "RNASkim,bestk={0:d},accuracy={1:.4f}".format(k_range[np.argmax(rnaskim_accuracies)],max(rnaskim_accuracies))


#     plt.plot(index, salmon_accuracies, color='b', label=salmon_label)
#     plt.plot(index, sailfish_accuracies, color='r', label=sailfish_label)
#     plt.plot(index, kallisto_accuracies, color='g', label=kallisto_label)
#     plt.plot(index, rnaskim_accuracies, color='y', label=rnaskim_label)


#     plt.xlabel('k value')
#     plt.ylabel('Accuracy')
#     title = "readlen={0}|error_rate={1:.2f}%|coverage={2}".format(str(readlen),error_rate*100.0,str(coverage))
#     plt.title(title)

#     plot_name = "all_readlen{0}_error_rate{1:.0f}_coverage{2}".format(str(readlen),error_rate*1000.0,str(coverage))
#     plt.xticks(index, k_range)
#     plt.legend()
#     plt.savefig(project_dir + "/" + plot_name)


def plot_accuracy_for_tool(tool_name, plot_type, labels, k_range, accuracy_matrix):
    colors = ['r','g','b','y','m','c','k', 'pink']
    plt.figure()
    n_groups = len(k_range)
    index = np.arange(n_groups)
    axes = plt.gca()
    # axes.set_xlim([xmin,xmax])
    # axes.set_ylim([0,1])
    plt.xlabel('k value')
    plt.ylabel('Accuracy')
    plot_title = "{0} {1} accuracy plot".format(tool_name, plot_type)
    plt.title(plot_title)

    for i in range(0,len(accuracy_matrix)):
        current_array = accuracy_matrix[i]
        besk_k = k_range[np.argmax(current_array)]
        best_accuracy = max(current_array)
        label = "{0}={1},bestK={2},accuracy={3:.4f}".format(plot_type,str(labels[i]),str(besk_k),best_accuracy)
        plt.plot(index, current_array, color = colors[i], ls='-', marker='o', label=label)
    
    plt.xticks(index, k_range)
    plt.legend()
    plt.savefig(project_dir + "/" + tool_name + "_" + plot_type + "_accuracy_plot")


def plot_runtime_for_tool(tool_name, plot_type, labels, k_range, runtime_matrix):
    colors = ['r','g','b','y','m','c','k', 'pink']
    plt.figure()
    n_groups = len(k_range)
    index = np.arange(n_groups)
    axes = plt.gca()
    # axes.set_xlim([xmin,xmax])
    # axes.set_ylim([0,1])
    plt.xlabel('k value')
    plt.ylabel('Runtime (ms)')
    plot_title = "{0} {1} runtime plot".format(tool_name, plot_type)
    plt.title(plot_title)

    for i in range(0,len(runtime_matrix)):
        current_array = runtime_matrix[i]
        label = "{0}={1}".format(plot_type,str(labels[i]))
        plt.plot(index, current_array, color = colors[i], ls='-', marker='o', label=label)
    
    plt.xticks(index, k_range)
    plt.legend()
    plt.savefig(project_dir + "/" + tool_name + "_" + plot_type + "_runtime_plot")


# def run_with_simulation_parameters_for_all(number_of_transcripts, readlen, error_rate, coverage):
#     print_and_log("Simulation settings:")
#     print_and_log('{:>50}  {:>12}'.format('number_of_transctipts:', str(number_of_transcripts)))
#     print_and_log('{:>50}  {:>12}'.format('readlen:', str(readlen)))
#     print_and_log('{:>50}  {:>12}'.format('error_rate:', str(error_rate)))
#     print_and_log('{:>50}  {:>12}'.format('coverage:', str(coverage)))
#     ground_truth_map = simulate_reads(simulation_script_path, number_of_transcripts, readlen, error_rate, coverage, project_dir)
#     print_and_log('{:>50}  {:>12}'.format('Total Number of Reads:', str(sum(ground_truth_map.values()))))
#     print_and_log("")

#     salmon_accuracies = []
#     sailfish_accuracies = []
#     kallisto_accuracies = []
#     rnaskim_accuracies = []

#     for k in k_range:
#         salmon_accuracy, sailfish_accuracy, kallisto_accuracy, rnaskim_accuracy = run_with_k_for_all(k, ground_truth_map, transcriptome_reference_file, simulated_reads_dir)
#         salmon_accuracies.append(salmon_accuracy)
#         sailfish_accuracies.append(sailfish_accuracy)
#         kallisto_accuracies.append(kallisto_accuracy)  
#         rnaskim_accuracies.append(rnaskim_accuracy)  
#     plot_result_all(readlen, error_rate, coverage, k_range, salmon_accuracies, sailfish_accuracies, kallisto_accuracies, rnaskim_accuracies)

#     return salmon_accuracies, sailfish_accuracies, kallisto_accuracies, rnaskim_accuracies


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
        accuracy, runtime_ms = run_with_k_for_tool(tool_name, k, ground_truth_map, transcriptome_reference_file, simulated_reads_dir)
        accuracies.append(accuracy) 
        runtimes.append(runtime_ms)
    # plot_result_all(readlen, error_rate, coverage, k_range, salmon_accuracies, sailfish_accuracies, kallisto_accuracies, rnaskim_accuracies)

    return accuracies, runtimes


def run_coverage_for_tool(tool_name, coverage_range, k_range):
    accuracy_matrix = []
    runtime_matrix = []

    for coverage in coverage_range:
        accuracies, runtimes = run_with_simulation_parameters_for_tool(tool_name, k_range, number_of_transcripts, default_readlen, default_error_rate, coverage)
        accuracy_matrix.append(accuracies)
        runtime_matrix.append(runtimes)
    plot_accuracy_for_tool(tool_name, "coverage", coverage_range, k_range, accuracy_matrix)
    plot_runtime_for_tool(tool_name, "coverage", coverage_range, k_range, runtime_matrix)


def run_error_rate_for_tool(tool_name, error_rate_range, k_range):
    accuracy_matrix = []
    runtime_matrix = []
    
    for error_rate in error_rate_range:
        accuracies, runtimes = run_with_simulation_parameters_for_tool(tool_name, k_range, number_of_transcripts, default_readlen, error_rate, default_coverage)
        accuracy_matrix.append(accuracies)
        runtime_matrix.append(runtimes)
    plot_accuracy_for_tool(tool_name, "error_rate", error_rate_range, k_range, accuracy_matrix)
    plot_runtime_for_tool(tool_name, "error_rate", error_rate_range, k_range, runtime_matrix)


def run_readlen_for_tool(tool_name, readlen_range, k_range):
    accuracy_matrix = []
    runtime_matrix = []
    
    for readlen in readlen_range:
        accuracies, runtimes = run_with_simulation_parameters_for_tool(tool_name, k_range, number_of_transcripts, readlen, default_error_rate, default_coverage)
        accuracy_matrix.append(accuracies)
        runtime_matrix.append(runtimes)
    plot_accuracy_for_tool(tool_name, "readlen", readlen_range, k_range, accuracy_matrix)
    plot_runtime_for_tool(tool_name, "readlen", readlen_range, k_range, runtime_matrix)


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


def main():
    init()
    ## test ranges
    # k_range = np.arange(29,33,2)
    # coverage_range = np.arange(20,40,10)
    # error_rate_range = np.arange(0.005,0.02,0.005)
    # readlen_range = np.arange(80,110,10)

    ## real ranges
    k_range = np.arange(21,32,2)
    coverage_range = np.arange(10,50,10)
    error_rate_range = np.arange(0.0,0.08,0.01)
    readlen_range = np.arange(70,130,10)

    tools = ["salmon","sailfish","kallisto"]
    for tool_name in tools:
        run_coverage_for_tool(tool_name, coverage_range, k_range)
        run_error_rate_for_tool(tool_name, error_rate_range, k_range)
        run_readlen_for_tool(tool_name, readlen_range, k_range)


if __name__ == "__main__":
    main()
