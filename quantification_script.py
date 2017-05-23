import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from general_utils import simulate_reads
from general_utils import get_average_accuracy
from general_utils import execute_command
from general_utils import remove_file_if_exists
import salmon_utils as salmon
import sailfish_utils as sailfish
import kallisto_utils as kallisto
import RNASkim_utils as rnaskim


OS = sys.platform

print("Running on platform: " + OS)
print("Setting global variables...")

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


################### Settings ############################

k_range = np.arange(29,33,2)
coverage_range = np.arange(20,40,10)
error_rate_range = np.arange(0.0,0.1,0.005)
readlen_range = np.arange(0.0,0.1,0.005)

number_of_transcripts = 10
default_readlen = 100
default_error_rate = 0.005
default_coverage = 20

##########################################################

def get_index_dir_by_toolname(tool_name):
    return project_dir + "/" + tool_name + "/index"


def get_output_dir_by_toolname(tool_name):
    return project_dir + "/" + tool_name + "/output"


def run_with_k(k, ground_truth_map, transcriptome_reference_file, simulated_reads_dir):
    print("quant with k=" + str(k) + "...")

    # salmon
    salmon_quantificatoin_map = salmon.run_salmon(k, transcriptome_reference_file, get_index_dir_by_toolname("salmon"), simulated_reads_dir, get_output_dir_by_toolname("salmon"))
    salmon_accuracy = get_average_accuracy(ground_truth_map, salmon_quantificatoin_map)
   
    # sailfish
    sailfish_quantificatoin_map = sailfish.run_sailfish(k, transcriptome_reference_file, get_index_dir_by_toolname("sailfish"), simulated_reads_dir, get_output_dir_by_toolname("sailfish"))
    sailfish_accuracy = get_average_accuracy(ground_truth_map, sailfish_quantificatoin_map)
       
    # kallisto, max allowed k=31
    kallisto_quantificatoin_map = kallisto.run_kallisto(k, transcriptome_reference_file, get_index_dir_by_toolname("kallisto"), simulated_reads_dir, get_output_dir_by_toolname("kallisto"))
    kallisto_accuracy = get_average_accuracy(ground_truth_map, kallisto_quantificatoin_map)
    
    # rnaskim
    rnaskim_quantificatoin_map = rnaskim.run_RNASkim(k, transcriptome_reference_file, get_index_dir_by_toolname("rnaskim"), simulated_reads_dir, get_output_dir_by_toolname("rnaskim"), 4)
    rnaskim_accuracy = get_average_accuracy(ground_truth_map, rnaskim_quantificatoin_map)
  
    print("** salmon_accuracy=\t" + str(salmon_accuracy))
    print("** sailfish_accuracy=\t" + str(sailfish_accuracy))
    print("** kallisto_accuracy=\t" + str(kallisto_accuracy))
    print("** rnaskim_accuracy=\t" + str(rnaskim_accuracy))

    return salmon_accuracy, sailfish_accuracy, kallisto_accuracy, rnaskim_accuracy


def plot_result_all(readlen, error_rate, coverage, k_range, salmon_accuracies, sailfish_accuracies, kallisto_accuracies, rnaskim_accuracies):
    plt.figure()
    n_groups = len(k_range)
    index = np.arange(n_groups)

    axes = plt.gca()
    # axes.set_xlim([xmin,xmax])
    # axes.set_ylim([0,1])
    salmon_label = "Salmon,bestk={0:d},accuracy={1:.4f}".format(k_range[np.argmax(salmon_accuracies)],max(salmon_accuracies))
    sailfish_label = "Sailfish,bestk={0:d},accuracy={1:.4f}".format(k_range[np.argmax(sailfish_accuracies)],max(sailfish_accuracies))
    kallisto_label = "Kallisto,bestk={0:d},accuracy={1:.4f}".format(k_range[np.argmax(kallisto_accuracies)],max(kallisto_accuracies))
    rnaskim_label = "RNASkim,bestk={0:d},accuracy={1:.4f}".format(k_range[np.argmax(rnaskim_accuracies)],max(rnaskim_accuracies))


    plt.plot(index, salmon_accuracies, color='b', label=salmon_label)
    plt.plot(index, sailfish_accuracies, color='r', label=sailfish_label)
    plt.plot(index, kallisto_accuracies, color='r', label=kallisto_label)
    plt.plot(index, rnaskim_accuracies, color='r', label=rnaskim_label)


    plt.xlabel('k value')
    plt.ylabel('Accuracy')
    title = "readlen={0}|error_rate={1:.2f}%|coverage={2}".format(str(readlen),error_rate*100.0,str(coverage))
    plt.title(title)

    plot_name = "all_readlen{0}_error_rate{1:.0f}_coverage{2}".format(str(readlen),error_rate*1000.0,str(coverage))
    plt.xticks(index, k_range)
    plt.legend()
    plt.savefig(project_dir + "/" + plot_name)


def plot_result_line_for_tool(tool_name, plot_type, labels, k_range, accuracy_matrix):
    colors = ['r','g','b','y','m','c','k']
    plt.figure()
    n_groups = len(k_range)
    index = np.arange(n_groups)
    axes = plt.gca()
    # axes.set_xlim([xmin,xmax])
    # axes.set_ylim([0,1])
    plt.xlabel('k value')
    plt.ylabel('Accuracy')
    plot_title = "{0} {1} plot".format(tool_name, plot_type)
    plt.title(plot_title)

    for i in range(0,len(accuracy_matrix)):
        current_array = accuracy_matrix[i]
        besk_k = k_range[np.argmax(current_array)]
        best_accuracy = max(current_array)
        label = "{0}={1},bestK={2},accuracy={3:.4f}%".format(plot_type,str(labels[i]),str(besk_k),best_accuracy)
        plt.plot(index, current_array, colors[i]+'o-', label=label)
    
    plt.xticks(index, k_range)
    plt.legend()
    plt.savefig(project_dir + "/" + tool_name + "_" + plot_type + "_plot")


def run_with_simulation_parameters(number_of_transcripts, readlen, error_rate, coverage):
    print("Simulation settings:")
    print("\treadlen = " + str(readlen))
    print("\terror_rate = " + str(error_rate))
    print("\tcoverage = " + str(coverage))

    ground_truth_map = simulate_reads(simulation_script_path, number_of_transcripts, readlen, error_rate, coverage, project_dir)
    salmon_accuracies = []
    sailfish_accuracies = []
    kallisto_accuracies = []
    rnaskim_accuracies = []
    for k in k_range:
        salmon_accuracy, sailfish_accuracy, kallisto_accuracy, rnaskim_accuracy = run_with_k(k, ground_truth_map, transcriptome_reference_file, simulated_reads_dir)
        salmon_accuracies.append(salmon_accuracy)
        sailfish_accuracies.append(sailfish_accuracy)
        kallisto_accuracies.append(kallisto_accuracy)  
        rnaskim_accuracies.append(rnaskim_accuracy)  
    plot_result_all(readlen, error_rate, coverage, k_range, salmon_accuracies, sailfish_accuracies, kallisto_accuracies, rnaskim_accuracies)

    return salmon_accuracies, sailfish_accuracies, kallisto_accuracies, rnaskim_accuracies


def main():
    # loop coverage
    salmon_accuracy_matrix = []
    kallisto_accuracy_matrix = []
    sailfish_accuracy_matrix = []
    rnaskim_accuracy_matrix = []

    for coverage in coverage_range:
        salmon_accuracies, sailfish_accuracies, kallisto_accuracies, rnaskim_accuracies = run_with_simulation_parameters(number_of_transcripts, default_readlen, default_error_rate, coverage)
        salmon_accuracy_matrix.append(salmon_accuracies)
        sailfish_accuracy_matrix.append(sailfish_accuracies)
        kallisto_accuracy_matrix.append(kallisto_accuracies)
        rnaskim_accuracy_matrix.append(rnaskim_accuracies)
    plot_result_line_for_tool("salmon", "coverage", coverage_range, k_range, salmon_accuracy_matrix)
    plot_result_line_for_tool("sailfish", "coverage", coverage_range, k_range, sailfish_accuracy_matrix)
    plot_result_line_for_tool("kallisto", "coverage", coverage_range, k_range, kallisto_accuracy_matrix)
    plot_result_line_for_tool("rnaskim", "coverage", coverage_range, k_range, rnaskim_accuracy_matrix)


    # # loop error_rate
    # for error_rate in error_rate_range:
    #     run_with_simulation_parameters(number_of_transcripts, default_readlen, error_rate, default_coverage)

    # # loop read_len
    # for readlen in readlen_range:
    #     run_with_simulation_parameters(number_of_transcripts, readlen, default_coverage, default_coverage)


if __name__ == "__main__":
    main()
