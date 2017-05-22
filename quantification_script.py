import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from general_utils import simulate_reads
from general_utils import get_average_accuracy
import salmon_utils as salmon
import sailfish_utils as sailfish
import kallisto_utils as kallisto


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
coverage_range = np.arange(20,30,10)
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
    
    # kallisto, max allowed k=31
    kallisto_quantificatoin_map = kallisto.run_kallisto(k, transcriptome_reference_file, get_index_dir_by_toolname("kallisto"), simulated_reads_dir, get_output_dir_by_toolname("kallisto"))
    kallisto_accuracy = get_average_accuracy(ground_truth_map, kallisto_quantificatoin_map)
    
    print("** salmon_accuracy=\t" + str(salmon_accuracy))
    print("** kallisto_accuracy=\t" + str(kallisto_accuracy))

    return salmon_accuracy, kallisto_accuracy


def plot_result(plot_name, k_range, salmon_accuracies, kallisto_accuracies):
    plt.figure()
    n_groups = len(k_range)
    index = np.arange(n_groups)

    axes = plt.gca()
    # axes.set_xlim([xmin,xmax])
    axes.set_ylim([0,1])
    bar_width = 0.2
    opacity = 0.7
    plt.bar(index, salmon_accuracies, bar_width, alpha=opacity, color='b', label='Salmon')
    plt.bar(index + bar_width, kallisto_accuracies, bar_width, alpha=opacity, color='r', label='Kallisto')

    plt.xlabel('k value')
    plt.ylabel('Accuracy')
    plt.title('Accuracy by k')

    plt.xticks(index + bar_width/2.0, k_range)
    plt.legend()
    plt.tight_layout()
    print("*** Plot save to: " + project_dir + "/" + plot_name)
    plt.savefig(project_dir + "/" + plot_name)


def run_with_simulation_parameters(number_of_transcripts, readlen, error_rate, coverage):
    print("Simulation settings:")
    print("\treadlen = " + str(readlen))
    print("\terror_rate = " + str(error_rate))
    print("\tcoverage = " + str(coverage))

    ground_truth_map = simulate_reads(simulation_script_path, number_of_transcripts, readlen, error_rate, coverage, project_dir)
    salmon_accuracies = []
    kallisto_accuracies = []
    for k in k_range:
        salmon_accuracy, kallisto_accuracy = run_with_k(k, ground_truth_map, transcriptome_reference_file, simulated_reads_dir)
        salmon_accuracies.append(salmon_accuracy)
        kallisto_accuracies.append(kallisto_accuracy)
    plot_name = "readlen"+str(readlen)+"_error"+str(error_rate).replace(".","-")+"_coverage"+str(coverage)
    
    plot_result(plot_name, k_range, salmon_accuracies, kallisto_accuracies)


def main():
    # loop coverage
    for coverage in coverage_range:
        run_with_simulation_parameters(number_of_transcripts, default_readlen, default_error_rate, coverage)

    # # loop error_rate
    # for error_rate in error_rate_range:
    #     run_with_simulation_parameters(number_of_transcripts, default_readlen, error_rate, default_coverage)

    # # loop read_len
    # for readlen in readlen_range:
    #     run_with_simulation_parameters(number_of_transcripts, readlen, default_coverage, default_coverage)


if __name__ == "__main__":
    main()
