import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from general_utils import simulate_reads
from general_utils import get_average_percentage_error
import salmon_utils as salmon
import sailfish_utils as sailfish
import kallisto_utils as kallisto


OS = "ubuntu"

if(OS=="ubuntu"):
    project_dir = "/home/ubuntu/cs229"
    simulation_script_path = project_dir + "/CS229S_Project/simulation_script.R"
else:
    project_dir = "/Users/liyuanqi/Google_Drive/UCLA_MSCS/Quarter3/CS229S/Project"
    simulation_script_path = project_dir + "/simulation_script.R"

####### Universal Settings ###########################
transcriptome_reference_file = project_dir + "/chr22_small.fa"
simulated_reads_dir = project_dir + "/simulated_reads"


####### Salmon Settings #############################
salmon_index_dir = project_dir + "/salmon/index"
salmon_output_dir = project_dir + "/salmon/output"


####### Sailfish Settings ###########################
sailfish_index_dir = project_dir + "/sailfish/index"
sailfish_output_dir = project_dir + "/sailfish/output"


####### Kallisto Settings ###########################
kallisto_index_dir = project_dir + "/kallisto/index"
kallisto_output_dir = project_dir + "/kallisto/output"



def run_with_k(k, ground_truth_map, transcriptome_reference_file, simulated_reads_dir):
    print("quant with k=" + str(k) + "...")
    # salmon
    salmon_quantificatoin_map = salmon.run_salmon(k, transcriptome_reference_file, salmon_index_dir, simulated_reads_dir, salmon_output_dir)
    salmon_error = get_average_percentage_error(ground_truth_map, salmon_quantificatoin_map)
    
    # kallisto, max allowed k=31
    kallisto_quantificatoin_map = kallisto.run_kallisto(k, transcriptome_reference_file, kallisto_index_dir, simulated_reads_dir, kallisto_output_dir)
    kallisto_error = get_average_percentage_error(ground_truth_map, kallisto_quantificatoin_map)
    
    return salmon_error, kallisto_error


def plot_result(k_range, salmon_errors, kallisto_errors):
    n_groups = len(k_range)
    index = np.arange(n_groups)

    axes = plt.gca()
    # axes.set_xlim([xmin,xmax])
    axes.set_ylim([0,1])
    bar_width = 0.2
    opacity = 0.7
    plt.bar(index, salmon_errors, bar_width, alpha=opacity, color='b', label='Salmon')
    plt.bar(index + bar_width, kallisto_errors, bar_width, alpha=opacity, color='r', label='Kallisto')

    plt.xlabel('k value')
    plt.ylabel('Error')
    plt.title('Error by k')

    plt.xticks(index + bar_width/2.0, k_range)
    plt.legend()
    plt.tight_layout()
    plt.savefig(project_dir + "/k")


def main():
    number_of_transcripts = 10
    readlen = 100
    error_rate = 0.001
    coverage = 20
    output_dir = project_dir

    ground_truth_map = simulate_reads(
        simulation_script_path,
        number_of_transcripts,
        readlen,
        error_rate,
        coverage,
        output_dir)

    # sailfish.quant_with_k(31, 
    #     simulated_reads_dir + "/sample_01_1.fasta",
    #     simulated_reads_dir + "/sample_01_2.fasta",
    #     sailfish_index_dir,
    #     sailfish_output_dir,
    #     transcriptome_reference_file)

    k_range = range(11,14,2)
    salmon_errors = []
    kallisto_errors = []
    for k in k_range:
        salmon_error, kallisto_error = run_with_k(k, ground_truth_map, transcriptome_reference_file, simulated_reads_dir)
        salmon_errors.append(salmon_error)
        kallisto_errors.append(kallisto_error)

    plot_result(k_range, salmon_errors, kallisto_errors)

if __name__ == "__main__":
    main()
