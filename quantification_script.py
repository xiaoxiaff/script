from general_utils import simulate_reads
# import salmon_utils as salmon
import sailfish_utils as sailfish


OS = "mac"

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

    print(ground_truth_map)
        # print("Building salmon index...")
    # salmon.build_index(transcriptome_reference_file, salmon_index_dir)


    # print("Quant with salmon, k=31...")
    # salmon.quant_with_k(
    #     31,
    #     simulated_reads_dir+"/sample_01_1.fasta",
    #     simulated_reads_dir+"/sample_01_2.fasta",
    #     salmon_index_dir,
    #     salmon_output_dir
    #     )


    # print("Quant with sailfish, k=31...")
    # sailfish.quant_with_k(
    #     31,
    #     simulated_reads_dir+"/sample_01_1.fasta",
    #     simulated_reads_dir+"/sample_01_2.fasta",
    #     sailfish_index_dir,
    #     sailfish_output_dir,
    #     transcriptome_reference_file,
    #     sailfish_index_dir
    #     )


if __name__ == "__main__":
    main()
