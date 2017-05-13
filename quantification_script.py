import salmon_utils as salmon
import sailfish_utils as sailfish


####### Universal Settings ###########################
transcriptome_reference_file = "/home/ubuntu/cs229/chr22_small.fa"
simulated_reads_dir = "/home/ubuntu/cs229/simulated_reads"


####### Salmon Settings #############################
salmon_index_dir = "/home/ubuntu/cs229/salmon/index"
salmon_output_dir = "/home/ubuntu/cs229/salmon/output"


####### Sailfish Settings ###########################
sailfish_index_dir = "/home/ubuntu/cs229/sailfish/index"
sailfish_output_dir = "/home/ubuntu/cs229/sailfish/output"


def main():
    print("Building salmon index...")
    salmon.build_index(transcriptome_reference_file, salmon_index_dir)


    print("Quant with salmon, k=31...")
    salmon.quant_with_k(
        31,
        simulated_reads_dir+"/sample_01_1.fasta",
        simulated_reads_dir+"/sample_01_2.fasta",
        salmon_index_dir,
        salmon_output_dir
        )


    print("Quant with sailfish, k=31...")
    sailfish.quant_with_k(
        31,
        simulated_reads_dir+"/sample_01_1.fasta",
        simulated_reads_dir+"/sample_01_2.fasta",
        sailfish_index_dir,
        sailfish_output_dir
        )


if __name__ == "__main__":
    main()
