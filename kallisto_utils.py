from general_utils import execute_command
# https://pachterlab.github.io/kallisto/manual

# k: default and maximum k is 31, must be odd
def build_index_with_k(transcriptome_reference_file, k, index_output_path):
	command = "kallisto index -i " \
		+ index_output_path \
		+ " -k " \
		+ str(k) \
		+ " " \
		+ transcriptome_reference_file
	
	execute_command(command, True)


def quant(index_dir, output_dir, sample_pair1, sample_pair2):
	command = "kallisto quant -i " \
		+ index_dir \
		+ " -o " \
		+ output_dir \
		+ " " \
		+ sample_pair1 \
		+ " " \
		+ sample_pair2 \

	execute_command(command, True)



build_index_with_k("chr22_small.fa", 29, "kallisto/index2")
quant("kallisto/index2", "kallisto/output2", "simulated_reads/sample_01_1.fasta", "simulated_reads/sample_01_2.fasta")