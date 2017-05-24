from general_utils import execute_command
import numpy as np
import time
# https://pachterlab.github.io/kallisto/manual


verbose = False
def set_verbose(v):
    global verbose
    verbose = v


# k: default and maximum k is 31, must be odd
def build_index_with_k(transcriptome_reference_file, k, index_output_path):
	command = "kallisto index -i " \
		+ index_output_path \
		+ " -k " \
		+ str(k) \
		+ " " \
		+ transcriptome_reference_file
	
	execute_command(command, verbose)


def quant(index_dir, output_dir, sample_pair1, sample_pair2):
	command = "kallisto quant -i " \
		+ index_dir \
		+ " -o " \
		+ output_dir \
		+ " " \
		+ sample_pair1 \
		+ " " \
		+ sample_pair2 \

	execute_command(command, verbose)


def get_result_dict(result_dir):       
    matrix = np.genfromtxt(
        result_dir + '/abundance.tsv',
        names = True,
        dtype=None,
        delimiter="\t")

    transcripts = matrix['target_id']
    num_reads = matrix['est_counts']
    res = dict()

    for i in range(0, len(transcripts)):
        res[transcripts[i]] = num_reads[i]

    return res


def run(k, transcriptome_reference_file, index_output_dir, sample_dir, output_dir):
    time1 = time.time()

    build_index_with_k(transcriptome_reference_file, k, index_output_dir)

    res_dict = dict()
    for i in range(1, 11):
        if i < 10:
            sample_pair1 = sample_dir + '/sample_0' + str(i)+ '_1.fasta'
            sample_pair2 = sample_dir + '/sample_0' + str(i) + '_2.fasta'
        else:
            sample_pair1 = sample_dir + '/sample_' + str(i)+ '_1.fasta'
            sample_pair2 = sample_dir + '/sample_' + str(i) + '_2.fasta'

        result_output_dir = output_dir + '/sample' + str(i) + '_result'
        quant(index_output_dir, result_output_dir, sample_pair1, sample_pair2)
        sample_dict = get_result_dict(result_output_dir)

        for key in sample_dict:
            if i == 1:
                res_dict[key] = sample_dict[key]
            else:
                res_dict[key] += sample_dict[key]

    for key in res_dict:
        res_dict[key] /= 10
  
    time2 = time.time()
    elapsed_ms = (time2-time1)*1000.0
    
    return res_dict, elapsed_ms


# run_kallisto(31, "chr22_small.fa", "test_results/kallisto/index", "simulated_reads", "test_results/kallisto/quant_results")


