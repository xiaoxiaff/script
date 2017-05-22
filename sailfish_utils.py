from general_utils import execute_command
import numpy as np


def build_index(k, transcriptome_reference_file, index_output_path):
    command = "sailfish index -t " \
        + transcriptome_reference_file \
        + " -o " \
        + index_output_path \
        + " -k " \
        + str(k)

    execute_command(command,True)


def quant_with_k(sample_pair1, sample_pair2, index_dir, output_dir):
    command = "sailfish quant -i " \
        + index_dir \
        + " -l " \
        + "T=PE:O=><:S=SA" \
        + " -1 " \
        + sample_pair1 \
        + " -2 " \
        + sample_pair2 \
        + " -p 8 -o " \
        + output_dir

    execute_command(command,True)

def get_result_dict(result_dir):       
    matrix = np.genfromtxt(
        result_dir + '/quant.sf',
        names = True,
        dtype=None,
        delimiter="\t",
        skip_header=4)

    # print(matrix.dtype.names)

    transcripts = matrix['Transcript']
    num_reads = matrix['EstimatedNumReads']
    res = dict()

    for i in range(0, len(transcripts)):
        res[transcripts[i]] = num_reads[i]

    return res


def run_sailfish(k, transcriptome_reference_file, index_dir, sample_dir, output_dir):
    build_index(k, transcriptome_reference_file, index_dir)
    
    res_dict = dict()
    for i in range(1, 11):
        if i < 10:
            sample_pair1 = sample_dir + '/sample_0' + str(i)+ '_1.fasta'
            sample_pair2 = sample_dir + '/sample_0' + str(i) + '_2.fasta'
        else:
            sample_pair1 = sample_dir + '/sample_' + str(i)+ '_1.fasta'
            sample_pair2 = sample_dir + '/sample_' + str(i) + '_2.fasta'

        result_output_dir = output_dir + '/sample' + str(i) + '_result'
        quant_with_k(sample_pair1, sample_pair2, index_dir, result_output_dir)

        sample_dict = get_result_dict(result_output_dir)

        for key in sample_dict:
            if i == 1:
                res_dict[key] = sample_dict[key]
            else:
                res_dict[key] += sample_dict[key]

    for key in res_dict:
        res_dict[key] /= 10
          

    print(res_dict)
    return res_dict


run_sailfish(31, "chr22_small.fa", "test_results/sailfish/index", "simulated_reads", "test_results/sailfish/results")

