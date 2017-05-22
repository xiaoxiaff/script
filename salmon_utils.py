from general_utils import execute_command
import numpy as np


def build_index(transcriptome_reference_file,index_output_path):
    command = "salmon index -t " \
        + transcriptome_reference_file \
        + " -i " \
        + index_output_path
        
    execute_command(command,True)


def quant_with_k(k,sample_pair1,sample_pair2, index_dir, output_dir):
    command = "salmon quant -i " \
        + index_dir \
        + " -l A" \
        + " -1 " \
        + sample_pair1 \
        + " -2 " \
        + sample_pair2 \
        + " -p 8 -o " \
        + output_dir \
        + " -k " \
        + str(k)

    execute_command(command,True)


def get_result_dict(result_dir):       
    matrix = np.genfromtxt(
        result_dir + '/quant.sf',
        names = True,
        dtype=None,
        delimiter="\t")
<<<<<<< HEAD
    
    # print(matrix.dtype.names)

    transcripts = matrix['Name']
    num_reads = matrix['NumReads']
    res = dict()

    for i in range(0, len(transcripts)):
        res[transcripts[i]] = num_reads[i]

    return res


def run_salmon(k, transcriptome_reference_file, index_output_dir, sample_dir, output_dir):
    # build_index(transcriptome_reference_file, index_output_dir)

    res_dict = dict()
    for i in range(1, 11):
        if i < 10:
            sample_pair1 = sample_dir + '/sample_0' + str(i)+ '_1.fasta'
            sample_pair2 = sample_dir + '/sample_0' + str(i) + '_2.fasta'
        else:
            sample_pair1 = sample_dir + '/sample_' + str(i)+ '_1.fasta'
            sample_pair2 = sample_dir + '/sample_' + str(i) + '_2.fasta'

        result_output_dir = output_dir + '/sample' + str(i) + '_result'
        # quant_with_k(k, sample_pair1, sample_pair2, index_output_dir, result_output_dir)
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





run_salmon(31, "chr22_small.fa", "test_results/salmon/index", "simulated_reads", "test_results/salmon/quant_results")


# build_index("chr22_small.fa", "test_results/salmon")
# quant_with_k(31, "simulated_reads/sample_01_1.fasta", "simulated_reads/sample_01_2.fasta", "test_results/salmon", "test_results/salmon")
# get_result_matrix("test_results/salmon")


=======
    # matrix[matrix[:,1].argsort()]
    # matrix.sort(axis='')
    print matrix.dtype.names
    print matrix



# get_result_matrix('/Users/liyuanqi/Google Drive/UCLA_MSCS/Quarter3/CS229S/Project/salmon/output')
>>>>>>> 170926fefb5d6e86e567c807d1fca65047e757b1
