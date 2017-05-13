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


def get_result_matrix(result_dir):       
    matrix = np.genfromtxt(
        result_dir + '/quant.sf',
        names = True,
        dtype=None,
        delimiter="\t")
    print matrix.dtype.names
    print matrix

# get_result_matrix('/Users/liyuanqi/Google Drive/UCLA_MSCS/Quarter3/CS229S/Project/salmon/output')
