from general_utils import execute_command


lib_dir = "LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ubuntu/cs229/sailfish/Sailfish-0.6.3-Linux_x86-64/lib"

def build_index(k, transcriptome_reference_file, index_output_path):
    command = lib_dir \
        + "sailfish index -t " \
        + transcriptome_reference_file \
        + " -i " \
        + index_output_path \
        + " -k " \
        + str(k)

    execute_command(command,True)


def quant_with_k(k, sample_pair1,sample_pair2, index_dir, output_dir, transcriptome_reference_file, index_output_path):
    build_index(k, transcriptome_reference_file, index_output_path)
    command = lib_dir \
        + "sailfish quant -i " \
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
