from general_utils import execute_command


def build_index(k, transcriptome_reference_file, index_output_path):
    command = "sailfish index -t " \
        + transcriptome_reference_file \
        + " -i " \
        + index_output_path \
        + " -k " \
        + str(k)

    execute_command(command,True)


def quant_with_k(k, sample_pair1,sample_pair2, index_dir, output_dir, transcriptome_reference_file, index_output_path):
    build_index(k, transcriptome_reference_file, index_output_path)
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
