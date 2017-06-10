from general_utils import execute_command
from general_utils import get_command_output
from general_utils import cleanup_dir
import numpy as np
import time


verbose = False
def set_verbose(v):
    global verbose
    verbose = v


# GLOG_logtostderr=1 ./rs_cluster  -gene_fasta=gene.fa -num_threads=4 -output=clustered.fa -rs_length=60
def cluster(transcriptome_reference_file, cluster_output, k, numthreads):
    command = "rs_cluster -gene_fasta=" \
        + transcriptome_reference_file \
        + " -num_threads=" \
        + str(numthreads) \
        + " -output=" \
        + cluster_output \
        + " -rs_length=" \
        + str(k)

    execute_command(command, verbose)



# GLOG_logtostderr=1  ./rs_index -transcript_fasta=clustered.fa -index_file=clustered_gene.fa.pb -rs_length=60 -num_threads 4
def build_index(cluster_output, index_file, k, numthreads):
    command = "rs_index -transcript_fasta=" \
        + cluster_output \
        + " -index_file=" \
        + index_file \
        + " -rs_length=" \
        + str(k) \
        + " -num_threads=" \
        + str(numthreads)

    execute_command(command, verbose)



# GLOG_logtostderr=1 ./rs_select -index_file=clustered_gene.fa.pb -selected_keys_file=clustered_gene.fa.sk
def select(index_file, selected_keys_file):
    command = "rs_select -index_file=" \
        + index_file \
        + " -selected_keys_file=" \
        + selected_keys_file
    
    execute_command(command, verbose)


# GLOG_logtostderr=1 ./rs_select -index_file=clustered_gene.fa.pb -selected_keys_file=clustered_gene.fa.sk  -rs_length=60
def select_with_k(index_file, selected_keys_file, k):
    command = "rs_select -index_file=" \
        + index_file \
        + " -selected_keys_file=" \
        + selected_keys_file \
        + " -rs_length=" \
        + str(k)
    
    execute_command(command, verbose)


# GLOG_logtostderr=1  ../src/rs_count  -selected_keys_file=clustered_gene.fa.sk -count_file=clustered_gene.fa.cf -read_files1=../test/test.fastq_1 -read_files2=../test/test.fastq_2 -num_threads=1
def count(selected_keys_file, count_file, sample_pair_1, sample_pair_2, numthreads):
    command = "rs_count -selected_keys_file=" \
        + selected_keys_file \
        + " -count_file=" \
        + count_file \
        + " -read_files1=" \
        + sample_pair_1 \
        + " -read_files2=" \
        + sample_pair_2 \
        + " -num_threads=" \
        + str(numthreads)

    execute_command(command, verbose)


def estimate(count_file, estimation_file):
    command = "rs_estimate -count_file=" \
        + count_file

    output = get_command_output(command, verbose)
    text_file = open(estimation_file, "w")
    text_file.write(str(output))
    text_file.close()
    #print("Command Output:\n")
    #print(output)
    #print("----------------------------")


def get_result_dict(result_dir):       
    matrix = np.genfromtxt(
        result_dir + '/estimation.txt',
        names = None,
        dtype=None,
        delimiter="\t")

    transcripts = [x[0] for x in matrix]
    num_reads = [x[2] for x in matrix]

    res = dict()

    for i in range(0, len(transcripts)):
        res[transcripts[i]] = num_reads[i]
    #print(result_dir + ": " + str(len(res)))

    return res


def run(k, transcriptome_reference_file, index_dir, sample_dir, output_dir):
    cleanup_dir(index_dir)
    time1 = time.time()
    numthreads = 4
    cluster_output = index_dir + "/clustered.fa"
    index_file = index_dir + "/clustered_gene.fa.pb"
    selected_keys_file = index_dir + "/clustered_gene.fa.sk"

    cluster(transcriptome_reference_file, cluster_output, k, numthreads)
    build_index(cluster_output, index_file, k, numthreads)
    select_with_k(index_file, selected_keys_file, k)


    res_dict = dict()
    for i in range(1, 11):
        if i < 10:
            sample_pair1 = sample_dir + '/sample_0' + str(i)+ '_1.fasta'
            sample_pair2 = sample_dir + '/sample_0' + str(i) + '_2.fasta'
        else:
            sample_pair1 = sample_dir + '/sample_' + str(i)+ '_1.fasta'
            sample_pair2 = sample_dir + '/sample_' + str(i) + '_2.fasta'
        result_output_dir = output_dir + '/sample' + str(i) + '_result'
        cleanup_dir(result_output_dir)
        count_file = result_output_dir + "/clustered_gene.fa.cf"
        estimation_file = result_output_dir + "/estimation.txt"

        count(selected_keys_file, count_file, sample_pair1, sample_pair2, numthreads)
        estimate(count_file, estimation_file)

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
    
    return res_dict, elapsed_ms, dict()


# run_RNASkim(60, "chr22_small.fa", "test_results/RNASkim/index", "simulated_reads", "test_results/RNASkim/results", 4)



# get_result_dict("/home/ubuntu/cs229/rnaskim")



# cluster("chr22_first4.fa", "RNASkim/src/test_results/clustered.fa", 60, 4)
# cluster("RNASkim/data/homo_sapiens/current/genes.fa", "RNASkim/src/test_results/clustered.fa", 60, 4)
# cluster("RNASkim/src/test_data/gene_fasta_example.fa", "RNASkim/src/test_results/clustered.fa", 60, 4)
# build_index("RNASkim/src/test_results/clustered.fa", "RNASkim/src/test_results/clustered_gene.fa.pb", 60, 4)
# select_with_k("RNASkim/src/test_results/clustered_gene.fa.pb", "RNASkim/src/test_results/clustered_gene.fa.sk", 60)
# count("RNASkim/src/test_results/clustered_gene.fa.sk", "RNASkim/src/test_results/clustered_gene.fa.cf", "RNASkim/src/test_data/fa_reader_test.fasta.1", "RNASkim/src/test_data/fa_reader_test.fasta.2", 1)
# estimate("RNASkim/src/test_results/clustered_gene.fa.cf", "RNASkim/src/test_results/estimation")


