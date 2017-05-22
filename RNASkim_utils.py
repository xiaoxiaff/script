from general_utils import execute_command

# GLOG_logtostderr=1 ./rs_cluster  -gene_fasta=gene.fa -num_threads=4 -output=clustered.fa -rs_length=60
def cluster(transcriptome_reference_file, cluster_output, k, numthreads):
	command = "GLOG_logtostderr=1 rs_cluster -gene_fasta=" \
		+ transcriptome_reference_file \
		+ " -num_threads=" \
		+ str(numthreads) \
		+ " -output=" \
		+ cluster_output \
		+ " -rs_length=" \
		+ str(k)

	execute_command(command, True)



# GLOG_logtostderr=1  ./rs_index -transcript_fasta=clustered.fa -index_file=clustered_gene.fa.pb -rs_length=60 -num_threads 4
def build_index(cluster_output, index_file, k, numthreads):
	command = "GLOG_logtostderr=1 rs_index -transcript_fasta=" \
		+ cluster_output \
		+ " -index_file=" \
		+ index_file \
		+ " -rs_length=" \
		+ str(k) \
		+ " -num_threads=" \
		+ str(numthreads)

	execute_command(command, True)



# GLOG_logtostderr=1 ./rs_select -index_file=clustered_gene.fa.pb -selected_keys_file=clustered_gene.fa.sk
def select(index_file, selected_keys_file):
	command = "GLOG_logtostderr=1 rs_select -index_file=" \
		+ index_file \
		+ " -selected_keys_file=" \
		+ selected_keys_file
	
	execute_command(command, True)

# GLOG_logtostderr=1 ./rs_select -index_file=clustered_gene.fa.pb -selected_keys_file=clustered_gene.fa.sk  -rs_length=60
def select_with_k(index_file, selected_keys_file, k):
	command = "GLOG_logtostderr=1 rs_select -index_file=" \
		+ index_file \
		+ " -selected_keys_file=" \
		+ selected_keys_file \
		+ " -rs_length=" \
		+ str(k)
	
	execute_command(command, True)


# GLOG_logtostderr=1  ../src/rs_count  -selected_keys_file=clustered_gene.fa.sk -count_file=clustered_gene.fa.cf -read_files1=../test/test.fastq_1 -read_files2=../test/test.fastq_2 -num_threads=1
def count(selected_keys_file, count_file, sample_pair_1, sample_pair_2, numthreads):
	command = "GLOG_logtostderr=1 rs_count -selected_keys_file=" \
		+ selected_keys_file \
		+ " -count_file=" \
		+ count_file \
		+ " -read_files1=" \
		+ sample_pair_1 \
		+ " -read_files2=" \
		+ sample_pair_2 \
		+ " -num_threads=" \
		+ str(numthreads)

	execute_command(command, True)


def estimate(count_file, estimation_file):
	command = "rs_estimate -count_file=" \
		+ count_file \
		+ " > " \
		+ estimation_file

	execute_command(command, True)



# cluster("chr22_first4.fa", "RNASkim/src/test_results/clustered.fa", 60, 4)
# cluster("RNASkim/data/homo_sapiens/current/genes.fa", "RNASkim/src/test_results/clustered.fa", 60, 4)
# cluster("RNASkim/src/test_data/gene_fasta_example.fa", "RNASkim/src/test_results/clustered.fa", 60, 4)
# build_index("RNASkim/src/test_results/clustered.fa", "RNASkim/src/test_results/clustered_gene.fa.pb", 60, 4)
# select_with_k("RNASkim/src/test_results/clustered_gene.fa.pb", "RNASkim/src/test_results/clustered_gene.fa.sk", 60)
# count("RNASkim/src/test_results/clustered_gene.fa.sk", "RNASkim/src/test_results/clustered_gene.fa.cf", "RNASkim/src/test_data/fa_reader_test.fasta.1", "RNASkim/src/test_data/fa_reader_test.fasta.2", 1)
# estimate("RNASkim/src/test_results/clustered_gene.fa.cf", "RNASkim/src/test_results/estimation")


