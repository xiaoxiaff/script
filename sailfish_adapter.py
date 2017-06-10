from general_utils import execute_command
from general_utils import cleanup_dir
import numpy as np
import time
import os
import sys

verbose = False
def set_verbose(v):
    global verbose
    verbose = v


def build_index(k, transcriptome_reference_file, index_output_path):
    command = "sailfish index -t " \
        + transcriptome_reference_file \
        + " -o " \
        + index_output_path \
        + " -k " \
        + str(k) \
        + " --force"

    print(command)

    execute_command(command,verbose)

def evaluate_result(index_dir):
    command = "analyze_kmer_tradeoff -c " \
        + index_dir+"/jf.counts_0" \
        + " -i " \
        + index_dir+"/transcriptome" \
        + " -b " \
        + index_dir+"/kmerEquivClasses" \
        + " -o " \
        + index_dir+"/duck" \
        + " -l " \
        + index_dir+"/transcriptome"

    execute_command(command,verbose)    

def quant_with_k(sample_pair1, sample_pair2, index_dir, output_dir):
    command = "sailfish quant -i " \
        + index_dir \
        + " -l " \
        + "'T=PE:O=><:S=SA'" \
        + " -1 " \
        + sample_pair1 \
        + " -2 " \
        + sample_pair2 \
        + " -p 4 -o " \
        + output_dir

    #execute_command(command,verbose)
    import subprocess
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()
    print (process.returncode)
    print("execute command finished")

    print(command)
    print(os.path.exists(output_dir))
    print(os.path.exists(output_dir+'/quant.sf'))

def get_result_dict(result_dir):
    print(result_dir)
    print(os.path.exists(result_dir))
    print(os.path.exists(result_dir+'/quant.sf'))
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


def run(k, transcriptome_reference_file, index_dir, sample_dir, output_dir):
    cleanup_dir(index_dir)
    time1 = time.time()
    duck_dict = dict()
    build_index(k, transcriptome_reference_file, index_dir)
    try:
        f = open('duck2.txt', "r+")
        outputStr = f.readline()
        print(outputStr)
        print(outputStr[outputStr.index("transcriptHash size")+19:])
        duck_dict['transcriptHashSize'] = outputStr[outputStr.index("transcriptHash size")+19:]
        outputStr = f.readline()
        outputStr = f.readline()
        print(outputStr)
        duck_dict['numOfTranscript'] = outputStr[outputStr.index("numTranscripts")+14:]
        outputStr = f.readline()
        print(outputStr)
        duck_dict['numOfGene'] = outputStr[outputStr.index("numGene")+7:]
        outputStr = f.readline()
        print(outputStr)
        duck_dict['numOfKmers'] = outputStr[outputStr.index("num of kmers")+12:]
        outputStr = f.readline()
        print(outputStr)
        duck_dict['numOfEquivClass'] = outputStr[outputStr.index("num of equivClass")+17:]
        outputStr = f.readline()
        print(outputStr)
        duck_dict['transcriptForKmerTableSize'] = outputStr[outputStr.index("tforktable size")+15:]
    except:
        print("Unexpected error:", sys.exc_info()[0])

    print(duck_dict)
    res_dict = dict()
          

    time2 = time.time()
    elapsed_ms = (time2-time1)*1000.0
    
    return res_dict, elapsed_ms, duck_dict


# run_sailfish(31, "chr22_small.fa", "test_results/sailfish/index", "simulated_reads", "test_results/sailfish/results")

