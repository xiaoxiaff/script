import sys
import os
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from general_utils import execute_command
from general_utils import remove_file_if_exists
from general_utils import get_average_accuracy
from general_utils import cleanup_dir
import datetime
from timeit import timeit


def plot_accuracy_for_tool(tool_name, plot_type, labels, k_range, accuracy_matrix):
    n_groups = len(k_range)
    index = np.arange(n_groups)

    tick_range = []
    tick_index_range = []
    if(len(k_range)>20):
        tick_range = k_range[::3]
        tick_index_range = index[::3]
    else:
        tick_range = k_range
        tick_index_range = index
    colors = ['r','g','b','y','m','c','k', 'pink']
    plt.figure()
    axes = plt.gca()
    # axes.set_xlim([xmin,xmax])
    # axes.set_ylim([0,1])
    plt.xlabel('k value')
    plt.ylabel('Accuracy')
    plot_title = "{0} {1} accuracy plot".format(tool_name, plot_type)
    plt.title(plot_title)
    for i in range(0,len(accuracy_matrix)):
        current_array = accuracy_matrix[i]
        besk_k = k_range[np.argmax(current_array)]
        best_accuracy = max(current_array)
        label = "{0}={1},bestK={2},accuracy={3:.4f}".format(plot_type,str(labels[i]),str(besk_k),best_accuracy)
        plt.plot(index, current_array, color = colors[i], ls='-', marker='o', label=label)
    
    plt.xticks(tick_index_range, tick_range)
    plt.legend()
    plt.subplots_adjust(left=0.15)
    plt.savefig(tool_name + "_" + plot_type + "_accuracy_plot")
    plt.close()

def plot_duck_for_tool(tool_name, plot_type, labels, k_range, duck_matrix):
    n_groups = len(k_range)
    index = np.arange(n_groups)

    tick_range = []
    tick_index_range = []
    if(len(k_range)>20):
        tick_range = k_range[::3]
        tick_index_range = index[::3]
    else:
        tick_range = k_range
        tick_index_range = index
    colors = ['r','g','b','y','m','c','k', 'pink']
    if tool_name == "salmon":
        for seq in range(2):
            plt.figure()
            axes = plt.gca()
            # axes.set_xlim([xmin,xmax])
            # axes.set_ylim([0,1])
            plt.xlabel('k value')
            if seq == 0:
                plt.ylabel('hash table size')
                plot_title = "{0} {1} hash table size".format(tool_name, plot_type)
            elif seq == 1:
                plt.ylabel('invalid kmer')
                plot_title = "{0} {1} invalid kmer counts".format(tool_name, plot_type)
            
            plt.title(plot_title)
            for i in range(0,len(labels)):
                current_array = duck_matrix[i]
                result_array = []
                for item in current_array:
                    if seq == 0:
                        result_array.append(item['khashSize'])
                    if seq == 1:
                        result_array.append(item['invalidTime'])
                label = "{0}={1}".format(plot_type,str(labels[i]))
                plt.plot(index, result_array, color = colors[i], ls='-', marker='o', label=label)
            
            plt.xticks(tick_index_range, tick_range)
            plt.legend()
            plt.subplots_adjust(left=0.15)
            if seq == 0:
                plt.savefig(tool_name + "_" + plot_type + "_table_size_plot")
            else:
                plt.savefig(tool_name + "_" + plot_type + "_invalid_count_plot")
            plt.close()
    if tool_name == "kallisto":
        for seq in range(2):
            plt.figure()
            axes = plt.gca()
            # axes.set_xlim([xmin,xmax])
            # axes.set_ylim([0,1])
            plt.xlabel('k value')
            if seq == 0:
                plt.ylabel('kmers count')
                plot_title = "{0} {1} kmers count".format(tool_name, plot_type)
            elif seq == 1:
                plt.ylabel('contigs count')
                plot_title = "{0} {1} contigs count".format(tool_name, plot_type)
            
            plt.title(plot_title)
            for i in range(0,len(labels)):
                current_array = duck_matrix[i]
                result_array = []
                for item in current_array:
                    if seq == 0:
                        result_array.append(item['kmers'])
                    if seq == 1:
                        result_array.append(item['contigs'])
                label = "{0}={1}".format(plot_type,str(labels[i]))
                plt.plot(index, result_array, color = colors[i], ls='-', marker='o', label=label)
            
            plt.xticks(tick_index_range, tick_range)
            plt.legend()
            plt.subplots_adjust(left=0.15)
            if seq == 0:
                plt.savefig(tool_name + "_" + plot_type + "_kmers_plot")
            else:
                plt.savefig(tool_name + "_" + plot_type + "_contigs_plot")
            plt.close()

    if tool_name == "sailfish":
        attrList = ['transcriptHashSize','numOfTranscript','numOfGene', 'numOfKmers', 'numOfEquivClass', 'transcriptForKmerTableSize']
        for seq in range(len(attrList)):
            plt.figure()
            axes = plt.gca()
            # axes.set_xlim([xmin,xmax])
            # axes.set_ylim([0,1])
            plt.xlabel('k value')

            plot_title = "{0}, {1}".format(tool_name, attrList[seq])
            
            plt.title(plot_title)
            current_array = duck_matrix[0]
            result_array = []
            for item in current_array:
                result_array.append(item[attrList[seq]])
            label = "{0}".format(tool_name)
            plt.plot(index, result_array, color = colors[0], ls='-', marker='o', label=label)
           
            plt.xticks(tick_index_range, tick_range)
            plt.legend()
            plt.subplots_adjust(left=0.15)
            plt.savefig(tool_name + "_" + attrList[seq] + "_plot")
            plt.close()


def plot_runtime_for_tool(tool_name, plot_type, labels, k_range, runtime_matrix):
    n_groups = len(k_range)
    index = np.arange(n_groups)

    tick_range = []
    tick_index_range = []
    if(len(k_range)>20):
        tick_range = k_range[::3]
        tick_index_range = index[::3]
    else:
        tick_range = k_range
        tick_index_range = index

    colors = ['r','g','b','y','m','c','k', 'pink']
    plt.figure()
    axes = plt.gca()
    # axes.set_xlim([xmin,xmax])
    # axes.set_ylim([0,1])
    plt.xlabel('k value')
    plt.ylabel('Runtime (ms)')
    plot_title = "{0} {1} runtime plot".format(tool_name, plot_type)
    plt.title(plot_title)

    for i in range(0,len(runtime_matrix)):
        current_array = runtime_matrix[i]
        label = "{0}={1}".format(plot_type,str(labels[i]))
        plt.plot(index, current_array, color = colors[i], ls='-', marker='o', label=label)
    
    plt.xticks(tick_index_range, tick_range)
    plt.legend()
    plt.subplots_adjust(left=0.15)
    plt.savefig(tool_name + "_" + plot_type + "_runtime_plot")
    plt.close()


def save_result_matrix_as_csv(tool_name, result_type, loop_type, k_range, loop_range, result_matrix):
    with open(tool_name + "_" + result_type + "_" + loop_type + ".csv", "w") as outfile:
        first_line = ""
        for k in k_range:
            first_line = first_line + "," + "k=" + str(k)
        outfile.write(first_line + "\n")


        for i in range(0,len(loop_range)):
            line = loop_type + "=" + str(loop_range[i])
            for element in result_matrix[i]:
                line = line + "," + str(element)
            outfile.write(line + "\n")
