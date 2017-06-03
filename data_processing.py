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
