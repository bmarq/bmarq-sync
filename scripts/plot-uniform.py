#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import unicode_literals
import matplotlib as mpl

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True

from matplotlib.font_manager import FontProperties
import subprocess
import numpy as np
from matplotlib import pyplot
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.ticker as ticker
import matplotlib.lines as lines
from matplotlib.lines import Line2D


# *********************** PLOTS ***********************************
def print_ci(data_a, data_b, figxlabel, Nr, Nc, subplot):

    print 'Subplot: %d' % (subplot)

    TON = 60
    delta = 0.80 * TON
    #delta = 0.80 * np.mean(data_a)
    data_std = np.std(data_a)
    mean_a = np.mean(data_a)

    fig.suptitle(r'$\alpha=0.875; \beta=100;$', fontsize=16, fontweight='bold')

    ax = fig.add_subplot(Nr, Nc, subplot)

    plt.yticks(size=13)  # , weight='bold')

    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    ax.set_axisbelow(True)
    ax.grid(True)

    if (subplot == 1):
        print data_a, np.min(data_a), np.max(data_a), np.mean(data_a), np.std(data_a)
        print data_b, np.min(data_b), np.max(data_b), np.mean(data_b), np.std(data_b)
        print

        ax.set_ylim([0.0, 0.5])
        ax.set_xlim([0.0, 3.0])

        weights = np.ones_like(data_b) / len(data_b)
        n, bins, patches = ax.hist(data_b, bins=5, weights=weights, facecolor="limegreen", alpha=0.75, label='Node 2')

        weights = np.ones_like(data_a) / len(data_a)
        n, bins, patches = ax.hist(data_a, bins=5, weights=weights, facecolor="r", alpha=0.75, label='Node 1')

        ax.legend(loc='upper right', fontsize='13')
        ax.set_ylabel('Occurrence', fontsize='14', weight='bold')
        ax.set_xlabel(r'a) Delay (s)', fontsize='14')

    elif (subplot == 2):

        print 'data\t min\t max\t mean\t std\t delta'
        print data_a, np.min(data_a), np.max(data_a), np.mean(data_a), data_std, delta
        print

        ax.set_ylim([0.0, 0.6])
        ax.set_xlim([45.0, 60.0])

        weights = np.ones_like(data_a) / len(data_a)
        n, bins, patches = ax.hist(data_a, bins=5, weights=weights, facecolor="blue", alpha=0.75)

        # ax.legend(loc='lower right', fontsize='14')
        line1 = [(delta, 0.0), (delta, 0.5)]
        (line1_xs, line1_ys) = zip(*line1)
        ax.add_line(Line2D(line1_xs, line1_ys, linewidth=2, color='red'))
        plt.text(delta - 1, -0.05, r'$\nabla$', color='red', fontsize=18)

        plt.text(-9.0, 0.39,
                 r'$E[T_{Sensors_{ON}}]=$' + str(round(mean_a, 1)) + '\n' + r'$Std[T_{Sensors_{ON}}]=$' + str(
                     round(data_std, 1)) + '\n' + r'$\Delta=$' + str(round(delta, 1)), fontsize=14)

        ax.set_ylabel('Occurrence', fontsize='14', weight='bold')
        ax.set_xlabel(r'b) $T_{Sensors_{ON}}$ (s)', fontsize='18')  # , weight='bold')

    elif (subplot == 3):

        print 'data\t min\t max\t mean\t std\t delta'
        print data_a, np.min(data_a), np.max(data_a), np.mean(data_a), data_std, delta
        print

        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        ax.set_ylabel(r'$T_{ON}$', fontsize='16', weight='bold')
        ax.boxplot(data_a, notch=False, sym='', whis=[5, 95], meanline=True, showmeans=True)

        ax.set_xlabel(r'c) $T_{Sensors_{ON}}$ (\% of $T_{ON}$)', fontsize='18')

        ax.annotate('%.2f' % mean_a, xy=(1.10, mean_a), xytext=(1.10, mean_a), color='black',
                    fontproperties=FontProperties(size=14, weight='bold'))

    elif (subplot == 4):

        print 'data\t min\t max\t mean\t std\t delta'
        print data_a, np.min(data_a), np.max(data_a), np.mean(data_a), data_std, delta
        print

        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
        ax.set_ylabel('Time', fontsize='14', weight='bold')
        ax.boxplot(data_a, notch=False, sym='', whis=[5, 95], meanline=True, showmeans=True)

        ax.set_xlabel(r'd) $\beta \cdot |\delta_{k,n}|$ (s)', fontsize='18')

        ax.annotate('%.3f' % mean_a, xy=(1.10, mean_a), xytext=(1.10, mean_a), color='black',
                    fontproperties=FontProperties(size=14, weight='bold'))


'''---8<------8<------8<------8<------8<------8<------8<------8<------8<--- '''
# Main block of the program
subprocess.call('clear', shell=True)  # clearing stdio
'''---8<------8<------8<------8<------8<------8<------8<------8<------8<--- '''

print 'Generating graphs ! Wait !...'

# ****************************************************************************
fig = plt.figure(figsize=(17, 6.5), dpi=200)
fig.subplots_adjust(left=0.07, bottom=0.15, right=0.99, top=0.90, wspace=0.55, hspace=0.5)
Nr = 1
Nc = 4

figxlabel = 'Uniform'

# ************************* Ploting Constant Delay Graphs
# ************************* Delay Histogram
data_a = np.loadtxt('../results/uniform/data-delay-node1-tcycle_uniform-delay_uniform-sim_1.txt', dtype=np.float,
                    comments='#', delimiter=',',
                    converters=None,
                    skiprows=0, usecols=None, unpack=False)
data_b = np.loadtxt('../results/uniform/data-delay-node2-tcycle_uniform-delay_uniform-sim_1.txt', dtype=np.float,
                    comments='#', delimiter=',',
                    converters=None,
                    skiprows=0, usecols=None, unpack=False)
subplot = 1
print_ci(data_a, data_b, figxlabel, Nr, Nc, subplot)

# ************************* T_Sensors_ON Histogram
data_a = np.loadtxt(
    '../results/uniform/data-tsensors_on-tcycle_uniform-delay_uniform-gamma_0.8-alpha_0.125-beta_1.0-sim_1.txt',
    dtype=np.float, comments='#', delimiter=',',
                    converters=None,
                    skiprows=0, usecols=None, unpack=False)
data_b = []

subplot = 2
print_ci(data_a, data_b, figxlabel, Nr, Nc, subplot)

# ************************* T_Sensors_ON in % of TON
data_a = np.loadtxt(
    '../results/uniform/data-tsensors_on-percent-tcycle_uniform-delay_uniform-gamma_0.8-alpha_0.125-beta_1.0-sim_1.txt',
    dtype=np.float, comments='#', delimiter=',',
                    converters=None,
                    skiprows=0, usecols=None, unpack=False)
data_b = []

subplot = 3
print_ci(data_a, data_b, figxlabel, Nr, Nc, subplot)

# ************************* Sleeping offset - b.|dk|
data_a = np.loadtxt(
    '../results/uniform/data-sleepOffset-node1-tcycle_uniform-delay_uniform-gamma_0.8-alpha_0.125-beta_1.0-sim_1.txt',
    dtype=np.float, comments='#', delimiter=',',
                    converters=None,
                    skiprows=0, usecols=None, unpack=False)
data_b = []

subplot = 4
print_ci(data_a, data_b, figxlabel, Nr, Nc, subplot)

# ****************************************************************************
plt.savefig('../graphs/plot-uniform-gama_08-alpha-0125_beta-1.png')

print 'All Done!'