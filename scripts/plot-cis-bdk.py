#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division

import subprocess

from matplotlib.font_manager import FontProperties
from pylab import *
import numpy as np
import scipy as sp
import scikits.bootstrap as bootstrap
import matplotlib.ticker as ticker

# *********************** BOXPLOT ***************************************
def print_ci(data, figxlabel, Nr, Nc, subplot):
    print data_a
    print np.min(data_a), np.max(data_a)
    print np.mean(data_a)
    print np.std(data_a)

    mean_a = np.mean(data)

    ax = fig.add_subplot(Nr, Nc, subplot)

    #plt.setp(ax.get_yticklabels(), visible=False)
    ax.yaxis.set_ticklabels([])
    plt.yticks(size=16)  # , weight='bold')

    plt.setp(ax.get_xticklabels(), visible=False)

    ax.set_ylim([0.00, 3.50])
    ax.set_axisbelow(True)
    ax.grid(True)
    ax.set_xlabel(figxlabel, fontproperties=FontProperties(size=18)) #, weight='bold'))

    #ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f')) #%.1f'

    annotate('%.2f' % mean_a, xy=(1.10, mean_a), xytext=(1.10, mean_a), color='black',
             fontproperties=FontProperties(size=14, weight='bold'))

    if (subplot == 1):
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f')) #%.1f'
        plt.ylabel(r'$\beta \cdot | \delta k,n|$ (in s)', fontproperties=FontProperties(size=24)) #, weight='bold'

    #ax.boxplot(data, notch=False, sym='', whis=[5, 95], meanline=True, showmeans=True, showcaps=True, showfliers=True)
    ax.boxplot(data, notch=False, sym='', whis=[5, 95], meanline=True, showmeans=True)



# ****************************************************************************
# Main block of the program
subprocess.call('clear', shell=True)  # clearing stdio
print '---> Generating graphs ! <---'

# ****************************************************************************
fig = plt.figure(figsize=(9,6), dpi=200)
#fig = plt.figure(figsize=(7, 4), dpi=200)
fig.subplots_adjust(left=0.10, bottom=0.07, right=0.99, top=0.97, wspace=0.15, hspace=0.15)
Nr = 1
Nc = 3
figname = 'plot-ci-bdk'

# ****************************************************************************
# ************** Constant Distribution Delay
# ****************************************************************************
#data_a = np.loadtxt('../stats/data-constant_bdk.txt', dtype=np.float, comments='#', delimiter=',',
#                    converters=None,
#                    skiprows=0, usecols=None, unpack=False)

#figxlabel = 'Constant'
#subplot = 1
#print_ci(data_a, figxlabel, Nr, Nc, subplot)

# ****************************************************************************
# ************** Variable Uniform  Distribution Delay
# ****************************************************************************
data_a = np.loadtxt('../stats/data-uniform_bdk.txt', dtype=np.float, comments='#', delimiter=',',
                    converters=None,
                    skiprows=0, usecols=None, unpack=False)

figxlabel = 'Uniform'
subplot = 1
print_ci(data_a, figxlabel, Nr, Nc, subplot)

# ****************************************************************************
# ************** Gaussian Distribution Delay
# ****************************************************************************
data_a = np.loadtxt('../stats/data-gauss_bdk.txt', dtype=np.float, comments='#', delimiter=',', converters=None,
                    skiprows=0, usecols=None, unpack=False)

figxlabel = 'Gaussian'
subplot = 2
print_ci(data_a, figxlabel, Nr, Nc, subplot)

# ****************************************************************************
# ************** Exponential Distribution Delay
# ****************************************************************************
data_a = np.loadtxt('../stats/data-exponential_bdk.txt', dtype=np.float, comments='#', delimiter=',',
                    converters=None,
                    skiprows=0, usecols=None, unpack=False)

figxlabel = 'Exponential'
subplot = 3
print_ci(data_a, figxlabel, Nr, Nc, subplot)

# ***************************************************************************
#plt.tight_layout()
#plt.savefig(figname + '.jpg')
plt.savefig(figname + '.eps') # was .pdf
ion()
plt.draw()
plt.show()
# ***************************************************************************
