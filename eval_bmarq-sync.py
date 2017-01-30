#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Copyright (c) 2017, Bruno Marques, INESC TEC, IPV/ESTGV, https://github.com/bmarq/bmarq-sync
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.
'''

from __future__ import division
import os
import subprocess
import argparse
import time

from pylab import *
import numpy as np
import scipy as sp
from pylab import *

__version__ = "eval_bmarq.py v0.2, (C) 2017, Bruno Marques, INESC TEC, IPV/ESTGV"

global TON, TOFF, TCYCLE, TMAX, tMaxCycle, tMinCycle
global rnd_tDelayDist, tDelayDist, rnd_tCycleDist, tCycleDist
global nSim
global pDiscard
global n1, n2, n3, n4
global delay_1, delay_2, delay_3, delay_4
global delay1, delay2, delay3, delay4
global expected1_tn, expected2_tn, expected3_tn, expected4_tn
global sampled1_tn, sampled2_tn, sampled3_tn, sampled4_tn
global sampled1_tn_1, sampled2_tn_1, sampled3_tn_1, sampled4_tn_1
global delta1_k, delta2_k, delta3_k, delta4_k
global delta1_k_1, delta2_k_1, delta3_k_1
global tsleep1_k, tsleep2_k, tsleep3_k, tsleep4_k
global r1on, r2on, r3on, r4on
global r1off, r3off, r3off, r4off
global alpha, beta, gamma, sigma
global i_counter, hit
global lower, upper
global success, temp_success, delta_success, success_counter
global tsensors_on, tsensors_on_percent

# ***************************************************************************

parser = argparse.ArgumentParser(description='Evaluation of the bmarq-sync sycnhronization mechanism\n')
parser.add_argument('--version', action='version', version=__version__, help='Version number of eval_bmarq.py')
parser.add_argument('--nsim', type=int, default=1, help='Number of simulations to perform (default: 1)')
parser.add_argument('--maxcicles', type=int, default=1000,
                    help='Maximum number of cycles per simulation (default: 1000)')
parser.add_argument('--alpha', type=float, default=0.125, choices=[0.125, 0.50, 0.875],
                    help='Value for alpha parameter default: 0.125)')
parser.add_argument('--beta', type=float, default=10, choices=[1, 10, 50, 100],
                    help='Value for beta parameter (default: 10)')
parser.add_argument('--gamma', type=float, default=0.80, choices={0.5, 0.6, 0.7, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0},
                    help='%% of TON for TSensorsOn success (default: 0.80)')
parser.add_argument('--sigma', type=float, default=0.20, choices={0.0, 0.05, 0.10, 0.15, 0.20, 0.25},
                    help='Value of standard deviation for delays (default: 0.20)')
parser.add_argument('--ton', type=float, default=60, help='Value for TON (default: 60)')
parser.add_argument('--tmaxcycle', type=float, default=900, help='Maximum value for TCycle (default: 900)')
parser.add_argument('--tmincycle', type=float, default=120, help='Minimum value for TCycle (default: 120)')
parser.add_argument('--pdiscard', type=float, default=0.10,
                    help='Initial %% of cycles to discard for estability purposes [default: 0.10 (10%%)]')
parser.add_argument('--tdelaydist', type=str, default='uniform',
                    choices={'constant', 'uniform', 'normal', 'exponential', 'chisquare', 'poisson'},
                    help='Type of random distribution for delays (default: uniform)')
parser.add_argument('--tcycledist', type=str, default='uniform',
                    choices={'constant', 'uniform', 'normal', 'exponential', 'chisquare', 'poisson'},
                    help='Type of random distribution for TCycle (default: uniform). If the distribution is constant, the default value equals tMinCycle')

args = parser.parse_args()
nSim = args.nsim
TMAX = args.maxcicles
alpha = args.alpha  # test with 0.125; 0.50; 0.875
beta = args.beta  # test with 1, 10, 50, 100
gamma = args.gamma  # % of TON for TsensorsOn success
sigma = args.sigma  # value for delays standard deviation
TON = args.ton
tMaxCycle = args.tmaxcycle
tMinCycle = args.tmincycle
pDiscard = args.pdiscard
tDelayDist = args.tdelaydist
tCycleDist = args.tcycledist

'''---8<------8<------8<------8<------8<------8<------8<------8<------8<--- '''
# Main block of the program

subprocess.call('clear', shell=True)  # clearing stdio
start_t = time.clock()

delay_1 = 0.5
delay_2 = 1
delay_3 = 2
delay_4 = 5

expected1_tn = expected2_tn = expected3_tn = expected4_tn = 0.000  # t'k for nodes 1, 2, 3 and 4 (expected reception time for packet k)
sampled1_tn = sampled2_tn = sampled3_tn = sampled4_tn = 0.000  # tk -1 for nodes 1, 2, 3 and 4 (real reception time for packet k-1 [actual (k)])
sampled1_tn_1 = sampled2_tn_1 = sampled3_tn_1 = sampled4_tn_1 = 0.000  # tk -1 for nodes 1, 2, 3 and 4 (real reception time for packet k-1 [anterior (k-1)])

delta1_k = delta2_k = delta3_k = delta4_k = 0.000
delta1_k_1 = delta2_k_1 = delta3_k_1 = delta4_k_1 = 0.000

tsleep1_k = tsleep2_k = tsleep3_k = tsleep4_k = 0.000

r1on = r2on = r3on = r4on = 0.000
r1off = r2off = r3off = r3off = TON

n1 = 1  # Node 1
n2 = 2  # Node 2
n3 = 3  # Node 3
n4 = 4  # Node 4

i_counter = 0
hit = 0
success = success_counter = temp_success = 0
delta_success = gamma * TON  # gamma % of TON

print 'Started processing at: %f ...\n' % (start_t)

file_eval = open('../results/evaluation-readme.txt', 'w')
f = (
      'number of simulations: %d\ntcycle dist: %s\ndelay dist: %s\nalpha = %.3f\nbeta = %.3f\ngamma = %.3f\nsigma = %.3f\nTON = %.3f\nminimum time for success = %.3f\n') % (
      nSim, tCycleDist, tDelayDist, alpha, beta, gamma, sigma, TON, delta_success)
file_eval.write(f)
file_eval.close()
print f

'''---8<------8<------8<------8<------8<------8<------8<------8<------8<--- '''
for j in range(1, int(nSim) + 1):
  print 'running cycle %d:' % (j)
  success_counter += 1

  # generate random delays between tMinCycle and tMaxCycle with 'tDelayDist' tDelayDistution

  # use the folowing two lines if you want to reproduce the experiment and get the same results
  np.random.seed(j)
  np.random.RandomState(j)

  if (tCycleDist == 'uniform'):
    rnd_tCycleDist = 'np.random.' + tCycleDist + '(tMinCycle, tMaxCycle)'

  if (tCycleDist == 'normal'):
    rnd_tCycleDist = 'np.random.' + tCycleDist + '(tMinCycle, sigma * tMaxCycle)'

  if (tCycleDist == 'exponential'):
    rnd_tCycleDist = 'np.random.' + tCycleDist + '(tMinCycle)'

  if (tCycleDist == 'chisquare'):
    rnd_tCycleDist = 'np.random.' + tCycleDist + '(tMinCycle)'

  if (tCycleDist == 'poisson'):
    rnd_tCycleDist = 'np.random.' + tCycleDist + '(tMinCycle)'

  if (tCycleDist == 'constant'):
    TCYCLE = tMinCycle
  else:
    TCYCLE = eval(rnd_tCycleDist)

  TOFF = TCYCLE - TON

  file_n = open('../results/data-nCycles.txt', 'w')
  file_all = open(
    '../results/results-tcycle-dist-' + tCycleDist + '-delay_' + str(tDelayDist) + '-sim_' + str(j) + '.txt', 'w')
  file_delay1 = open(
    '../results/data-delay-node1-tcycle_' + tCycleDist + '-delay_' + str(tDelayDist) + '-sim_' + str(j) + '.txt', 'w')
  file_delay2 = open(
    '../results/data-delay-node2-tcycle_' + tCycleDist + '-delay_' + str(tDelayDist) + '-sim_' + str(j) + '.txt', 'w')
  file_delay3 = open(
    '../results/data-delay-node3-tcycle_' + tCycleDist + '-delay_' + str(tDelayDist) + '-sim_' + str(j) + '.txt', 'w')
  file_delay4 = open(
    '../results/data-delay-node4-tcycle_' + tCycleDist + '-delay_' + str(tDelayDist) + '-sim_' + str(j) + '.txt', 'w')
  file_tSensorsOn = open(
    '../results/data-tsensors_on-tcycle_' + tCycleDist + '-delay_' + str(tDelayDist) + '-sim_' + str(j) + '.txt', 'w')
  file_tSensorsOn_percent = open(
    '../results/data-tsensors_on-percent-tcycle_' + tCycleDist + '-delay_' + str(tDelayDist) + '-sim_' + str(
      j) + '.txt', 'w')
  file_bdk = open('../results/data-bdk-tcycle_' + tCycleDist + '-delay_' + str(tDelayDist) + '-sim_' + str(j) + '.txt',
                  'w')
  file_tSensorsOn_success = open(
    '../results/data-tsensors_on-success-tcycle_' + '-delay_' + str(tDelayDist) + '-sim_' + str(j) + '.txt', 'w')
  file_tcycle = open('../results/data-tcycle_' + tCycleDist + '-delay_' + str(tDelayDist) + '-sim_' + str(j) + '.txt',
                     'w')

  a = 'cycle\tnode\tt\'k\ttk\tdelay\tt\'k-tk\tdk\tr_on\tr_off\ttsleep\tb.|dk|\ttcycle\tton\ttoff\t\tsensors_on(s)\ttsensors_on(%)\n'
  a = a.expandtabs(8)
  file_all.write(a)

  for n in range(0, int(TMAX) + 1):
    # print n

    # initialization for node n1
    # use the following two lines if you want to reproduce the experiment and get the same results
    np.random.seed(j + n)
    np.random.RandomState(j + n)

    if (tDelayDist == 'uniform'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '((1 - sigma) * delay_1, (1 + sigma) * delay_1)'

    if (tDelayDist == 'normal'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_1, sigma * delay_1)'

    if (tDelayDist == 'exponential'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_1)'

    if (tDelayDist == 'chisquare'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_1)'

    if (tDelayDist == 'poisson'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_1)'

    if (tDelayDist == 'constant'):
      delay1 = delay_1
    else:
      delay1 = eval(rnd_tDelayDist)

    # print 'delay for node 1 in cycle %d: %f3' %(n, delay1)

    sampled1_tn_1 = sampled1_tn
    sampled1_tn = n * TCYCLE + delay1
    expected1_tn = sampled1_tn_1 + TON + TOFF
    delta1_k_1 = delta1_k

    # initialization for node n2
    # use the following two lines if you want to reproduce the experiment and get the same results
    np.random.seed(j + n + 1)
    np.random.RandomState(j + n + 1)

    if (tDelayDist == 'constant'):
      rnd_tDelayDist = delay_2

    if (tDelayDist == 'uniform'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '((1 - sigma) * delay_2, (1 + sigma) * delay_2)'

    if (tDelayDist == 'normal'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_2, sigma * delay_2)'

    if (tDelayDist == 'exponential'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_2)'

    if (tDelayDist == 'chisquare'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_2)'

    if (tDelayDist == 'poisson'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_2)'

    if (tDelayDist == 'constant'):
      delay2 = delay_2
    else:
      delay2 = eval(rnd_tDelayDist)

    # print 'delay for node 2 in cycle %d: %f3' %(n, delay2)

    sampled2_tn_1 = sampled2_tn
    sampled2_tn = n * TCYCLE + delay2
    expected2_tn = sampled2_tn_1 + TON + TOFF
    delta2_k_1 = delta2_k

    # initialization for node n3
    # use the following two lines if you want to reproduce the experiment and get the same results
    np.random.seed(j + n + 2)
    np.random.RandomState(j + n + 2)

    if (tDelayDist == 'constant'):
      rnd_tDelayDist = delay_3

    if (tDelayDist == 'uniform'):
      rnd_tDelayDist == 'np.random.' + tDelayDist + '((1 - sigma) * delay_3, (1 + sigma) * delay_3)'

    if (tDelayDist == 'normal'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_3, sigma * delay_3)'

    if (tDelayDist == 'exponential'):
      rnd_tDelayDist == 'np.random.' + tDelayDist + '(delay_3)'

    if (tDelayDist == 'chisquare'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_3)'

    if (tDelayDist == 'poisson'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_3)'

    if (tDelayDist == 'constant'):
      delay3 = delay_3
    else:
      delay3 = eval(rnd_tDelayDist)

    # print 'delay for node 3 in cycle %d: %f3' %(n, delay3)

    sampled3_tn_1 = sampled3_tn
    sampled3_tn = n * TCYCLE + delay3
    expected3_tn = sampled3_tn_1 + TON + TOFF
    delta3_k_1 = delta3_k

    # initialization for node n4
    # use the following two lines if you want to reproduce the experiment and get the same results
    np.random.seed(j + n + 3)
    np.random.RandomState(j + n + 3)

    if (tDelayDist == 'constant'):
      rnd_tDelayDist = delay_4

    if (tDelayDist == 'uniform'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '((1 - sigma) * delay_4, (1 + sigma) * delay_4)'

    if (tDelayDist == 'normal'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_4, sigma * delay_4)'

    if (tDelayDist == 'exponential'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_4)'

    if (tDelayDist == 'chisquare'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_4)'

    if (tDelayDist == 'poisson'):
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_4)'

    if (tDelayDist == 'constant'):
      delay4 = delay_4
    else:
      delay4 = eval(rnd_tDelayDist)

    # print 'delay for node 4 in cycle %d: %f3' %(n, delay4)

    sampled4_tn_1 = sampled4_tn
    sampled4_tn = n * TCYCLE + delay4
    expected4_tn = sampled4_tn_1 + TON + TOFF
    delta4_k_1 = delta4_k

    if (n == 0):
      delta1_k = delta2_k = delta3_k = delta4_k = 0
      r1on = r2on = r3on = r4on = 0
    else:
      delta1_k = (1.0 - alpha) * delta1_k_1 + alpha * (expected1_tn - sampled1_tn)
      r1on = expected1_tn - beta * abs(delta1_k)

      delta2_k = (1.0 - alpha) * delta2_k_1 + alpha * (expected2_tn - sampled2_tn)
      r2on = expected2_tn - beta * abs(delta2_k)

      delta3_k = (1.0 - alpha) * delta3_k_1 + alpha * (expected3_tn - sampled3_tn)
      r3on = expected3_tn - beta * abs(delta3_k)

      delta4_k = (1.0 - alpha) * delta4_k_1 + alpha * (expected4_tn - sampled4_tn)
      r4on = expected4_tn - beta * abs(delta4_k)

    tsleep1_k = TOFF - beta * abs(delta1_k)
    tsleep2_k = TOFF - beta * abs(delta2_k)
    tsleep3_k = TOFF - beta * abs(delta3_k)
    tsleep4_k = TOFF - beta * abs(delta4_k)

    r1off = r1on + TON
    r2off = r2on + TON
    r3off = r3on + TON
    r4off = r4on + TON

    lower = (r1on, r2on, r3on, r4on)
    upper = (r1off, r2off, r3off, r4off)
    tsensors_on = (min(upper) - max(lower))
    tsensors_on_percent = ((min(upper) - max(lower)) / TON) * 100

    if (n >= int(TMAX * pDiscard)):
      # if (n > 0):
      i_counter += 1
      # Node 1
      a = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\n'.format(str(n),
                                                                                              str(n1),
                                                                                              str(
                                                                                                round(expected1_tn, 3)),
                                                                                              str(
                                                                                                round(sampled1_tn, 3)),
                                                                                              str(round(delay1, 3)),
                                                                                              str(
                                                                                                round(
                                                                                                  expected1_tn - sampled1_tn,
                                                                                                  3)),
                                                                                              str(
                                                                                                round(
                                                                                                  delta1_k,
                                                                                                  3)),
                                                                                              str(round(r1on,
                                                                                                        3)),
                                                                                              str(round(r1off,
                                                                                                        3)),
                                                                                              str(round(
                                                                                                tsleep1_k,
                                                                                                3)),
                                                                                              str(round(beta * abs(
                                                                                                delta1_k),
                                                                                                        3)),
                                                                                              str(round(TCYCLE, 3)),
                                                                                              str(round(TON, 3)),
                                                                                              str(round(TOFF, 3)),
                                                                                              str(
                                                                                                round(
                                                                                                  tsensors_on,
                                                                                                  3)),
                                                                                              str(round(
                                                                                                tsensors_on_percent,
                                                                                                3)))

      a = a.expandtabs(8)
      file_all.write(a)

      b1 = '{0}\n'.format(str(round(delay1, 3)))
      b1 = b1.expandtabs(8)
      file_delay1.write(b1)

      # Node 2
      a = '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\n'.format(str(n2),
                                                                                          str(
                                                                                            round(expected2_tn, 3)),
                                                                                          str(
                                                                                            round(sampled2_tn, 3)),
                                                                                          str(round(delay2, 3)),
                                                                                          str(
                                                                                            round(
                                                                                              expected2_tn - sampled2_tn,
                                                                                              3)),
                                                                                          str(
                                                                                            round(
                                                                                              delta2_k,
                                                                                              3)),
                                                                                          str(round(r2on,
                                                                                                    3)),
                                                                                          str(round(r2off,
                                                                                                    3)),
                                                                                          str(round(
                                                                                            tsleep2_k,
                                                                                            3)),
                                                                                          str(round(beta * abs(
                                                                                            delta2_k), 3)),
                                                                                          str(round(TCYCLE, 3)),
                                                                                          str(round(TON, 3)),
                                                                                          str(round(TOFF, 3)),
                                                                                          str(
                                                                                            round(
                                                                                              tsensors_on,
                                                                                              3)),
                                                                                          str(round(
                                                                                            tsensors_on_percent,
                                                                                            3)))
      a = a.expandtabs(8)
      file_all.write(a)

      b2 = '{0}\n'.format(str(round(delay2, 3)))
      b2 = b2.expandtabs(8)
      file_delay2.write(b2)

      # Node 3
      a = '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\n'.format(str(n3),
                                                                                          str(
                                                                                            round(expected3_tn, 3)),
                                                                                          str(
                                                                                            round(sampled3_tn, 3)),
                                                                                          str(round(delay3, 3)),
                                                                                          str(
                                                                                            round(
                                                                                              expected3_tn - sampled3_tn,
                                                                                              3)),
                                                                                          str(
                                                                                            round(
                                                                                              delta3_k,
                                                                                              3)),
                                                                                          str(round(r3on,
                                                                                                    3)),
                                                                                          str(round(r3off,
                                                                                                    3)),
                                                                                          str(round(
                                                                                            tsleep3_k,
                                                                                            3)),
                                                                                          str(round(beta * abs(
                                                                                            delta3_k), 3)),
                                                                                          str(round(TCYCLE, 3)),
                                                                                          str(round(TON, 3)),
                                                                                          str(round(TOFF, 3)),
                                                                                          str(
                                                                                            round(
                                                                                              tsensors_on,
                                                                                              3)),
                                                                                          str(round(
                                                                                            tsensors_on_percent,
                                                                                            3)))

      a = a.expandtabs(8)
      file_all.write(a)

      b3 = '{0}\n'.format(str(round(delay3, 3)))
      b3 = b3.expandtabs(8)
      file_delay3.write(b3)

      # Node 4
      a = '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\n\n'.format(str(n4),
                                                                                                        str(
                                                                                                          round(
                                                                                                            expected4_tn,
                                                                                                            3)),
                                                                                                        str(
                                                                                                          round(
                                                                                                            sampled4_tn,
                                                                                                            3)),
                                                                                                        str(
                                                                                                          round(delay4,
                                                                                                                3)),
                                                                                                        str(
                                                                                                          round(
                                                                                                            expected4_tn - sampled4_tn,
                                                                                                            3)),
                                                                                                        str(
                                                                                                          round(
                                                                                                            delta4_k,
                                                                                                            3)),
                                                                                                        str(round(r4on,
                                                                                                                  3)),
                                                                                                        str(round(r4off,
                                                                                                                  3)),
                                                                                                        str(round(
                                                                                                          tsleep4_k,
                                                                                                          3)),
                                                                                                        str(round(
                                                                                                          beta * abs(
                                                                                                            delta4_k),
                                                                                                          3)),
                                                                                                        str(
                                                                                                          round(TCYCLE,
                                                                                                                3)),
                                                                                                        str(round(TON,
                                                                                                                  3)),
                                                                                                        str(round(TOFF,
                                                                                                                  3)),
                                                                                                        str(
                                                                                                          round(
                                                                                                            tsensors_on,
                                                                                                            3)),
                                                                                                        str(round(
                                                                                                          tsensors_on_percent,
                                                                                                          3)))

      a = a.expandtabs(8)
      file_all.write(a)

      b4 = '{0}\n'.format(str(round(delay4, 3)))
      b4 = b4.expandtabs(8)
      file_delay4.write(b4)

      c1 = '{0}\n'.format(str(round(tsensors_on, 3)))
      c1 = c1.expandtabs(8)
      file_tSensorsOn.write(c1)

      c2 = '{0}\n'.format(str(round(tsensors_on_percent, 3)))
      c2 = c2.expandtabs(8)
      file_tSensorsOn_percent.write(c2)

      d = '{0}\n'.format(str(n))
      d = d.expandtabs(8)
      file_n.write(d)

      e1 = '{0}\n'.format(str(round(beta * abs(delta1_k), 3)))
      e1 = e1.expandtabs(8)
      file_bdk.write(e1)

      e2 = '{0}\n'.format(str(round(beta * abs(delta2_k), 3)))
      e2 = e2.expandtabs(8)
      file_bdk.write(e2)

      e3 = '{0}\n'.format(str(round(beta * abs(delta3_k), 3)))
      e3 = e3.expandtabs(8)
      file_bdk.write(e3)

      e4 = '{0}\n'.format(str(round(beta * abs(delta4_k), 3)))
      e4 = e4.expandtabs(8)
      file_bdk.write(e4)

      if ((min(upper) - max(lower)) >= delta_success):
        hit += 1

      success = (hit / i_counter) * 100  # in percentage

    if (n >= int(TMAX * pDiscard)):
      a = '{0}\n'.format(str(round(success, 3)))
      a = a.expandtabs(8)
      file_tSensorsOn_success.write(a)

      a = '{0}\n'.format(str(round(TCYCLE, 3)))
      a = a.expandtabs(8)
      file_tcycle.write(a)

  a = '{0}\t{1}\n'.format('Total hit success', str(round(success, 3)))
  a = a.expandtabs(8)
  file_all.write(a)

  print 'Success hit for cycle %d: %.1f\n' % (j, success)
  temp_success = temp_success + success

mean_success = temp_success / success_counter
print 'Mean Success: %.1f\n' % (mean_success)

file_n.close()
file_all.close()
file_delay1.close()
file_delay2.close()
file_delay3.close()
file_tSensorsOn.close()
file_tSensorsOn_percent.close()
file_bdk.close()
file_tSensorsOn_success.close()

stop_t = time.clock()
print 'Stop processing at: %.3f' % (stop_t)
print 'total processing time: %.3f\n' % (stop_t - start_t)
# ***************************************************************************
