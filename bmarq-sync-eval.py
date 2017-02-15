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
'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
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
import subprocess
import argparse
import time
from pylab import *

__version__ = 'bmarq-sync-eval.py v1.0, (C)2017, Bruno Marques, INESC TEC, IPV/ESTGV, (bmarq@estgv.ipv.pt)'

global TON, TOFF, TCYCLE, MAXCICLES, TMAXCYCLE, TMINCYCLE, RNDSEED, NUMSIM
global tCycle1_n, tCycle2_n, tCycle3_n, tCycle4_n
global tCycle1_n_1, tCycle2_n_1, tCycle3_n_1, tCycle4_n_1
global tOn1_n, tOn2_n, tOn3_n, tOn4_n
global tOn1_n_1, tOn2_n_1, tOn3_n_1, tOn4_n_1
global tOff1_n, tOff2_n, tOff3_n, tOff4_n
global tOff1_n_1, tOff2_n_1, tOff3_n_1, tOff4_n_1
global rnd_delayDist, delayDist, rnd_cycleDist, cycleDist
global discard
global n1, n2, n3, n4
global delay_1, delay_2, delay_3, delay_4
global delay1, delay2, delay3, delay4
global delay1_n, delay2_n, delay3_n, delay4_n
global delay1_n_1, delay2_n_1, delay3_n_1, delay4_n_1
global expected1_n, expected2_n, expected3_n, expected4_n
global real1_n, real2_n, real3_n, real4_n
global real1_n_1, real2_n_1, real3_n_1, real4_n_1
global delta1_n, delta2_n, delta3_n, delta4_n
global delta1_n_1, delta2_n_1, delta3_n_1, delta4_n_1
global DELTA1_n, DELTA2_n, DELTA3_n, DELTA4_n
global DELTA1_n_1, DELTA2_n_1, DELTA3_n_1, DELTA4_n_1
global tsleep1_n, tsleep2_n, tsleep3_n, tsleep4_n
global r1on_n, r2on_n, r3on_n, r4on_n
global r1on_n_1, r2on_n_1, r3on_n_1, r4on_n_1
global r1off_n, r2off_n, r3off_n, r4off_n
global r1off_n_1, r2off_n_1, r3off_n_1, r4off_n_1
global t1i_n, t2i_n, t3i_n, t4i_n
global t1i_n_1, t2i_n_1, t3i_n_1, t4i_n_1
global t1f_n, t2f_n, t3f_n, t4f_n
global t1f_n_1, t2f_n_1, t3f_n_1, t4f_n_1
global sleepOffset1_n, sleepOffset2_n, sleepOffset3_n, sleepOffset4_n
global sleepOffset1_n_1, sleepOffset2_n_1, sleepOffset3_n_1, sleepOffset4_n_1
global status1_n, status2_n, status3_n, status4_n
global alpha, beta, gamma, sigma
global n_counter, hit
global lower, upper
global success, temp_success, delta_success, success_counter
global tsensors_on, tsensors_on_percent, tSensorsOnSuccess

# ***************************************************************************
parser = argparse.ArgumentParser(description='Evaluation of the bmarq-sync sycnhronization mechanism\n')
parser.add_argument('--version', action='version', version=__version__, help='Version number of eval_bmarq.py')
parser.add_argument('--numsim', type=int, default=1, help='Number of simulations to perform (default is 1)')
parser.add_argument('--maxcycles', type=long, default=1000,
                    help='Maximum number of cycles per simulation (default is 1000)')
parser.add_argument('--alpha', type=float, default=0.125, choices=sorted({0.01, 0.125, 0.25, 0.50, 0.75, 0.875, 0.99}),
                    help='Value for alpha parameter (defaultis 0.125)')
parser.add_argument('--beta', type=float, default=1, choices=sorted({1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5}),
                    help='Value for beta parameter (default is 1)')
parser.add_argument('--gamma', type=float, default=0.80, choices=sorted({0.5, 0.75, 0.80, 0.85, 0.90, 0.95}),
                    help='%% of TON for TSensorsOn success (default is 0.80)')
parser.add_argument('--sigma', type=float, default=0.20, choices=sorted({0.0, 0.01, 0.05, 0.10, 0.15, 0.20, 0.25}),
                    help='Value of standard deviation for the generated delays (default is 0.20)')
parser.add_argument('--ton', type=float, default=60, help='Value for TON (defaultis 60)')
parser.add_argument('--tmaxcycle', type=float, default=3600, help='Maximum value for TCycle (default is 3600)')
parser.add_argument('--tmincycle', type=float, default=120, help='Minimum value for TCycle (default is 120)')
parser.add_argument('--discard', type=float, default=0.10,
                    help='Initial %% of cycles to discard for estability purposes (default is 0.10 (10%%))')
parser.add_argument('--delaydist', type=str, default='uniform',
                    choices=sorted(
                      {'constant', 'uniform', 'normal', 'exponential', 'chisquare', 'poisson', 'pareto', 'weibull'}),
                    help='Type of random distribution for delays (default is uniform)')
parser.add_argument('--cycledist', type=str, default='uniform',
                    choices=sorted(
                      {'constant', 'uniform', 'normal', 'exponential', 'chisquare', 'poisson', 'pareto', 'weibull'}),
                    help='Type of random distribution for TCycle (default is uniform). If the distribution is constant, the default value equals TMINCYCLE')
parser.add_argument('--rndseed', type=str, default='T',
                    choices=sorted({'True', 'T', 'False', 'F'}),
                    help='Use predefined random seeds to reproduce experiments? (T)rue/(F)alse (default is (T)rue)')

args = parser.parse_args()
NUMSIM = args.numsim
MAXCICLES = args.maxcycles
alpha = args.alpha
beta = args.beta
gamma = args.gamma
sigma = args.sigma
TMAXCYCLE = args.tmaxcycle
TMINCYCLE = args.tmincycle
discard = args.discard
delayDist = args.delaydist
cycleDist = args.cycledist
TON = args.ton
RNDSEED = args.rndseed

# Initial values declaration
n1 = 1  # Node 1
n2 = 2  # Node 2
n3 = 3  # Node 3
n4 = 4  # Node 4

delta_success = gamma * TON

# ***************************************************************************
# Main block of the program
subprocess.call('clear', shell=True)  # clearing stdio
start_t = time.clock()

# if not os.path.exists('results'):
#  os.makedirs('results')

# if not os.path.exists('graphs'):
#  os.makedirs('graphs')

# ***************************************************************************
print 'Started processing at: %f ...\n' % (start_t)

# file_eval = open('./results/evaluation-parameters.txt', 'w')
f = (
      'number of simulations: %d\ncycle dist: %s\ndelay dist: %s\nalpha = %.3f\nbeta = %.1f\ngamma = %.2f\nsigma = %.2f\ntMinCycle = %.1f\ntMaxCycle = %.1f\nTON = %.1f\nMinimum time for simultaneous Sensors On for success consideration = %.2f (gamma * TON) \nnumber of cycles: %d\nIgnore first %.2f * 100 %%\n') % (
      NUMSIM, cycleDist, delayDist, alpha, beta, gamma, sigma, TMINCYCLE, TMAXCYCLE, TON, delta_success, MAXCICLES,
      discard)
# file_eval.write(f)
# file_eval.close()

print f

# ***************************************************************************
for j in range(1, int(NUMSIM) + 1):
  # Initial values declaration for the start of each simulation
  status1_n = status2_n = status3_n = status4_n = tSensorsOnSuccess = 'FAIL'

  delay1_n = delay_1 = delay1 = 0.5  # Initial delay for Node 1
  delay2_n = delay_2 = delay2 = 1.0  # Initial delay for Node 2
  delay3_n = delay_3 = delay3 = 2.0  # Initial delay for Node 3
  delay4_n = delay_4 = delay4 = 5.0  # Initial delay for Node 4

  delay1_n_1 = delay1_n
  delay2_n_1 = delay2_n
  delay3_n_1 = delay3_n
  delay4_n_1 = delay4_n

  delta_success = gamma * TON  # gamma % of TON
  temp_success = 0.00

  expected1_n = expected2_n = expected3_n = expected4_n = 0.000  # t'n for nodes 1, 2, 3 and 4 (expected reception time for packet n)
  expected1_n_1 = expected2_n_1 = expected3_n_1 = expected4_n_1 = 0.000

  real1_n = real2_n = real3_n = real4_n = 0.000  # tn-1 for nodes 1, 2, 3 and 4 (real reception time for packet n-1 [actual (n)])
  real1_n_1 = real2_n_1 = real3_n_1 = real4_n_1 = 0.000  # tn-1 for nodes 1, 2, 3 and 4 (real reception time for packet n-1 [anterior (n-1)])

  delta1_n = delta2_n = delta3_n = delta4_n = 0.000
  delta1_n_1 = delta2_n_1 = delta3_n_1 = delta4_n_1 = 0.000

  DELTA1_n = DELTA2_n = DELTA3_n = DELTA4_n = 0.000
  DELTA1_n_1 = DELTA2_n_1 = DELTA3_n_1 = DELTA4_n_1 = 0.000

  tsleep1_n = tsleep2_n = tsleep3_n = tsleep4_n = 0.000

  r1on_n = r2on_n = r3on_n = r4on_n = 0.000
  r1on_n_1 = r2on_n_1 = r3on_n_1 = r4on_n_1 = 0.000

  r1off_n = r2off_n = r3off_n = r4off_n = TON
  r1off_n_1 = r2off_n_1 = r3off_n_1 = r4off_n_1 = TON

  sleepOffset1_n = sleepOffset2_n = sleepOffset3_n = sleepOffset4_n = 0.000
  sleepOffset1_n_1 = sleepOffset2_n_1 = sleepOffset3_n_1 = sleepOffset4_n_1 = 0.000

  tOff1_n = tOff2_n = tOff3_n = tOff4_n = 0.000
  tOff1_n_1 = tOff2_n_1 = tOff3_n_1 = tOff4_n_1 = 0.000

  tCycle1_n = tCycle2_n = tCycle3_n = tCycle4_n = 0.000
  tCycle1_n_1 = tCycle2_n_1 = tCycle3_n_1 = tCycle4_n_1 = 0.000

  t1i_n = t2i_n = t3i_n = t4i_n = 0.0000
  t1i_n_1 = t2i_n_1 = t3i_n_1 = t4i_n_1 = 0.0000
  t1f_n = t2f_n = t3f_n = t4f_n = 0.0000
  t1f_n_1 = t2f_n_1 = t3f_n_1 = t4f_n_1 = 0.0000

  tsensors_on = tsensors_on_percent = 0.00
  success = success_counter = 0
  n_counter = 0
  hit = 0

  # ***************************************************************************
  file_n = open('./results/data-nCycles.txt', 'w')

  # file_all = open(
  #  './results/results-tcycle-dist_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
  #    gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  # file_sim_log = open(
  #  './results/data-log-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(gamma) + '-alpha_' + str(
  #    alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_delay1 = open(
    './results/data-delay-node1-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-sim_' + str(j) + '.txt', 'w')

  file_delay2 = open(
    './results/data-delay-node2-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-sim_' + str(j) + '.txt', 'w')

  file_delay3 = open(
    './results/data-delay-node3-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-sim_' + str(j) + '.txt', 'w')

  file_delay4 = open(
    './results/data-delay-node4-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-sim_' + str(j) + '.txt', 'w')

  file_tSensorsOn = open(
    './results/data-tsensors_on-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_tSensorsOn_percent = open(
    './results/data-tsensors_on-percent-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(
      j) + '.txt', 'w')

  file_tSensorsOn_success = open(
    './results/data-tsensors_on-success-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(
      j) + '.txt', 'w')

  file_tcycle = open(
    './results/data-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-sim_' + str(j) + '.txt', 'w')

  file_t_expected_node1 = open(
    './results/data-expected-node1-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_t_expected_node2 = open(
    './results/data-expected-node2-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_t_expected_node3 = open(
    './results/data-expected-node3-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_t_expected_node4 = open(
    './results/data-expected-node4-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_t_real_node1 = open(
    './results/data-real-node1-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_t_real_node2 = open(
    './results/data-real-node2-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_t_real_node3 = open(
    './results/data-real-node3-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_t_real_node4 = open(
    './results/data-real-node4-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_bdk_node1 = open(
    './results/data-sleepOffset-node1-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_bdk_node2 = open(
    './results/data-sleepOffset-node2-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_bdk_node3 = open(
    './results/data-sleepOffset-node3-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_bdk_node4 = open(
    './results/data-sleepOffset-node4-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_t_delta_node1 = open(
    './results/data-delta_n-node1-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_t_delta_node2 = open(
    './results/data-delta_n-node2-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_t_delta_node3 = open(
    './results/data-delta_n-node3-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_t_delta_node4 = open(
    './results/data-delta_n-node4-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_t_DELTA_node1 = open(
    './results/data-DELTA__n-node1-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_t_DELTA_node2 = open(
    './results/data-DELTA__n-node2-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_t_DELTA_node3 = open(
    './results/data-DELTA__n-node3-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  file_t_DELTA_node4 = open(
    './results/data-DELTA__n-node4-tcycle_' + cycleDist + '-delay_' + str(delayDist) + '-gamma_' + str(
      gamma) + '-alpha_' + str(alpha) + '-beta_' + str(beta) + '-sim_' + str(j) + '.txt', 'w')

  # ***************************************************************************
  a = 'Simulation\t%d\n' % (j)
  # a = a.expandtabs(4)
  # file_all.write(a)

  a = 'Cycle\tNode\tDelay\tReal\tExpect.\tDELTA\tdn\trOn\trOff\tTSleep\tb.|dn|\tTcycle\tTON\tTOFF\tStatus\tTSensorsOn\tTSensorsOn(%)\n'
  # a = a.expandtabs(4)
  # file_all.write(a)

  a = 'TON\t%.3f\tTSensorsOnSuccess\t%.3f\n\n' % (round(TON, 3), round((gamma * TON), 3))
  # a = a.expandtabs(4)
  # file_sim_log.write(a)

  # ***************************************************************************
  success_counter += 1
  n_counter = 1

  for n in range(0, int(MAXCICLES) + 1):
    a = 'Cycle\t%d\n' % (n)
    # a = a.expandtabs(4)
    # file_sim_log.write(a)

    # ***************************************************************************
    # generate delays for nodes for next cycle
    if (RNDSEED == 'T' or RNDSEED == 'True'):
      np.random.seed(j + n)
      np.random.RandomState(j + n)

    if (delayDist == 'uniform'):
      # node 1
      rnd_delayDist = 'np.random.' + delayDist + '((1 - sigma) * delay_1, (1 + sigma) * delay_1)'
      delay1 = eval(rnd_delayDist)

      # node 2
      rnd_delayDist = 'np.random.' + delayDist + '((1 - sigma) * delay_2, (1 + sigma) * delay_2)'
      delay2 = eval(rnd_delayDist)

      # node 3
      rnd_delayDist = 'np.random.' + delayDist + '((1 - sigma) * delay_3, (1 + sigma) * delay_3)'
      delay3 = eval(rnd_delayDist)

      # node 4
      rnd_delayDist = 'np.random.' + delayDist + '((1 - sigma) * delay_4, (1 + sigma) * delay_4)'
      delay4 = eval(rnd_delayDist)

    if (delayDist == 'normal'):
      # node 1
      rnd_delayDist = 'np.random.' + delayDist + '(delay_1, (1 + sigma) * delay_1)'
      delay1 = eval(rnd_delayDist)

      # node 2
      rnd_delayDist = 'np.random.' + delayDist + '(delay_2, (1 + sigma) * delay_2)'
      delay2 = eval(rnd_delayDist)

      # node 3
      rnd_delayDist = 'np.random.' + delayDist + '(delay_3, (1 + sigma) * delay_3)'
      delay3 = eval(rnd_delayDist)

      # node 4
      rnd_delayDist = 'np.random.' + delayDist + '(delay_4, (1 + sigma) * delay_4)'
      delay4 = eval(rnd_delayDist)

    if (delayDist == 'exponential'):
      # node 1
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_1)'
      delay1 = eval(rnd_delayDist)

      # node 2
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_2)'
      delay2 = eval(rnd_delayDist)

      # node 3
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_3)'
      delay3 = eval(rnd_delayDist)

      # node 4
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_4)'
      delay4 = eval(rnd_delayDist)

    if (delayDist == 'chisquare'):
      # node 1
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_1)'
      delay1 = eval(rnd_delayDist)

      # node 2
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_2)'
      delay2 = eval(rnd_delayDist)

      # node 3
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_3)'
      delay3 = eval(rnd_delayDist)

      # node 4
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_4)'
      delay4 = eval(rnd_delayDist)

    if (delayDist == 'poisson'):
      # node 1
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_1)'
      delay1 = eval(rnd_delayDist)

      # node 2
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_2)'
      delay2 = eval(rnd_delayDist)

      # node 3
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_3)'
      delay3 = eval(rnd_delayDist)

      # node 4
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_4)'
      delay4 = eval(rnd_delayDist)

    if (delayDist == 'pareto'):
      # node 1
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_1)'
      delay1 = eval(rnd_delayDist)

      # node 2
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_2)'
      delay2 = eval(rnd_delayDist)

      # node 3
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_3)'
      delay3 = eval(rnd_delayDist)

      # node 4
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_4)'
      delay4 = eval(rnd_delayDist)

    if (delayDist == 'weibull'):
      # node 1
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_1)'
      delay1 = eval(rnd_delayDist)

      # node 2
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_2)'
      delay2 = eval(rnd_delayDist)

      # node 3
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_3)'
      delay3 = eval(rnd_delayDist)

      # node 4
      rnd_delayDist = 'np.random.' + delayDist + '((1 + sigma) * delay_4)'
      delay4 = eval(rnd_delayDist)

    if (delayDist == 'constant'):
      # nodes 1, 2, 3, 4
      delay1 = delay_1
      delay2 = delay_2
      delay3 = delay_3
      delay4 = delay_4

    # ***************************************************************************
    # generate random TCycle next cycle
    if (RNDSEED == 'T' or RNDSEED == 'True'):
      np.random.seed(j + n)  # (n + 1)
      np.random.RandomState(j + n)  # (n + 1)

    if (cycleDist == 'uniform'):
      rnd_cycleDist = 'np.random.' + cycleDist + '((1 - sigma) * TMINCYCLE, (1 + sigma) * TMAXCYCLE)'
      TCYCLE = eval(rnd_cycleDist)

    if (cycleDist == 'normal'):
      rnd_cycleDist = 'np.random.' + cycleDist + '((1 - sigma) * TMINCYCLE, (1 + sigma) * TMAXCYCLE)'
      TCYCLE = eval(rnd_cycleDist)

    if (cycleDist == 'exponential'):
      rnd_cycleDist = 'np.random.' + cycleDist + '((1 + sigma) * TMAXCYCLE)'
      TCYCLE = eval(rnd_cycleDist)

    if (cycleDist == 'chisquare'):
      rnd_cycleDist = 'np.random.' + cycleDist + '((1 + sigma) * TMAXCYCLE)'
      TCYCLE = eval(rnd_cycleDist)

    if (cycleDist == 'poisson'):
      rnd_cycleDist = 'np.random.' + cycleDist + '((1 + sigma) * TMAXCYCLE)'
      TCYCLE = eval(rnd_cycleDist)

    if (cycleDist == 'pareto'):
      rnd_cycleDist = 'np.random.' + cycleDist + '((1 + sigma) * TMAXCYCLE)'
      TCYCLE = eval(rnd_cycleDist)

    if (cycleDist == 'weibull'):
      rnd_cycleDist = 'np.random.' + cycleDist + '((1 + sigma) * TMAXCYCLE)'
      TCYCLE = eval(rnd_cycleDist)

    if (cycleDist == 'constant'):
      TCYCLE = TMINCYCLE

    TOFF = TCYCLE - TON

    tOn1_n = tOn2_n = tOn3_n = tOn4_n = TON
    tOff1_n = tOff2_n = tOff3_n = tOff4_n = TOFF

    tCycle1_n = tCycle2_n = tCycle3_n = tCycle4_n = TCYCLE
    tCycle1_n_1 = tCycle2_n_1 = tCycle3_n_1 = tCycle4_n_1 = TCYCLE

    delay1_n = delay1
    delay2_n = delay2
    delay3_n = delay3
    delay4_n = delay4

    if (n == 0):
      tOn1_n_1 = tOn2_n_1 = tOn3_n_1 = tOn4_n_1 = 0
      r1off_n = tOn1_n
      r2off_n = tOn2_n
      r3off_n = tOn3_n
      r4off_n = tOn4_n
      t1i_n = t2i_n = t3i_n = t3i_n = 0
      t1f_n = tCycle1_n
      t2f_n = tCycle2_n
      t3f_n = tCycle3_n
      t4f_n = tCycle4_n

    else:
      # **********************************************************
      # Node 1
      t1i_n = t1i_n_1 + tCycle1_n_1
      real1_n = t1i_n + delay1
      expected1_n = t1i_n + delta1_n_1

      DELTA1_n = expected1_n - real1_n
      delta1_n = (1.0 - alpha) * delta1_n_1 + alpha * abs(DELTA1_n)
      sleepOffset1_n = beta * abs(delta1_n)
      tsleep1_n = tOff1_n - sleepOffset1_n

      r1on_n = t1i_n - sleepOffset1_n_1
      t1f_n = r1on_n + tOn1_n
      r1off_n = t1f_n

      if (r1on_n < real1_n or r1on_n < expected1_n):
        status1_n = 'PASS'
      else:
        status1_n = 'FAIL'

      # **********************************************************
      # Node 2
      t2i_n = t2i_n_1 + tCycle2_n_1
      real2_n = t2i_n + delay2
      expected2_n = t2i_n + delta2_n_1

      DELTA2_n = expected2_n - real2_n
      delta2_n = (1.0 - alpha) * delta2_n_1 + alpha * abs(DELTA2_n)
      sleepOffset2_n = beta * abs(delta2_n)
      tsleep2_n = tOff2_n - sleepOffset2_n

      r2on_n = t2i_n - sleepOffset2_n_1
      t2f_n = r2on_n + tOn2_n
      r2off_n = t2f_n

      if (r2on_n < real2_n or r2on_n < expected2_n):
        status2_n = 'PASS'
      else:
        status2_n = 'FAIL'

      # **********************************************************
      # Node 3
      t3i_n = t3i_n_1 + tCycle3_n_1
      real3_n = t3i_n + delay3
      expected3_n = t3i_n + delta3_n_1

      DELTA3_n = expected3_n - real3_n
      delta3_n = (1.0 - alpha) * delta3_n_1 + alpha * abs(DELTA3_n)
      sleepOffset3_n = beta * abs(delta3_n)
      tsleep3_n = tOff3_n - sleepOffset3_n

      r3on_n = t3i_n - sleepOffset3_n_1
      t3f_n = r3on_n + tOn3_n
      r3off_n = t3f_n

      if (r3on_n < real3_n or r3on_n < expected3_n):
        status3_n = 'PASS'
      else:
        status3_n = 'FAIL'

      # **********************************************************
      # Node 4
      t4i_n = t4i_n_1 + tCycle4_n_1
      real4_n = t4i_n + delay4
      expected4_n = t4i_n + delta4_n_1

      DELTA4_n = expected4_n - real4_n
      delta4_n = (1.0 - alpha) * delta4_n_1 + alpha * abs(DELTA4_n)
      sleepOffset4_n = beta * abs(delta4_n)
      tsleep4_n = tOff4_n - sleepOffset4_n

      r4on_n = t4i_n - sleepOffset4_n_1
      t4f_n = r4on_n + tOn4_n
      r4off_n = t4f_n

      if (r4on_n < real4_n or r4on_n < expected4_n):
        status4_n = 'PASS'
      else:
        status4_n = 'FAIL'

    # ***************************************************************************
    # Node 1
    delta1_n_1 = delta1_n
    DELTA1_n_1 = DELTA1_n
    real1_n_1 = real1_n
    expected1_n_1 = expected1_n
    sleepOffset1_n_1 = sleepOffset1_n
    tCycle1_n_1 = tCycle1_n
    tOn1_n_1 = tOn1_n
    tOff1_n_1 = tOff1_n
    r1on_n_1 = r1on_n
    r1off_n_1 = r1off_n
    delay1_n_1 = delay1_n
    t1i_n_1 = t1i_n
    t1f_n_1 = t1f_n

    # ***************************************************************************
    # Node 2
    delta2_n_1 = delta2_n
    DELTA2_n_1 = DELTA2_n
    real2_n_1 = real2_n
    expected2_n_1 = expected2_n
    sleepOffset2_n_1 = sleepOffset2_n
    tCycle2_n_1 = tCycle2_n
    tOn2_n_1 = tOn2_n
    tOff2_n_1 = tOff2_n
    r2on_n_1 = r2on_n
    r2off_n_1 = r2off_n
    delay2_n_1 = delay2_n
    t2i_n_1 = t2i_n
    t2f_n_1 = t2f_n

    # ***************************************************************************
    # Node 3
    delta3_n_1 = delta3_n
    DELTA3_n_1 = DELTA3_n
    real3_n_1 = real3_n
    expected3_n_1 = expected3_n
    sleepOffset3_n_1 = sleepOffset3_n
    tCycle3_n_1 = tCycle3_n
    tOn3_n_1 = tOn3_n
    tOff3_n_1 = tOff3_n
    r3on_n_1 = r3on_n
    r3off_n_1 = r3off_n
    delay3_n_1 = delay3_n
    t3i_n_1 = t3i_n
    t3f_n_1 = t3f_n

    # ***************************************************************************
    # Node 4
    delta4_n_1 = delta4_n
    DELTA4_n_1 = DELTA4_n
    real4_n_1 = real4_n
    expected4_n_1 = expected4_n
    sleepOffset4_n_1 = sleepOffset4_n
    tCycle4_n_1 = tCycle4_n
    tOn4_n_1 = tOn4_n
    tOff4_n_1 = tOff4_n
    r4on_n_1 = r4on_n
    r4off_n_1 = r4off_n
    delay4_n_1 = delay4_n
    t4i_n_1 = t4i_n
    t4f_n_1 = t4f_n

    # ***************************************************************************
    lower = (r1on_n, r2on_n, r3on_n, r4on_n)
    upper = (r1off_n, r2off_n, r3off_n, r4off_n)

    tsensors_on = (min(upper) - max(lower))
    tsensors_on_percent = (tsensors_on / TON) * 100

    if ((tsensors_on >= (gamma * TON)) and (tsensors_on <= TON)):
      tSensorsOnSuccess = 'Success'
    else:
      tSensorsOnSuccess = 'Fail'

    a = 'TCYCLE=\t%.3f\tTOFF\t%.3f\n' % (round(TCYCLE, 3), round(TOFF, 3))
    # a = a.expandtabs(4)
    # file_sim_log.write(a)

    a = 'r1on_n\t%.3f\tr1off_n\t%.3f\n' % (round(r1on_n, 3), round(r1off_n, 3))
    # file_sim_log.write(a)

    a = 'r2on_n\t%.3f\tr2off_n\t%.3f\n' % (round(r2on_n, 3), round(r2off_n, 3))
    # a = a.expandtabs(4)
    # file_sim_log.write(a)

    a = 'r3on_n\t%.3f\tr3off_n\t%.3f\n' % (round(r3on_n, 3), round(r3off_n, 3))
    # a = a.expandtabs(4)
    # file_sim_log.write(a)

    a = 'r4on_n\t%.3f\tr4off_n\t%.3f\n' % (round(r4on_n, 3), round(r4off_n, 3))
    # a = a.expandtabs(4)
    # file_sim_log.write(a)

    a = 'lower_t\t%.3f\tupper_t\t%.3f\ttSensorsOn(s)\t%.3f\ttSensorsOn(%%)\t%.3f\t%s\n\n' % (
      round(max(lower), 3), round(min(upper), 3), round(tsensors_on, 3), round(tsensors_on_percent, 3),
      tSensorsOnSuccess)
    # a = a.expandtabs(4)
    # file_sim_log.write(a)

    if (n >= int(MAXCICLES * discard)):
      n_counter += 1

      if ((min(upper) - max(lower)) >= delta_success):
        hit += 1

      success = (hit / (n_counter - 1)) * 100  # in percentage

      a = '{0}\n'.format(str(round(success, 3)))
      # a = a.expandtabs(4)
      file_tSensorsOn_success.write(a)

      # Node 1
      a = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\n'.format(str(n),
                                                                                                    str(n1),
                                                                                                    str(
                                                                                                      round(delay1_n,
                                                                                                            3)),
                                                                                                    str(
                                                                                                      round(real1_n,
                                                                                                            3)),
                                                                                                    str(
                                                                                                      round(
                                                                                                        expected1_n,
                                                                                                        3)),
                                                                                                    str(
                                                                                                      round(DELTA1_n,
                                                                                                            3)),
                                                                                                    str(
                                                                                                      round(
                                                                                                        delta1_n,
                                                                                                        3)),
                                                                                                    str(
                                                                                                      round(r1on_n, 3)),
                                                                                                    str(round(r1off_n,
                                                                                                              3)),
                                                                                                    str(round(
                                                                                                      tsleep1_n,
                                                                                                      3)),
                                                                                                    str(
                                                                                                      round(
                                                                                                        sleepOffset1_n,
                                                                                                        3)),
                                                                                                    str(
                                                                                                      round(TCYCLE, 3)),
                                                                                                    str(round(TON, 3)),
                                                                                                    str(round(TOFF, 3)),
                                                                                                    status1_n)

      # a = a.expandtabs(4)
      # file_all.write(a)

      b1 = '{0}\n'.format(str(round(delay1_n, 3)))
      # b1 = b1.expandtabs(4)
      file_delay1.write(b1)

      c1 = '{0}\n'.format(str(round(sleepOffset1_n, 3)))
      # c1 = c1.expandtabs(4)
      file_bdk_node1.write(c1)

      d1 = '{0}\n'.format(str(round(delta1_n, 3)))
      # d1 = d1.expandtabs(4)
      file_t_delta_node1.write(d1)

      e1 = '{0}\n'.format(str(round(DELTA1_n, 3)))
      # e1 = e1.expandtabs(4)
      file_t_DELTA_node1.write(e1)

      f1 = '{0}\n'.format(str(round(expected1_n, 3)))
      # f1 = f1.expandtabs(4)
      file_t_expected_node1.write(f1)

      g1 = '{0}\n'.format(str(round(real1_n, 3)))
      # g1 = g1.expandtabs(4)
      file_t_real_node1.write(g1)


      # Node 2
      a = '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n'.format(str(n2),
                                                                                                str(round(delay2, 3)),
                                                                                                str(
                                                                                                  round(real2_n, 3)),
                                                                                                str(
                                                                                                  round(expected2_n,
                                                                                                        3)),
                                                                                                str(
                                                                                                  round(DELTA2_n,
                                                                                                        3)),
                                                                                                str(
                                                                                                  round(
                                                                                                    delta2_n,
                                                                                                    3)),
                                                                                                str(round(r2on_n,
                                                                                                          3)),
                                                                                                str(round(r2off_n,
                                                                                                          3)),
                                                                                                str(round(
                                                                                                  tsleep2_n,
                                                                                                  3)),
                                                                                                str(
                                                                                                  round(sleepOffset2_n,
                                                                                                        3)),
                                                                                                str(round(TCYCLE, 3)),
                                                                                                str(round(TON, 3)),
                                                                                                str(round(TOFF, 3)),
                                                                                                status2_n)

      # a = a.expandtabs(4)
      # file_all.write(a)

      b2 = '{0}\n'.format(str(round(delay2, 3)))
      # b2 = b2.expandtabs(4)
      file_delay2.write(b2)

      c2 = '{0}\n'.format(str(round(sleepOffset2_n, 3)))
      # c2 = c2.expandtabs(4)
      file_bdk_node2.write(c2)

      d2 = '{0}\n'.format(str(round(delta2_n, 3)))
      # d2 = d2.expandtabs(4)
      file_t_delta_node2.write(d2)

      e2 = '{0}\n'.format(str(round(DELTA2_n, 3)))
      # e2 = e2.expandtabs(4)
      file_t_DELTA_node2.write(e2)

      f2 = '{0}\n'.format(str(round(expected2_n, 3)))
      # f2 = f2.expandtabs(4)
      file_t_expected_node2.write(f2)

      g2 = '{0}\n'.format(str(round(real2_n, 3)))
      # g2 = g2.expandtabs(4)
      file_t_real_node2.write(g2)


      # Node 3
      a = '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n'.format(str(n3),
                                                                                                str(round(delay3, 3)),
                                                                                                str(round(real3_n, 3)),
                                                                                                str(round(expected3_n,
                                                                                                          3)),
                                                                                                str(
                                                                                                  round(DELTA3_n, 3)),
                                                                                                str(round(delta3_n, 3)),
                                                                                                str(round(r3on_n, 3)),
                                                                                                str(round(r3off_n, 3)),
                                                                                                str(
                                                                                                  round(tsleep3_n, 3)),
                                                                                                str(
                                                                                                  round(sleepOffset3_n,
                                                                                                        3)),
                                                                                                str(round(TCYCLE, 3)),
                                                                                                str(round(TON, 3)),
                                                                                                str(round(TOFF, 3)),
                                                                                                status3_n)

      # a = a.expandtabs(4)
      # file_all.write(a)

      b3 = '{0}\n'.format(str(round(delay3, 3)))
      # b3 = b3.expandtabs(4)
      file_delay3.write(b3)

      c3 = '{0}\n'.format(str(round(sleepOffset3_n, 3)))
      # c3 = c3.expandtabs(4)
      file_bdk_node3.write(c3)

      d3 = '{0}\n'.format(str(round(delta3_n, 3)))
      # d3 = d3.expandtabs(4)
      file_t_delta_node3.write(d3)

      e3 = '{0}\n'.format(str(round(DELTA3_n, 3)))
      # e3 = e3.expandtabs(4)
      file_t_DELTA_node3.write(e3)

      f3 = '{0}\n'.format(str(round(expected3_n, 3)))
      # f3 = f3.expandtabs(4)
      file_t_expected_node3.write(f3)

      g3 = '{0}\n'.format(str(round(real3_n, 3)))
      # g3 = g3.expandtabs(4)
      file_t_real_node3.write(g3)

      # Node 4
      a = '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\n\n'.format(str(n4),
                                                                                                              str(
                                                                                                                round(
                                                                                                                  delay4,
                                                                                                                  3)),
                                                                                                              str(
                                                                                                                round(
                                                                                                                  real4_n,
                                                                                                                  3)),
                                                                                                              str(
                                                                                                                round(
                                                                                                                  expected4_n,
                                                                                                                  3)),
                                                                                                              str(
                                                                                                                round(
                                                                                                                  DELTA4_n,
                                                                                                                  3)),
                                                                                                              str(
                                                                                                                round(
                                                                                                                  delta4_n,
                                                                                                                  3)),
                                                                                                              str(round(
                                                                                                                r4on_n,
                                                                                                                3)),
                                                                                                              str(round(
                                                                                                                r4off_n,
                                                                                                                3)),
                                                                                                              str(round(
                                                                                                                tsleep4_n,
                                                                                                                3)),
                                                                                                              str(round(
                                                                                                                sleepOffset4_n,
                                                                                                                3)),
                                                                                                              str(
                                                                                                                round(
                                                                                                                  TCYCLE,
                                                                                                                  3)),
                                                                                                              str(
                                                                                                                round(
                                                                                                                  TON,
                                                                                                                  3)),
                                                                                                              str(round(
                                                                                                                TOFF,
                                                                                                                3)),
                                                                                                              status4_n,
                                                                                                              str(
                                                                                                                round(
                                                                                                                  tsensors_on,
                                                                                                                  3)),
                                                                                                              str(round(
                                                                                                                tsensors_on_percent,
                                                                                                                3)))

      # a = a.expandtabs(4)
      # file_all.write(a)

      b4 = '{0}\n'.format(str(round(delay4, 3)))
      # b4 = b4.expandtabs(4)
      file_delay4.write(b4)

      d4 = '{0}\n'.format(str(round(sleepOffset4_n, 3)))
      # d4 = d4.expandtabs(4)
      file_bdk_node4.write(d4)

      e4 = '{0}\n'.format(str(round(delta4_n, 3)))
      # e4 = e4.expandtabs(4)
      file_t_delta_node4.write(e4)

      f4 = '{0}\n'.format(str(round(DELTA4_n, 3)))
      # f4 = f4.expandtabs(4)
      file_t_DELTA_node4.write(f4)

      f4 = '{0}\n'.format(str(round(expected4_n, 3)))
      # f4 = f4.expandtabs(4)
      file_t_expected_node4.write(f4)

      g4 = '{0}\n'.format(str(round(real4_n, 3)))
      # g4 = g4.expandtabs(4)
      file_t_real_node4.write(g4)

      # ***************************************************************************
      f = '{0}\n'.format(str(n))
      f = f.expandtabs(4)
      file_n.write(f)

      a = '{0}\n'.format(str(round(TCYCLE, 3)))
      # a = a.expandtabs(4)
      file_tcycle.write(a)

      c1 = '{0}\n'.format(str(round(tsensors_on, 3)))
      # c1 = c1.expandtabs(4)
      file_tSensorsOn.write(c1)

      c2 = '{0}\n'.format(str(round(tsensors_on_percent, 3)))
      # c2 = c2.expandtabs(4)
      file_tSensorsOn_percent.write(c2)

  a = '{0}\t{1}\n'.format('Total hit success', str(round(success, 3)))
  # a = a.expandtabs(4)
  # file_all.write(a)

  print 'Simulation %d: hit = %d, n. cycles = %d, success = %.3f' % (j, hit, (n_counter - 1), success)

  file_n.close()
  # file_sim_log.close()
  # file_all.close()

  file_t_expected_node1.close()
  file_t_expected_node2.close()
  file_t_expected_node2.close()
  file_t_expected_node2.close()

  file_t_real_node1.close()
  file_t_real_node2.close()
  file_t_real_node3.close()
  file_t_real_node4.close()

  file_delay1.close()
  file_delay2.close()
  file_delay3.close()
  file_delay4.close()

  file_tSensorsOn.close()
  file_tSensorsOn_percent.close()
  file_tSensorsOn_success.close()

  file_bdk_node1.close()
  file_bdk_node2.close()
  file_bdk_node3.close()
  file_bdk_node4.close()

  file_t_delta_node1.close()
  file_t_delta_node2.close()
  file_t_delta_node3.close()
  file_t_delta_node4.close()

  file_t_DELTA_node1.close()
  file_t_DELTA_node2.close()
  file_t_DELTA_node3.close()
  file_t_DELTA_node4.close()

stop_t = time.clock()
print 'Stop processing at: %.3f' % (stop_t)
print 'Total processing time: %.3f\n' % (stop_t - start_t)
# ***************************************************************************
