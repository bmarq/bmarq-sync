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
import subprocess
import argparse
import time

from pylab import *

__version__ = "eval_bmarq.py v1.0, (C) 2017, Bruno Marques, INESC TEC, IPV/ESTGV"

global TON, TOFF, TCYCLE, MAXCICLES, TMAXCYCLE, TMINCYCLE
global tCycle1_n, tCycle2_n, tCycle3_n, tCycle4_n
global tCycle1_n_1, tCycle2_n_1, tCycle3_n_1, tCycle4_n_1
global tOff1_n_1, tOff2_n_1, tOff3_n_1, tOff4_n_1
global rnd_tDelayDist, tDelayDist, rnd_tCycleDist, tCycleDist
global NSIM
global pDiscard
global n1, n2, n3, n4
global delay_1, delay_2, delay_3, delay_4
global delay1, delay2, delay3, delay4
global expected1_tn, expected2_tn, expected3_tn, expected4_tn
global real1_tn, real2_tn, real3_tn, real4_tn
global real1_tn_1, real2_tn_1, real3_tn_1, real4_tn_1
global delta1_n, delta2_n, delta3_n, delta4_n
global delta1_n_1, delta2_n_1, delta3_n_1, delta4_n_1
global DELTA1_n, DELTA2_n, DELTA3_n, DELTA4_n
global DELTA1_n_1, DELTA2_n_1, DELTA3_n_1, DELTA4_n_1
global tsleep1_n, tsleep2_n, tsleep3_n, tsleep4_n
global r1on_n, r2on_n, r3on_n, r4on_n
global r1on_n_1, r2on_n_1, r3on_n_1, r4on_n_1
global r1off_n, r2off_n, r3off_n, r4off_n
global r1off_n_1, r2off_n_1, r3off_n_1, r4off_n_1
global sleepOffset1_n, sleepOffset2_n, sleepOffset3_n, sleepOffset4_n
global sleepOffset1_n_1, sleepOffset2_n_1, sleepOffset3_n_1, sleepOffset4_n_1
global status1_n, status2_n, status3_n, status4_n
global alpha, beta, gamma, sigma
global n_counter, hit
global lower, upper
global success, temp_success, delta_success, success_counter
global tsensors_on, tsensors_on_percent

# ***************************************************************************
parser = argparse.ArgumentParser(description='Evaluation of the bmarq-sync sycnhronization mechanism\n')
parser.add_argument('--version', action='version', version=__version__, help='Version number of eval_bmarq.py')
parser.add_argument('--nsim', type=int, default=1, help='Number of simulations to perform (default: 1)')
parser.add_argument('--maxcicles', type=int, default=100,
                    help='Maximum number of cycles per simulation (default: 10)')
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
                    help='Type of random distribution for TCycle (default: uniform). If the distribution is constant, the default value equals TMINCYCLE')

args = parser.parse_args()
NSIM = args.nsim
MAXCICLES = args.maxcicles
alpha = args.alpha  # test with 0.125; 0.50; 0.875
beta = args.beta  # test with 1, 10, 50, 100
gamma = args.gamma  # % of TON for TsensorsOn success
sigma = args.sigma  # value for delays standard deviation
TMAXCYCLE = args.tmaxcycle
TMINCYCLE = args.tmincycle
pDiscard = args.pdiscard
tDelayDist = args.tdelaydist
tCycleDist = args.tcycledist
TON = args.ton

# Initial values declaration
n1 = 1  # Node 1
n2 = 2  # Node 2
n3 = 3  # Node 3
n4 = 4  # Node 4

status1_n = status2_n = status3_n = status4_n = 'FAIL'

delay1 = delay_1 = 0.500  # Initial delay for Node 1
delay2 = delay_2 = 1.000  # Initial delay for Node 2
delay3 = delay_3 = 2.000  # Initial delay for Node 3
delay4 = delay_4 = 5.000  # Initial delay for Node 4

delta_success = gamma * TON  # gamma % of TON

temp_success = 0.00

'''
# ***************************************************************************
expected1_tn = expected2_tn = expected3_tn = expected4_tn = 0.000  # t'k for nodes 1, 2, 3 and 4 (expected reception time for packet k)
expected1_tn_1 = expected2_tn_1 = expected3_tn_1 = expected4_tn_1 = 0.000

real1_tn = real2_tn = real3_tn = real4_tn = 0.000  # tk -1 for nodes 1, 2, 3 and 4 (real reception time for packet k-1 [actual (k)])
real1_tn_1 = real2_tn_1 = real3_tn_1 = real4_tn_1 = 0.000  # tk -1 for nodes 1, 2, 3 and 4 (real reception time for packet k-1 [anterior (k-1)])

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

tOff1_n_1 = tOff2_n_1 = tOff3_n_1 = tOff4_n_1 = 0.000

tCycle1_n = tCycle2_n = tCycle3_n = tCycle4_n = 0.000
tCycle1_n_1 = tCycle2_n_1 = tCycle3_n_1 = tCycle4_n_1 = 0.000

tsensors_on = tsensors_on_percent = 0.00

TCYCLE = TMINCYCLE
TOFF = TCYCLE - TON

success = success_counter = temp_success = 0
n_counter = 0
hit = 0
'''

# ***************************************************************************
# Main block of the program
subprocess.call('clear', shell=True)  # clearing stdio
start_t = time.clock()

# ***************************************************************************
print 'Started processing at: %f ...\n' % (start_t)

file_eval = open('../results/evaluation-readme.txt', 'w')
f = (
      'number of simulations: %d\ntcycle dist: %s\ndelay dist: %s\nalpha = %.3f\nbeta = %.3f\ngamma = %.3f\nsigma = %.3f\nTON = %.3f\nMinimum time for simultaneous Sensors On for success consideration = %.3f (gamma * alpha) \nnumber of cycles: %d\nIgnore first %.3f * 100 %%\n') % (
      NSIM, tCycleDist, tDelayDist, alpha, beta, gamma, sigma, TON, delta_success, MAXCICLES, pDiscard)
file_eval.write(f)
file_eval.close()

print f

file_tSensorsOn_success = open(
  '../results/data-tsensors_on-success-tcycle_' + tCycleDist + '-delay_' + str(tDelayDist) + '.txt', 'w')

# ***************************************************************************
for j in range(1, int(NSIM) + 1):
  # ***************************************************************************
  expected1_tn = expected2_tn = expected3_tn = expected4_tn = 0.000  # t'k for nodes 1, 2, 3 and 4 (expected reception time for packet k)
  expected1_tn_1 = expected2_tn_1 = expected3_tn_1 = expected4_tn_1 = 0.000

  real1_tn = real2_tn = real3_tn = real4_tn = 0.000  # tk -1 for nodes 1, 2, 3 and 4 (real reception time for packet k-1 [actual (k)])
  real1_tn_1 = real2_tn_1 = real3_tn_1 = real4_tn_1 = 0.000  # tk -1 for nodes 1, 2, 3 and 4 (real reception time for packet k-1 [anterior (k-1)])

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

  tOff1_n_1 = tOff2_n_1 = tOff3_n_1 = tOff4_n_1 = 0.000

  tCycle1_n = tCycle2_n = tCycle3_n = tCycle4_n = 0.000
  tCycle1_n_1 = tCycle2_n_1 = tCycle3_n_1 = tCycle4_n_1 = 0.000

  tsensors_on = tsensors_on_percent = 0.00

  TCYCLE = TMINCYCLE
  TOFF = TCYCLE - TON

  success = success_counter = 0
  n_counter = 0
  hit = 0

  # ***************************************************************************
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

  file_tcycle = open(
    '../results/data-tcycle-tcycle_' + tCycleDist + '-delay_' + str(tDelayDist) + '-sim_' + str(j) + '.txt',
    'w')

  # file_bdk = open('../results/data-bdk-tcycle_' + tCycleDist + '-delay_' + str(tDelayDist) + '-sim_' + str(j) + '.txt',
  #                'w')

  file_bdk_node1 = open(
    '../results/data-sleepOffset-node1-tcycle_' + tCycleDist + '-delay_' + str(tDelayDist) + '-sim_' + str(j) + '.txt',
    'w')

  file_bdk_node2 = open(
    '../results/data-sleepOffset-node2-tcycle_' + tCycleDist + '-delay_' + str(tDelayDist) + '-sim_' + str(j) + '.txt',
    'w')

  file_bdk_node3 = open(
    '../results/data-sleepOffset-node3-tcycle_' + tCycleDist + '-delay_' + str(tDelayDist) + '-sim_' + str(j) + '.txt',
    'w')

  file_bdk_node4 = open(
    '../results/data-sleepOffset-node4-tcycle_' + tCycleDist + '-delay_' + str(tDelayDist) + '-sim_' + str(j) + '.txt',
    'w')

  # ***************************************************************************
  a = 'Cycle\tNode\tReal\tExpect.\tDelay\tDELTA\tdn\trOn\trOff\tTSleep\tb.|dn|\tTcycle\tTON\tTOFF\tStatus\tTSensorsOn\tTSensorsOn(%)\n'
  a = a.expandtabs(8)
  file_all.write(a)

  # ***************************************************************************
  success_counter += 1
  n_counter = 1

  for n in range(0, int(MAXCICLES) + 1):
    print
    print 'Simulation\t%d\tCycle\t%d:' % (j, n)

    if (n == 0):
      # initialization for nodes
      # Node 1
      r1on_n = 0.000
      r1off_n = TON
      real1_tn = r1on_n + delay1
      expected1_tn = real1_tn
      tsleep1_n = TOFF - sleepOffset1_n
      tCycle1_n_1 = TCYCLE
      tOff1_n_1 = TOFF

      # Node 2
      r2on_n = 0.000
      r2off_n = TON
      real2_tn = r2on_n + delay2
      expected2_tn = real2_tn
      tsleep2_n = TOFF - sleepOffset2_n
      tCycle2_n_1 = TCYCLE
      tOff2_n_1 = TOFF

      # Node 3
      r3on_n = 0.000
      r3off_n = TON
      real3_tn = r3on_n + delay3
      expected3_tn = real3_tn
      tsleep3_n = TOFF - sleepOffset3_n
      tCycle3_n_1 = TCYCLE
      tOff3_n_1 = TOFF

      # Node 4
      r4on_n = 0.000
      r4off_n = TON
      real4_tn = r4on_n + delay4
      expected4_tn = real4_tn
      tsleep4_n = TOFF - sleepOffset4_n
      tCycle4_n_1 = TCYCLE
      tOff4_n_1 = TOFF

    else:
      # adusting time and cicles for nodes
      # Node 1
      tCycle1_n = TCYCLE + sleepOffset1_n_1
      r1on_n = r1on_n_1 + tCycle1_n_1 - sleepOffset1_n_1
      r1off_n = r1on_n + TON

      real1_tn = r1on_n + delay1
      expected1_tn = r1on_n_1 + tCycle1_n_1

      DELTA1_n = expected1_tn - real1_tn
      delta1_n = (1.0 - alpha) * delta1_n_1 + alpha * DELTA1_n
      sleepOffset1_n = beta * abs(delta1_n)
      tsleep1_n = TOFF - sleepOffset1_n

      delta1_n_1 = delta1_n
      DELTA1_n_1 = DELTA1_n
      real1_tn_1 = real1_tn
      expected1_tn_1 = expected1_tn
      sleepOffset1_n_1 = sleepOffset1_n
      tCycle1_n_1 = tCycle1_n
      tOff1_n_1 = TOFF
      r1on_n_1 = r1on_n

      if (r1on_n < real1_tn or r1on_n <= expected1_tn):
        status1_n = 'PASS'
      else:
        status1_n = 'FAIL'


      # Node 2
      tCycle2_n = TCYCLE + sleepOffset2_n_1
      r2on_n = r2on_n_1 + tCycle2_n_1 - sleepOffset2_n_1
      r2off_n = r2on_n + TON

      real2_tn = r2on_n + delay2
      expected2_tn = r2on_n_1 + tCycle2_n_1

      DELTA2_n = expected2_tn - real2_tn
      delta2_n = (1.0 - alpha) * delta2_n_1 + alpha * DELTA2_n
      sleepOffset2_n = beta * abs(delta2_n)
      tsleep2_n = TOFF - sleepOffset2_n

      delta2_n_1 = delta2_n
      DELTA2_n_1 = DELTA2_n
      real2_tn_1 = real2_tn
      expected2_tn_1 = expected2_tn
      sleepOffset2_n_1 = sleepOffset2_n
      tCycle2_n_1 = tCycle2_n
      tOff2_n_1 = TOFF
      r2on_n_1 = r2on_n

      if (r2on_n < real2_tn or r2on_n <= expected2_tn):
        status2_n = 'PASS'
      else:
        status2_n = 'FAIL'


      # Node 3
      tCycle3_n = TCYCLE + sleepOffset3_n_1
      r3on_n = r3on_n_1 + tCycle3_n_1 - sleepOffset3_n_1
      r3off_n = r3on_n + TON

      real3_tn = r3on_n + delay3
      expected3_tn = r3on_n_1 + tCycle3_n_1

      DELTA3_n = expected3_tn - real3_tn
      delta3_n = (1.0 - alpha) * delta3_n_1 + alpha * DELTA3_n
      sleepOffset3_n = beta * abs(delta3_n)
      tsleep3_n = TOFF - sleepOffset3_n

      delta3_n_1 = delta3_n
      DELTA3_n_1 = DELTA3_n
      real3_tn_1 = real3_tn
      expected3_tn_1 = expected3_tn
      sleepOffset3_n_1 = sleepOffset3_n
      tCycle3_n_1 = tCycle3_n
      tOff3_n_1 = TOFF
      r3on_n_1 = r3on_n

      if (r3on_n < real3_tn or r3on_n <= expected3_tn):
        status3_n = 'PASS'
      else:
        status3_n = 'FAIL'

      # Node 4
      tCycle4_n = TCYCLE + sleepOffset4_n_1
      r4on_n = r4on_n_1 + tCycle4_n_1 - sleepOffset4_n_1
      r4off_n = r4on_n + TON

      real4_tn = r4on_n + delay4
      expected4_tn = r4on_n_1 + tCycle4_n_1

      DELTA4_n = expected4_tn - real4_tn
      delta4_n = (1.0 - alpha) * delta4_n_1 + alpha * DELTA4_n
      sleepOffset4_n = beta * abs(delta4_n)
      tsleep4_n = TOFF - sleepOffset4_n

      delta4_n_1 = delta4_n
      DELTA4_n_1 = DELTA4_n
      real4_tn_1 = real4_tn
      expected4_tn_1 = expected4_tn
      sleepOffset4_n_1 = sleepOffset4_n

      # tCycle4_n_1 = TCYCLE
      tCycle4_n_1 = tCycle4_n

      tOff4_n_1 = TOFF
      r4on_n_1 = r4on_n

      if (r4on_n < real4_tn or r4on_n <= expected4_tn):
        status4_n = 'PASS'
      else:
        status4_n = 'FAIL'

    lower = (r1on_n, r2on_n, r3on_n, r4on_n)
    upper = (r1off_n, r2off_n, r3off_n, r4off_n)
    tsensors_on = (min(upper) - max(lower))
    tsensors_on_percent = (tsensors_on / TON) * 100

    # print 'TON = %.3f, TOFF= %.3f, TCYCLE = %.3f' % (round(TON, 3), round(TOFF, 3), round(TCYCLE, 3))
    print 'r1on_n\t%.3f\tr1off_n\t%.3f' % (round(r1on_n, 3), round(r1off_n, 3))
    print 'r2on_n\t%.3f\tr2off_n\t%.3f' % (round(r2on_n, 3), round(r2off_n, 3))
    print 'r3on_n\t%.3f\tr3off_n\t%.3f' % (round(r3on_n, 3), round(r3off_n, 3))
    print 'r4on_n\t%.3f\tr4off_n\t%.3f' % (round(r4on_n, 3), round(r4off_n, 3))
    print
    # print round(max(lower), 3), round(min(upper), 3), round(tsensors_on, 3), round(tsensors_on_percent, 3)

    if (n >= int(MAXCICLES * pDiscard)):
      n_counter += 1

      if ((min(upper) - max(lower)) >= delta_success):
        hit += 1

      success = (hit / (n_counter - 1)) * 100  # in percentage

      # Node 1
      a = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\n'.format(str(n),
                                                                                                    str(n1),
                                                                                                    str(
                                                                                                      round(real1_tn,
                                                                                                            3)),
                                                                                                    str(
                                                                                                      round(
                                                                                                        expected1_tn,
                                                                                                        3)),
                                                                                                    str(
                                                                                                      round(delay1, 3)),
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

      a = a.expandtabs(8)
      file_all.write(a)

      b1 = '{0}\n'.format(str(round(delay1, 3)))
      b1 = b1.expandtabs(8)
      file_delay1.write(b1)

      c1 = '{0}\n'.format(str(round(sleepOffset1_n, 3)))
      c1 = c1.expandtabs(8)
      file_bdk_node1.write(c1)

      # Node 2
      a = '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n'.format(str(n2),
                                                                                                str(
                                                                                                  round(real2_tn, 3)),
                                                                                                str(
                                                                                                  round(expected2_tn,
                                                                                                        3)),
                                                                                                str(round(delay2, 3)),
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

      a = a.expandtabs(8)
      file_all.write(a)

      b2 = '{0}\n'.format(str(round(delay2, 3)))
      b2 = b2.expandtabs(8)
      file_delay2.write(b2)

      d2 = '{0}\n'.format(str(round(sleepOffset2_n, 3)))
      d2 = d2.expandtabs(8)
      file_bdk_node2.write(d2)

      # Node 3
      a = '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n'.format(str(n3),
                                                                                                str(round(real3_tn, 3)),
                                                                                                str(round(expected3_tn,
                                                                                                          3)),
                                                                                                str(round(delay3, 3)),
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

      a = a.expandtabs(8)
      file_all.write(a)

      b3 = '{0}\n'.format(str(round(delay3, 3)))
      b3 = b3.expandtabs(8)
      file_delay3.write(b3)

      d3 = '{0}\n'.format(str(round(sleepOffset3_n, 3)))
      d3 = d3.expandtabs(8)
      file_bdk_node3.write(d3)

      # Node 4
      a = '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\n\n'.format(str(n4),
                                                                                                              str(
                                                                                                                round(
                                                                                                                  real4_tn,
                                                                                                                  3)),
                                                                                                              str(
                                                                                                                round(
                                                                                                                  expected4_tn,
                                                                                                                  3)),
                                                                                                              str(
                                                                                                                round(
                                                                                                                  delay4,
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

      e4 = '{0}\n'.format(str(round(sleepOffset4_n, 3)))
      e4 = e4.expandtabs(8)
      file_bdk_node4.write(e4)

      a = '{0}\n'.format(str(round(TCYCLE, 3)))
      a = a.expandtabs(8)
      file_tcycle.write(a)

    # ***************************************************************************
    # generate delays for nodes for next cycle
    # use the folowing two lines if you want to reproduce the experiment and get the same results
    # np.random.seed(j + n)
    # np.random.RandomState(j + n)

    if (tDelayDist == 'uniform'):
      # node 1
      rnd_tDelayDist = 'np.random.' + tDelayDist + '((1 - sigma) * delay_1, (1 + sigma) * delay_1)'
      delay1 = eval(rnd_tDelayDist)

      # node 2
      rnd_tDelayDist = 'np.random.' + tDelayDist + '((1 - sigma) * delay_2, (1 + sigma) * delay_2)'
      delay2 = eval(rnd_tDelayDist)

      # node 3
      rnd_tDelayDist = 'np.random.' + tDelayDist + '((1 - sigma) * delay_3, (1 + sigma) * delay_3)'
      delay3 = eval(rnd_tDelayDist)

      # node 4
      rnd_tDelayDist = 'np.random.' + tDelayDist + '((1 - sigma) * delay_4, (1 + sigma) * delay_4)'
      delay4 = eval(rnd_tDelayDist)

    if (tDelayDist == 'normal'):
      # node 1
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_1, sigma * delay_1)'
      delay1 = eval(rnd_tDelayDist)

      # node 2
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_2, sigma * delay_2)'
      delay2 = eval(rnd_tDelayDist)

      # node 3
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_3, sigma * delay_3)'
      delay3 = eval(rnd_tDelayDist)

      # node 4
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_4, sigma * delay_4)'
      delay4 = eval(rnd_tDelayDist)

    if (tDelayDist == 'exponential'):
      # node 1
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_1)'
      delay1 = eval(rnd_tDelayDist)

      # node 2
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_2)'
      delay1 = eval(rnd_tDelayDist)

      # node 3
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_3)'
      delay3 = eval(rnd_tDelayDist)

      # node 4
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_4)'
      delay4 = eval(rnd_tDelayDist)

    if (tDelayDist == 'chisquare'):
      # node 1
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_1)'
      delay1 = eval(rnd_tDelayDist)

      # node 2
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_2)'
      delay2 = eval(rnd_tDelayDist)

      # node 3
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_3)'
      delay3 = eval(rnd_tDelayDist)

      # node 4
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_4)'
      delay4 = eval(rnd_tDelayDist)

    if (tDelayDist == 'poisson'):
      # node 1
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_1)'
      delay1 = eval(rnd_tDelayDist)

      # node 2
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_2)'
      delay2 = eval(rnd_tDelayDist)

      # node 3
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_3)'
      delay3 = eval(rnd_tDelayDist)

      # node 4
      rnd_tDelayDist = 'np.random.' + tDelayDist + '(delay_4)'
      delay4 = eval(rnd_tDelayDist)

    if (tDelayDist == 'constant'):
      # nodes 1, 2, 3, 4
      delay1 = delay_1
      delay2 = delay_2
      delay3 = delay_3
      delay4 = delay_4

    # ***************************************************************************
    # generate random TCycle next cycle
    # use the folowing two lines if you want to reproduce the experiment and get the same results
    # np.random.seed(j + n)
    # np.random.RandomState(j + n)

    if (tCycleDist == 'uniform'):
      rnd_tCycleDist = 'np.random.' + tCycleDist + '(TMINCYCLE, TMAXCYCLE)'
      TCYCLE = eval(rnd_tCycleDist)

    if (tCycleDist == 'normal'):
      rnd_tCycleDist = 'np.random.' + tCycleDist + '(TMINCYCLE, sigma * TMAXCYCLE)'
      TCYCLE = eval(rnd_tCycleDist)

    if (tCycleDist == 'exponential'):
      rnd_tCycleDist = 'np.random.' + tCycleDist + '(TMINCYCLE)'
      TCYCLE = eval(rnd_tCycleDist)

    if (tCycleDist == 'chisquare'):
      rnd_tCycleDist = 'np.random.' + tCycleDist + '(TMINCYCLE)'
      TCYCLE = eval(rnd_tCycleDist)

    if (tCycleDist == 'poisson'):
      rnd_tCycleDist = 'np.random.' + tCycleDist + '(TMINCYCLE)'
      TCYCLE = eval(rnd_tCycleDist)

    if (tCycleDist == 'constant'):
      TCYCLE = TMINCYCLE

    TOFF = TCYCLE - TON

  a = '{0}\t{1}\n'.format('Total hit success', str(round(success, 3)))
  a = a.expandtabs(8)
  file_all.write(a)

  a = '{0}\n'.format(str(round(success, 3)))
  a = a.expandtabs(8)
  file_tSensorsOn_success.write(a)

  temp_success += success
  print 'Simulation %d: hit = %d, n. cycles = %d, success = %.3f' % (j, hit, (n_counter - 1), success)

  file_n.close()
  file_all.close()
  file_delay1.close()
  file_delay2.close()
  file_delay3.close()
  file_tSensorsOn.close()
  file_tSensorsOn_percent.close()
  # file_bdk.close()
  file_bdk_node1.close()
  file_bdk_node2.close()
  file_bdk_node3.close()
  file_bdk_node4.close()

mean_success = temp_success / NSIM
print '\nMean Success: %.2f\n' % (mean_success)

file_tSensorsOn_success.close()

stop_t = time.clock()
print 'Stop processing at: %.3f' % (stop_t)
print 'Total processing time: %.3f\n' % (stop_t - start_t)
# ***************************************************************************
