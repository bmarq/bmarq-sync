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


#
# App to run several times using multiple processor cores the script to evaluate the bmarq-sync mechanism
#

from __future__ import division
import argparse
import subprocess
import time
import os
import sys
from pylab import *

global alpha, beta, gamma, counter

__version__ = 'multiRun.py v1.0, (C)2017, INESC TEC, IPV/ESTGV, ' + \
              'Bruno Marques bmarq@estgv.ipv.pt'

# ***************************************************************************
parser = argparse.ArgumentParser(
  description='Run several times the script bmarq-sync.py to evaluate the synchronization mechanism using all possible combinations of the mechanism parameters. Type <bmarq-sync-eval.py --help> to know all parameters available to retouch me!')
parser.add_argument('--progname', default='./bmarq-sync-eval.py',
                    help='Name of python script (default is <bmarq-sync-eval.py>)')
parser.add_argument('--runs', type=int, default=25, dest='runs', metavar="N",
                    help='Run <PROGNAME>  RUNS times (default is 25)')
parser.add_argument('--cycles', type=long, default=10E5,
                    help='Maximum number of cycles per simulation (default is 10E5)')
parser.add_argument('--version', action='version', version=__version__, help='Version number of multiRun.py')

args = parser.parse_args()

# ***************************************************************************
alpha = sorted({0.125, 0.5, 0.875})
beta = sorted({1, 2, 3, 4, 5})
gamma = sorted({0.5, 0.6, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95})
cycledist = sorted({'constant', 'uniform', 'normal', 'exponential', 'chisquare', 'poisson'})
delaydist = sorted({'constant', 'uniform', 'normal', 'exponential', 'chisquare', 'poisson'})
counter = 0

# ***************************************************************************
# Main block of the program
subprocess.call('clear', shell=True)  # clearing stdio
start_t = time.clock()

# ***************************************************************************
print 'Started processing at: %f ...\n' % (start_t)

if not os.path.isfile(args.progname):
  print "File not found:", args.progname
  sys.exit(1)

if not os.path.exists('results'):
  os.makedirs('results')

if not os.path.exists('graphs'):
  os.makedirs('graphs')

for t in cycledist:
  for d in delaydist:
    for g in gamma:
      for a in alpha:
        for b in beta:
          cmd = args.progname + ' --numsim=' + str(args.runs) + ' --maxcycles=' + str(args.cycles) + ' --cycledist=' + str(t) + ' --delaydist=' + str(d) + ' --gamma=' + str(g) + ' --alpha=' + str(a) + ' --beta=' + str(b)
          #print cmd
          os.system(cmd)
          counter += 1

stop_t = time.clock()
print 'Stop processing at: %.3f' % (stop_t)
print 'Total processing time: %.3f\n' % (stop_t - start_t)
print 'total of possible combinations: %d' % (counter)
# ***************************************************************************
