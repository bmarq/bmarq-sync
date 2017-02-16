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
# Run all scripts: simulation and plots from results extracted from the bmarq-sync-eval.py script
#

from __future__ import division
import argparse
import subprocess
import time
import os


__version__ = 'runAll.py v1.0, (C)2017, INESC TEC, IPV/ESTGV, ' + \
              'Bruno Marques bmarq@estgv.ipv.pt'

# ***************************************************************************
parser = argparse.ArgumentParser(
  description='Run all scripts: simulation and plots of results extracted from the bmarq-sync-eval.py script')
parser.add_argument('--version', action='version', version=__version__, help='Version number of runAll.py')

args = parser.parse_args()

# ***************************************************************************

# ***************************************************************************
# Main block of the program
subprocess.call('clear', shell=True)  # clearing stdio
start_t = time.clock()

# ***************************************************************************
print 'Started processing at: %f ...\n' % (start_t)

cmd = './multiRun.py --cycles=100'
os.system(cmd)

cmd = './scripts/multiPlot.py'
os.system(cmd)

stop_t = time.clock()
print 'Process startde at: %f ...\n' % (start_t)
print 'End at: %.3f' % (stop_t)
print 'Total processing time: %.3f\n' % (stop_t - start_t)
# ***************************************************************************
