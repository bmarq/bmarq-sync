#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import subprocess
import time
from pylab import *
import numpy as np
import scipy as sp
from pylab import *

global TON, TOFF, TCYCLE, TMAX
global t
global p_discard
global n1, n2, n3
global delay_1, delay_2, delay_3
global delay1, delay2, delay3
global expected1_tk, expected2_tk, expected3_tk
global sampled1_tk, sampled2_tk, sampled3_tk
global sampled1_tk_1, sampled2_tk_1, sampled3_tk_1
global delta1_k, delta2_k, delta3_k
global delta1_k_1, delta2_k_1, delta3_k_1
global tsleep1_k, tsleep2_k, tsleep3_k
global r1on, r2on, r3on
global r1off, r3off, r3off
global alpha, beta
global i_counter, hit
global lower, upper
global success, delta_success
global ton_width, ton_width_percent

# ***************************************************************************

'''---8<------8<------8<------8<------8<------8<------8<------8<------8<--- '''
# Main block of the program

subprocess.call('clear', shell=True)  # clearing stdio
start_t = time.clock()

TMAX = 10E5  # 10E5
TON = 60
TOFF = 900
TCYCLE = TON + TOFF

t = 0.0
p_discard = 0.10  # ignore 1st 10% processing time

delay_1 = 1/2
delay_2 = 1
delay_3 = 2

expected1_tk = expected2_tk = expected3_tk = 0.0  # t'k for nodes 1, 2 and 3 (expected reception time for packet k)
sampled1_tk = sampled2_tk = sampled3_tk = 0.0  # tk -1 for nodes 1, 2 and 3 (real reception time for packet k-1 [actual (k)])
sampled1_tk_1 = sampled2_tk_1 = sampled3_tk_1 = 0.0  # tk -1 for nodes 1, 2 and 3 (real reception time for packet k-1 [anterior (k-1)])

delta1_k = delta2_k = delta3_k = 0.0
delta1_k_1 = delta2_k_1 = delta3_k_1 = 0.0

tsleep1_k = tsleep2_k = tsleep3_k = 0.0

r1on = r2on = r3on = 0.0
r1off = r2off = r3off = TON

n1=1    # Node 1
n2=2    # Node 2
n3=3    # Node 3

alpha = 0.875  # test with 0.125; 0.50; 0.875
beta = 100  # test with 1, 10, 50, 100

i_counter = 0
hit = 0
success = 0
delta_success = TON * 0.8 # 80% of TON

file_i = open('data-uniform_delay_i.txt', 'w')
file_all = open('result-uniform_delay_deviation.txt', 'w')
file_delay1 = open('data-uniform_delay-1.txt', 'w')
file_delay2 = open('data-uniform_delay-2.txt', 'w')
file_delay3 = open('data-uniform_delay-3.txt', 'w')
file_ton_width = open('data-uniform_ton_width.txt', 'w')
file_ton_width_percent = open('data-uniform_ton_width_percent.txt', 'w')
file_bdk = open('data-uniform_bdk.txt', 'w')

a = 'i\tnode\tt\'k\ttk\tdelay\tt\'k-tk\tdk\tr_on\tr_off\ttsleep\tb.|dk|\tton_width\t%\n'
file_all.write(a)

print 'Started processing at: %f \n ...\n' % (start_t)

'''---8<------8<------8<------8<------8<------8<------8<------8<------8<--- '''
for i in range(0, int(TMAX) + 1):
    print i

    # initialization for node n1
    delay1 = np.random.uniform(0.80 * delay_1, 1.20 * delay_1) # +/ 20%  delay 1
    print delay1

    sampled1_tk_1 = sampled1_tk
    sampled1_tk = i * TCYCLE + delay1
    expected1_tk = sampled1_tk_1 + TON + TOFF
    delta1_k_1 = delta1_k

    # initialization for node n2
    delay2 = np.random.uniform(0.80 * delay_2, 1.20 * delay_2) # +/ 20%  delay 2
    print delay2

    sampled2_tk_1 = sampled2_tk
    sampled2_tk = i * TCYCLE + delay2
    expected2_tk = sampled2_tk_1 + TON + TOFF
    delta2_k_1 = delta2_k

    # initialization for node n3
    delay3 = np.random.uniform(0.80 * delay_3, 1.20 * delay_3) # +/ 20%  delay 3
    print delay3

    sampled3_tk_1 = sampled3_tk
    sampled3_tk = i * TCYCLE + delay3
    expected3_tk = sampled3_tk_1 + TON + TOFF
    delta3_k_1 = delta3_k

    if (i == 0):
        delta1_k = delta2_k = delta3_k = 0
        r1on = r2on = r3on = 0
    else:
        delta1_k = (1.0 - alpha) * delta1_k_1 + alpha * (expected1_tk - sampled1_tk)
        r1on = expected1_tk - beta * abs(delta1_k)

        delta2_k = (1.0 - alpha) * delta2_k_1 + alpha * (expected2_tk - sampled2_tk)
        r2on = expected2_tk - beta * abs(delta2_k)

        delta3_k = (1.0 - alpha) * delta3_k_1 + alpha * (expected3_tk - sampled3_tk)
        r3on = expected3_tk - beta * abs(delta3_k)

    tsleep1_k = TOFF - beta * abs(delta1_k)
    tsleep2_k = TOFF - beta * abs(delta2_k)
    tsleep3_k = TOFF - beta * abs(delta3_k)

    r1off = r1on + TON
    r2off = r2on + TON
    r3off = r3on + TON

    lower = (r1on, r2on, r3on)
    upper = (r1off, r2off, r3off)
    ton_width = (min(upper) - max(lower))
    ton_width_percent = ((min(upper) - max(lower)) / TON) * 100

    if (i >= int(TMAX * p_discard)):
        #if (i > 0):
        i_counter += 1
        # Node 1
        a = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n'.format(str(i),
                                                                               str(n1),
                                                                    str(
                                                                        round(expected1_tk, 3)),
                                                                    str(
                                                                        round(sampled1_tk, 3)),
                                                                    str(round(delay1, 3)),
                                                                    str(
                                                                        round(expected1_tk - sampled1_tk, 3)),
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
                                                                              3)))

        file_all.write(a)

        b1 = '{0}\n'.format(str(round(delay1, 3)))
        file_delay1.write(b1)

        # Node 2
        a = '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(str(n2),
                                                                    str(
                                                                        round(expected2_tk, 3)),
                                                                    str(
                                                                        round(sampled2_tk, 3)),
                                                                    str(round(delay2, 3)),
                                                                    str(
                                                                        round(expected2_tk - sampled2_tk, 3)),
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
                                                                                  delta2_k),
                                                                              3)))
        file_all.write(a)

        b2 = '{0}\n'.format(str(round(delay2, 3)))
        file_delay2.write(b2)

        # Node 3
        a = '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n\n'.format(str(n3),
                                                                                       str(
                                                                                           round(expected3_tk, 3)),
                                                                                       str(
                                                                                           round(sampled3_tk, 3)),
                                                                                       str(round(delay3, 3)),
                                                                                       str(
                                                                                           round(
                                                                                               expected3_tk - sampled3_tk,
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
                                                                                                     delta3_k),
                                                                                                 3)),
                                                                                       str(round(ton_width, 3)),
                                                                                       str(round(ton_width_percent, 3)))
        file_all.write(a)

        b3 = '{0}\n'.format(str(round(delay3, 3)))
        file_delay3.write(b3)

        c1 = '{0}\n'.format(str(round(ton_width, 3)))
        file_ton_width.write(c1)

        c2 = '{0}\n'.format(str(round(ton_width_percent, 3)))
        file_ton_width_percent.write(c2)

        d = '{0}\n'.format(str(i))
        file_i.write(d)

        e1 = '{0}\n'.format(str(round(beta * abs(delta1_k), 3)))
        file_bdk.write(e1)

        e2 = '{0}\n'.format(str(round(beta * abs(delta2_k), 3)))
        file_bdk.write(e2)

        e3 = '{0}\n'.format(str(round(beta * abs(delta3_k), 3)))
        file_bdk.write(e3)

        if((min(upper) - max(lower))>= delta_success):
            hit += 1

success = (hit / i_counter) * 100 # in percentage

a = '{0}\t{1}\n'.format('total hit success',str(success))
file_all.write(a)

file_i.close()
file_all.close()
file_delay1.close()
file_delay2.close()
file_delay3.close()
file_ton_width.close()
file_ton_width_percent.close()
file_bdk.close()

print 'total hit success: %d \n' %(success)

stop_t = time.clock()
print 'Stop processing at: %f \n' % (stop_t)
print 'total processing time: %f \n' % (stop_t - start_t)
# ***************************************************************************
