#!/usr/bin/python3
#
# draw E2/M1 map by reading summary_foobar
#
# ./map_transit.py summary_foobar E2
# ./map_transit.py summary_foobar M1 
#
#
#

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

thd = 0.    # threashold to show for E2
#thd = 100.

# scaling of width for E2/M1
w_scale = 1./400.  # 1 # 2

shft_np = 0.03 # shift J value for negative/positive parity

import sys
from math import *
import pylab 
import numpy as np


def read_summary_file(fn):
    
    fp = open(fn, 'r')
    for i in range(5):
        fp.readline()

    levels = []
    trns   = []
    JPN_to_n = dict()
    while True:
        arr = fp.readline().split()
        if not arr: break
        isshow = True
        if arr[0] == '#':
            isshow = False 
            del arr[0]
        N, p, N_Jp = int(arr[0]), arr[2], int(arr[3])

        if len(levels) != N-1: raise 'error number list'
        
        J = float(arr[1][:-2])*0.5 if '/' in arr[1] else int(arr[1])
        Ex = float(arr[6])
        
        levels.append( (J, p, N_Jp, Ex, isshow) )
        JPN_to_n[ (J, p, N_Jp) ] = N - 1

    while True:
        line = fp.readline()
        if not line: break
        if len(line)<6 or not line[:2] == 'B(':
            continue

        mpl = line[2:4]

        for i in range(2): line = fp.readline()
        
        while True:
            line = fp.readline()
            if len(line)<=1: break

            l = line[:4]
            J1 = float(l[:-2])*0.5 if '/' in l else int(l)
            N1 = int(line[6:8])
            p1 = line[4]
                       
            l = line[17:21]
            J2 = float(l[:-2])*0.5 if '/' in l else int(l)
            N2 = int(line[23:25])
            p2 = line[21]

            ni = JPN_to_n[ (J1, p1, N1) ]
            nj = JPN_to_n[ (J2, p2, N2) ]

            v  = float(line[41:49])
            trns.append( (ni, nj, mpl, v)  )

    fp.close()

    return levels, trns



def str_JJ(jj):
    if jj % 1 == 0: 
        return str(int(jj))
    else:
        return str(int(jj*2)) + '/2'


def main(fn, prty, mtpl, levels, trns):

    w_scale_l = w_scale
    if mtpl == 'M1': w_scale_l = 1000 * w_scale 

    i_prty = { '+': 0, '-' : 1}

    prty_list = ('+', '-')
    if prty == '-':
        prty_list = ('-',)
    elif prty == '+':
        prty_list = ('+',)

    
    x, y = [[], []], [[], []]
    for J, p, N_Jp, Ex, isshow in levels:
        # if p != prty: continue
        if not isshow: continue
        i = i_prty[p]
        x[i].append( J )
        y[i].append( Ex )


    pylab.rcParams['font.family'] = 'Times New Roman'


    pylab.title( fn + ' : ' + mtpl )
    pylab.xlim( [min(x[0]+x[1])-0.5, max(x[0]+x[1])+0.5]  )
    pylab.ylim( [-0.1, max(y[0]+y[1])*1.05] )

    xlist = sorted(list(set( x[0] + x[1] )))
    pylab.xticks(xlist, [str_JJ(j) for j in xlist ])

    # def plot_trn(dJ, prty, cs, alpha, dashes):
    def plot_trn(dJ, prty, cs, args):

        jjee = {}
        
        for ni, nj, mpl, v in trns:
            if mpl != mtpl: continue
            if v < thd: continue
            # print '# transition line : ', ni, nj, mpl, v
            J1, p1, Njp1, Ex1, is1 = levels[ni]
            J2, p2, Njp2, Ex2, is2  = levels[nj]
            if not is1 or not is2: continue
            if p1 != prty or p2 != prty: continue
            if abs(J1-J2) != dJ: continue

            lw = v * w_scale_l
            jjee[lw] = (J1, J2), (Ex1, Ex2)

        for lw in sorted(jjee.keys(), reverse=True):
            (J1, J2), (Ex1, Ex2) = jjee[lw]

            shft = shft_np if prty=='+' else -shft_np
            J1, J2 = J1 + shft, J2 + shft

            pylab.plot( [J1, J2], [Ex1, Ex2], cs, linewidth=lw,
                        **args)
            if 'label' in args: del args['label']


    def blabel(dJ, prty):
        return r'B(%s ; $J^{%s}\pm %s \rightarrow J^{%s}$)' % (mtpl, prty, dJ, prty)

    # line properties, color, dash, label, alpha (transparency), ... 
    if '+' in prty_list:
        # positive parity, dJ=2
        args = {'alpha':1., 'dashes': (None, None), 'label': blabel('2', '+') }
        plot_trn(2, '+', 'r-', args)

        # positive parity, dJ=1
        args = {'alpha':1., 'dashes': (8,2), 'label': blabel('1', '+') }
        plot_trn(1, '+', 'r-', args)

        # positive parity, dJ=0
        args = {'alpha':1., 'dashes': (2,1), 'label': blabel('0', '+') }
        plot_trn(0, '+', 'r-', args)
        
    if '-' in prty_list:
        # negative parity, dJ=2
        args = {'alpha':1., 'dashes': (None,None), 'label': blabel('2', '-') }
        plot_trn(2, '-', 'c-', args)

        # negative parity, dJ=1
        args = {'alpha':1., 'dashes': (8,2), 'label': blabel('1', '-') }
        plot_trn(1, '-', 'c-', args)

        # negative parity, dJ=1
        args = {'alpha':1., 'dashes': (4,1), 'label': blabel('0', '-') }
        plot_trn(0, '-', 'c-', args)
    
    if '+' in prty_list:
        xs = [ a+shft_np for a in x[0] ]
        # pylab.plot(x[0], y[0], 'ko', markerfacecolor='w', label=r'$J^+$ state')
        pylab.plot(xs, y[0], 'ko', markerfacecolor='w', label=r'$J^+$ state')
    if '-' in prty_list:
        xs = [ a-shft_np for a in x[1] ]
        # pylab.plot(x[1], y[1], 'b^', alpha=0.7, label=r'$J^-$ state')
        pylab.plot(xs, y[1], 'b^', alpha=0.7, label=r'$J^-$ state')

    lab = r'$J^+/J^-$'
    if len(prty_list)==1:
        lab = r'$J^{%s}$' % (prty_list[0])
    pylab.xlabel(lab,   fontsize=20)
    pylab.ylabel('$E_x$ (MeV)', fontsize=20)
    

    pylab.legend(loc='lower right')

    plot_exp(prty)

    pylab.tight_layout()

#    pylab.savefig('E2map.eps')

    pylab.show()


def plot_exp(prty):
    fn = './exp.dat'
    import os
    if not os.path.exists(fn):
        return
    fp = open(fn, 'r')
    x, y = [], []
    for line in fp.readlines():
        if not line: break
        if line[0] == "#": continue
        arr = line.split()
        J, p, n, Ex = int(arr[0]), arr[1], int(arr[2]), float(arr[3])
        if prty != p: continue
        x.append( J+0.05 )
        y.append( Ex )
    
    # pylab.plot(x, y, 'bx', markersize=10, linewidth=10)
    pylab.scatter(x, y, s=40, c='b', marker='x', linewidth=2)
    
    
        
        


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print ("usage : map_transit.py summary_foobar +(parity)  M1/E2")
        sys.exit(1)
    
    fname = sys.argv[1]
    prty = None
    mtpl = 'E2'

    for asc in sys.argv[2:]:
        if asc in ('+',  '-' ): prty = asc
        if asc in ('M1', 'm1'): mtpl = 'M1'
        if asc in ('E1', 'e1'): mtpl = 'E1'

    levels, trns = read_summary_file(fname)
    main(fname, prty, mtpl, levels, trns)

