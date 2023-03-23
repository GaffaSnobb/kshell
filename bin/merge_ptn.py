#!/usr/bin/env python
#
# usage: merge_ptn.py foo1.ptn foo2.ptn ...
#
# Thanks to Y. Tsunoda
#

import sys

ptn_pn = set()
for fn in sys.argv[1:]:
    inputfile = open(fn,'r')
    fn_snt = ""
    arr = inputfile.readline().split()
    while arr[0][0] in ['!','#']:
        if len(arr)>=5 and arr[1:4]==["partition","file","of"]:
            fn_snt = arr[4]
        arr = inputfile.readline().split()
    n_ferm = [int(arr[0]),int(arr[1])]
    iprty = int(arr[2])
    arr = inputfile.readline().split()
    while arr[0][0] in ['!','#']:
        arr=inputfile.readline().split()
    n_idp = int(arr[0])
    n_idn = int(arr[1])
    arr = inputfile.readline().split()
    while arr[0][0] in ['!','#']:
        arr=inputfile.readline().split()

    ptn_p = {}
    while True:
        if int(arr[0]) in ptn_p:
            print >> sys.stderr, "duplicated number for proton partition"
            exit()
        ptn_p[ int(arr[0]) ] = tuple([int(x) for x in arr[1:]])
        arr = inputfile.readline().split()
        if arr[0][0] in ['!', '#']: break
    while arr[0][0] in ['!', '#']:
        arr=inputfile.readline().split()

    ptn_n={}
    while True:
        if int(arr[0]) in ptn_n:
            print >> sys.stderr, "duplicated number for neutron partition"
            exit()
        ptn_n[ int(arr[0]) ] = tuple([int(x) for x in arr[1:]])
        arr = inputfile.readline().split()
        if arr[0][0] in ['!','#']: break
    while arr[0][0] in ['!','#']:
        arr = inputfile.readline().split()

    n_idpn = int(arr[0])
    while True:
        arr = inputfile.readline().split()
        if not arr: break
        if int(arr[0]) in ptn_p and int(arr[1]) in ptn_n:
            ptn_pn.add((ptn_p[int(arr[0])], ptn_n[int(arr[1])]))
    inputfile.close()

ptn_p2 = sorted(list(set([x[0] for x in ptn_pn])))
ptn_n2 = sorted(list(set([x[1] for x in ptn_pn])))
ptn_pn2 = [ (ptn_p2.index(x), ptn_n2.index(y)) \
            for x,y in sorted(list(ptn_pn)) ]

print "# partition file of %s  Z=%d  N=%d  parity=%+d" \
    % (fn_snt, n_ferm[0], n_ferm[1], iprty)
print " %d %d %d" % (n_ferm[0], n_ferm[1], iprty)
print "# num. of  proton partition, neutron partition"
print " %d %d" % (len(ptn_p2),len(ptn_n2))
print "# proton partition"
for i, arr in enumerate(ptn_p2):
    print " %5d   " % (i+1),
    for a in arr:
        print"%2d" % a,
    print
print "# neutron partition"
for i, arr in enumerate(ptn_n2):
    print " %5d   " % (i+1),
    for a in arr:
        print"%2d" % a,
    print
print "# partition of proton and neutron"
print "%d" % len(ptn_pn2)
for i,j in ptn_pn2:
    print "%5d %5d" % (i+1, j+1)
