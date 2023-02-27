#!/usr/bin/env python
#
# compare wave function 
# usage: comp_wave.py  foo.wav bar.wav

import sys, struct

wavfile1 = open(sys.argv[1],'rb')
neig = struct.unpack('<i',wavfile1.read(4))[0]
wavfile1.seek(4+12*neig,1)

wavfile2 = open(sys.argv[2],'rb')
wavfile2.seek(8+12*neig,1)

xl = 0.0
while True:
    x = wavfile1.read(8)
    y = wavfile2.read(8)
    if (not x) and (not y): break
    if (not x) or (not y): 
        print '\n different file size \n'
        sys.exit(1)
    x = struct.unpack('<d',x)[0]
    y = struct.unpack('<d',y)[0]
    xy = abs(x-y)
    xl = max(xl, xy)
    if xy > 3e-5:
        print "%25s" % repr(x), "%25s" % repr(y), "%25s" % repr(x-y)

print
print 'finish compare.    max. diff : %25s' % repr(xl)
print
