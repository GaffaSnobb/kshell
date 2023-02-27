#!/usr/bin/env python
#
# covert binary format of the KSHELL wave function
# ( binary format is autodetected )
#
# ./byteswap_wav.py input.wav output.wav
#

kwf = 8   # wave function in double precision
# kwf = 4  # wave function in single precision

# nchunk = 65536
nchunk = 2**20


#-----------------------------

import struct, sys, os
import array
import time

from collect_logs import str_JJ


dt_kwf = 'd' if kwf==8 else 'f'

ei, eo = '>', '<' # try big endian (sparc) to little endian (intel)

fn_in, fn_ot = sys.argv[1], sys.argv[2]

if os.path.exists(fn_ot):
    sys.exit('\n ERROR: output file exist '+fn_ot+'\n')

fp_in = open(fn_in, 'rb')
fp_ot = open(fn_ot, 'wb')

# header 
neig, mtotal = struct.unpack( ei + 'ii', fp_in.read(8) )
if neig >= 2**24:
    neig, mtotal = struct.unpack(ei+'ii', struct.pack(eo+'ii', neig, mtotal))
    ei, eo = eo, ei

e_val = struct.unpack( ei + '%id'%neig, fp_in.read(8*neig) )
jj    = struct.unpack( ei + '%ii'%neig, fp_in.read(4*neig) )

if ei == '>':
    print '\n *** convert big endian (e.g. SPARC) to little endian (e.g. x86) *** \n'
else:
    print '\n *** convert little endian (e.g. x86) to big endian (e.g. SPARC) *** \n'

print '  neig,     J,     Energy        M= %s' % (str_JJ(mtotal),)
for i in range(neig):
    print ' %5i %6s   %10.3f' % (i+1, str_JJ(jj[i]), e_val[i])
print
    
fp_ot.write( struct.pack( eo + 'ii', neig, mtotal ) )
fp_ot.write( struct.pack( eo + '%id'%neig, *e_val ) )
fp_ot.write( struct.pack( eo + '%ii'%neig, *jj) )


# w.f.
dim = ( os.path.getsize( fn_in ) - 8 - 8*neig - 4*neig ) / kwf
d_mb = dim*kwf/1024./1024.

start = time.time()

for i in range((dim-1)/nchunk+1):
    nb = min(nchunk, dim-i*nchunk)
    v = array.array(dt_kwf)
    v.fromfile(fp_in, nb)
    v.byteswap()
    # import numpy
    # v = numpy.fromfile( fp_in, dtype = dt_kwf, count = nb )
    # v = v.byteswap()
    v.tofile( fp_ot )

    d = (i*nchunk + nb)*kwf/1024./1024.
    e_time = time.time() - start
    est_time = e_time * (d_mb - d) / d
    sys.stdout.write( '\r %9.1f MB /%9.1f MB:  %5.1f ' % (d, d_mb, d/d_mb*100.)
                      + '%, ' +
                      'Elapsed %9.1f sec, Remain %9.1f sec ' % ( e_time, est_time )
    )
print 
print 


# v = struct.unpack( ei + str(dim)+dt_kwf, fp_in.read(kwf*dim) )
# fp_ot.write( struct.pack( eo + str(dim)+dt_kwf, *v ) )

fp_in.close()
fp_ot.close()
