#!/usr/bin/env python
#
# usage: extract_wf.py foo.snt bar.ptn input.wav output.wav n
#
#  extract n-th  wave function
#

import sys, struct, itertools
from math import *
import numpy

from count_dim import read_snt, read_ptn, \
    mp_add, mp_product, mps_product
from ratio_configuration_ptn import kwf, edn 

dt_wf = 'd' if kwf==8 else 'f'

chunksize = 10**7

if __name__ == "__main__":

    if len(sys.argv) < 6:
        print( '\nusage: extract_wf.py foo.snt bar.ptn input.wav output.wav n-m k l ...' )
        print(   '       n-m k l ... n-th state to copy\n' )
        sys.exit(1)

    fn_snt, fn_ptn, fn_wav, fn_out = sys.argv[1:5]
    outputlist = []
    for n in sys.argv[5:]:
        if '-' in n:
            n = n.split('-')
            outputlist += list( range(int(n[0])-1, int(n[1])) )
        elif ',' in n:
            outputlist += [ int(i)-1 for i in n.split(',')]
        else:
            outputlist += [ int(n)-1, ]
    
    
    n_jorb, n_core, norb, lorb, jorb, itorb = read_snt( fn_snt )
    n_ferm, iprty, ptn_p, ptn_n, ptn_pn = read_ptn( fn_ptn ) 

    print( "snt file: " + fn_snt )
    print( "ptn file: " + fn_ptn )
    print( "wav file: " + fn_wav )
    print( "Z=%3i, N=%3i" % (n_ferm[0]+n_core[0],n_ferm[1]+n_core[1]) )
    print( "parity=%2i" % iprty )

    # set M-scheme bits for single-j orbits and J^+ matrix element
    dim_jnm, mbit_jn, m_j_mbit, jup_j = {}, {}, {}, {}
    for j in set( jorb ):
        dim_jnm[j]  = [ {} for x in range(j+2) ]
        mbit_jn[j]  = [ [] for x in range(j+2) ]
        m_j_mbit[j] = []
        jup_j[j]    = []
        for mbit in range(2**(j+1)):
            n, m = 0, 0
            for i in range(j+1):
                n += (mbit&2**i)>>i
                m += ((mbit&2**i)>>i)*(2*i-j)
            dim_jnm[j][n][m] = dim_jnm[j][n].get(m, 0) + 1
            mbit_jn[j][n].append( mbit )
            m_j_mbit[j].append( m )
            jup_j[j].append( {} )
            for i in range(j):
                if (mbit&(3*2**i))>>i == 1:
                    jup_j[j][mbit][mbit+2**i] = sqrt( (j-i)*(i+1) )

    
    # dimension for proton
    dim_idp_mp = []
    for ptn in ptn_p:
        mps = []
        for i,n in enumerate(ptn):
            p = (-1)**(lorb[i]*n) 
            j = jorb[i]
            mps.append( dict( ( (m, p), d ) for m,d in list(dim_jnm[j][n].items()) ) )
        dim_idp_mp.append(mps_product(mps))

    # dimension for neutron
    dim_idn_mp = []
    for ptn in ptn_n:
        mps = []
        for i,n in enumerate(ptn):
            p = (-1)**( lorb[ n_jorb[0]+i ] * n )
            j = jorb[ n_jorb[0]+i ]
            mps.append( dict( ( (m, p), d ) for m,d in list(dim_jnm[j][n].items()) ) )
        dim_idn_mp.append( mps_product(mps) )

    # product dimensions of proton and neutron
    dim_mp = {}
    for idp,idn in ptn_pn:
        mp_add( dim_mp, mp_product(dim_idp_mp[idp], dim_idn_mp[idn]) )

    fp_wav = open(fn_wav,'rb')
    neig, mtotal = numpy.fromfile( fp_wav, dtype=edn+'i', count=2 )
    e_val        = numpy.fromfile( fp_wav, dtype=edn+'d', count=neig )
    jj           = numpy.fromfile( fp_wav, dtype=edn+'i', count=neig )
    print( "neig=%3i, M=%3i/2" % (neig,mtotal) )
    print( "  i     JP       E(MeV)" )
    for i in range(neig):
        print( "%3i %3i/2%s %12.5f"%(i+1,jj[i],"- +"[iprty+1], e_val[i]) )
    print()

    ndim = dim_mp[mtotal,iprty]
    print( "M-scheme dim.: ", ndim )
    print()
    

    neig2 = len(outputlist)

    if neig2 == 0: 
        print( "\n*** NO wave function ***\n" )
        sys.exit(0)

    fp_out = open(fn_out,'wb')

    print( "wav output file: " + fn_out )
    print( "neig=%3i, M=%3i/2" % (neig2,mtotal) )
    print( "  i     JP       E(MeV)" )


    
    for j,i in enumerate(outputlist):
        print( "%3i %3i/2%s %12.5f" % (j+1, jj[i], "- +"[iprty+1], e_val[i]) )
    print()

    fp_out.write( struct.pack( edn + 'ii', neig2, mtotal ) )
    fp_out.write( struct.pack( edn + '%id'%neig2, *[e_val[x] for x in outputlist] ) )
    fp_out.write( struct.pack( edn + '%ii'%neig2, *[jj[x] for x in outputlist] ) )

    for istate in range(neig):
        if not istate in outputlist:
            print( istate+1, "/", neig, " skip " )
            fp_wav.seek( kwf*ndim, 1 )
            continue

        print( istate+1, "/", neig, " copy " )
        
        for i in range(0, ndim, chunksize):
            ns = min( ndim - i, chunksize )
            wfdata  = numpy.fromfile( fp_wav, dtype=edn+dt_wf, count=ns )
            fp_out.write( struct.pack( edn + '%i'%(ns)+dt_wf, *wfdata ) )

    fp_wav.close()
    fp_out.close()
 
