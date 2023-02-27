#!/usr/bin/env python
#
# usage: mup_operate.py foo.snt bar.ptn input.wav output.wav
#
#  output.wav(M+2) = J^+ input.wav(M)
# 
# first ver.  by Y. Tsunoda 2017/04/20 
# second ver. by N. Shimizu 2017/06/22
#

import sys, struct, itertools, multiprocessing
from math import *
import numpy, scipy.sparse

from count_dim import read_snt, read_ptn, \
    mp_add, mp_product, mps_product
from ratio_configuration_ptn import kwf, edn 

dt_wf = 'd' if kwf==8 else 'f'
if kwf != 8: raise "not implemented kwf=4"


def wrapper_init_jup_p(x):
    return init_jup_p( x[0], **x[1] )


def init_jup_p(ptn, lorb, jorb, dim_jnm, mbit_jn, m_j_mbit, n_jorb, jup_j, ptn_p):
    print ptn_p.index(ptn)+1,"/",len(ptn_p)
    
    mps = []
    for i,n in enumerate(ptn):
        p = (-1)**(lorb[i]*n) 
        j = jorb[i]
        mps.append( dict( ( (m, p), d ) for m,d in dim_jnm[j][n].iteritems() ) )
    dim_mp = mps_product(mps)

    dim_m = dict( (m,x) for (m,p),x in dim_mp.iteritems() )
    
    mbit_m, idx_m_mbit, mbit2_m = {}, {}, {}
    for mbits in itertools.product(
            *[ mbit_jn[x][y] for x,y in zip(jorb[:n_jorb[0]], ptn) ][::-1] ):
        m = sum( m_j_mbit[x][y] for x,y 
                 in zip(jorb[:n_jorb[0]], mbits[::-1]) )
        if not m in mbit_m:
            mbit_m[m] = []
            idx_m_mbit[m] = {}
            mbit2_m[m] = []
        mbit_m[m].append( mbits[::-1] )
        idx_m_mbit[m][mbits[::-1]]=len(mbit_m[m])-1
        mbit2 = []
        for i,mbit in enumerate(mbits[::-1]):
            for j in range(jorb[i]+1):
                if mbit&2**j!=0:
                    mbit2.append(sum(jorb[:i])+i+j)
        mbit2_m[m].append(tuple(mbit2))

    jup_m = {}
    mmin = min(dim_m)
    jup_m[mmin-2] = scipy.sparse.csr_matrix( (dim_m[mmin], 0), dtype='d' )
    for m in range(mmin, -mmin+2, 2):
        dim1 = dim_m.get(m,0)
        dim2 = dim_m.get(m+2,0)
        jup_m[m] = scipy.sparse.lil_matrix( (dim2, dim1), dtype='d' )
        for idim1,mbitp1 in enumerate(mbit_m[m]):
            for i in range(n_jorb[0]):
                for mbit2, v in jup_j[jorb[i]][mbitp1[i]].items():
                    mbitp2 = mbitp1[:i] + (mbit2,) + mbitp1[i+1:]
                    jup_m[m][idx_m_mbit[m+2][mbitp2],idim1] = v
        jup_m[m] = scipy.sparse.csr_matrix(jup_m[m])

    return dim_mp, dim_m, mbit_m, idx_m_mbit, jup_m, mbit2_m


def wrapper_init_jup_n(x):
    return init_jup_n( x[0], **x[1] )


def init_jup_n(ptn, lorb, jorb, dim_jnm, mbit_jn, m_j_mbit, n_jorb, jup_j, ptn_n):

    print  ptn_n.index(ptn)+1, "/", len(ptn_n)

    mps = []
    for i,n in enumerate(ptn):
        p = (-1)**( lorb[ n_jorb[0]+i ] * n )
        j = jorb[ n_jorb[0]+i ]
        mps.append( dict( ( (m, p), d ) for m,d in dim_jnm[j][n].iteritems() ) )
    dim_mp = mps_product(mps)

    dim_m = dict( (m,x) for (m,p),x in dim_mp.items() )


    mbit_m, idx_m_mbit, mbit2_m = {}, {}, {}
    for mbits in itertools.product( *[ mbit_jn[x][y] for x,y 
                                       in zip( jorb[n_jorb[0]:], ptn) ][::-1] ):
        m = sum( m_j_mbit[x][y] for x,y in zip( jorb[n_jorb[0]:], mbits[::-1]) )
        if not m in mbit_m:
            mbit_m[m] = []
            idx_m_mbit[m] = {}
            mbit2_m[m] = []
        mbit_m[m].append(mbits[::-1])
        idx_m_mbit[m][mbits[::-1]] = len(mbit_m[m]) - 1
        mbit2 = []
        for i,mbit in enumerate(mbits[::-1]):
            for j in range(jorb[i+n_jorb[0]]+1):
                if mbit&2**j!=0:
                    mbit2.append(sum(jorb[n_jorb[0]:i+n_jorb[0]])+i+j)
        mbit2_m[m].append(tuple(mbit2))

    jup_m = {}
    mmin = min(dim_m)
    jup_m[mmin-2] = scipy.sparse.csr_matrix( (dim_m[mmin], 0), dtype='d' )
    for m in range(mmin,-mmin+2,2):
        dim1 = dim_m.get(m,0)
        dim2 = dim_m.get(m+2,0)
        jup_m[m] = scipy.sparse.lil_matrix( (dim2, dim1), dtype='d')
        for idim1,mbitn1 in enumerate(mbit_m[m]):
            for i in range(n_jorb[1]):
                for mbit2, v in jup_j[jorb[i+n_jorb[0]]][mbitn1[i]].items():
                    mbitn2 = mbitn1[:i] + (mbit2,) + mbitn1[i+1:]
                    jup_m[m][idx_m_mbit[m+2][mbitn2],idim1] = v
        jup_m[m] = scipy.sparse.csr_matrix(jup_m[m])

    return dim_mp, dim_m, mbit_m, idx_m_mbit, jup_m, mbit2_m


def operate_jup(wfdata,idp,idn,mtotal):
    wfdata_up={}
    for mtotp in range(min(dim_idp_m[idp]),max(dim_idp_m[idp])+2,2):
        mtotn = mtotal-mtotp
        dimp = dim_idp_m[idp][mtotp]
        dimn = dim_idn_m[idn].get(mtotn,0)
        if dimn==0: continue
        if dim_idp_m[idp].get(mtotp+2,0)!=0:
            if mtotp+2 in wfdata_up:
                wfdata_up[mtotp+2] += wfdata[mtotp]*jup_idp_m[idp][mtotp].T
            else:
                wfdata_up[mtotp+2] = wfdata[mtotp]*jup_idp_m[idp][mtotp].T
        if dim_idn_m[idn].get(mtotn+2,0) != 0:
            if mtotp in wfdata_up:
                wfdata_up[mtotp] += jup_idn_m[idn][mtotn]*wfdata[mtotp]
            else:
                wfdata_up[mtotp] = jup_idn_m[idn][mtotn]*wfdata[mtotp]
    return wfdata_up


def operate_jdown(wfdata,idp,idn,mtotal):
    wfdata_down={}
    for mtotp in range(min(dim_idp_m[idp]),max(dim_idp_m[idp])+2,2):
        mtotn = mtotal - mtotp
        dimp = dim_idp_m[idp][mtotp]
        dimn = dim_idn_m[idn].get(mtotn,0)
        if dimn == 0: continue
        if dim_idp_m[idp].get(mtotp-2,0) != 0:
            if mtotp-2 in wfdata_down:
                wfdata_down[mtotp-2] += wfdata[mtotp]*jup_idp_m[idp][mtotp-2]
            else:
                wfdata_down[mtotp-2] = wfdata[mtotp]*jup_idp_m[idp][mtotp-2]
        if dim_idn_m[idn].get(mtotn-2,0)!=0:
            if mtotp in wfdata_down:
                wfdata_down[mtotp] += jup_idn_m[idn][mtotn-2].T*wfdata[mtotp]
            else:
                wfdata_down[mtotp] = jup_idn_m[idn][mtotn-2].T*wfdata[mtotp]
    return wfdata_down


if __name__ == "__main__":

    if len(sys.argv) < 5:
        print '\nincrease M of the KSHELL wave function by laddger operator\n'
        print 'usage: mup_operate.py foo.snt bar.ptn input_m0p.wav output_m2p.wav\n'
        sys.exit(1)

    fn_snt, fn_ptn, fn_wav, fn_out = sys.argv[1:5]

    n_jorb,n_core,norb,lorb,jorb,itorb = read_snt( fn_snt )
    n_ferm,iprty,ptn_p,ptn_n,ptn_pn = read_ptn( fn_ptn ) 

    print "snt file: " + fn_snt
    print "ptn file: " + fn_ptn
    print "wav file: " + fn_wav
    print "Z=%3i, N=%3i" % (n_ferm[0]+n_core[0],n_ferm[1]+n_core[1])
    print "parity=%2i" % iprty

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

    
    ######################################
    p = multiprocessing.Pool()

    print "initialize j operator for proton"
    xargs = { "lorb":lorb, "jorb":jorb, "dim_jnm":dim_jnm,
              "mbit_jn":mbit_jn, "m_j_mbit":m_j_mbit, "n_jorb":n_jorb, 
              "jup_j":jup_j, "ptn_p":ptn_p }
    dim_idp_mp, dim_idp_m, mbit_idp_m, idx_idp_m_mbit, jup_idp_m, mbit2_idp_m \
        = zip( *p.map(wrapper_init_jup_p, [ (x, xargs) for x in ptn_p ] ))


    print "initialize j operator for neutron"
    xargs = { "lorb":lorb, "jorb":jorb, "dim_jnm":dim_jnm,
              "mbit_jn":mbit_jn, "m_j_mbit":m_j_mbit, "n_jorb":n_jorb, 
              "jup_j":jup_j, "ptn_n":ptn_n }
    dim_idn_mp, dim_idn_m, mbit_idn_m, idx_idn_m_mbit, jup_idn_m, mbit2_idn_m \
        = zip( *p.map(wrapper_init_jup_n,  [ (x, xargs) for x in ptn_n ] ) )

    p.close()
    ######################################
    
    
    dim_idpn_mp = []
    dim_mp = {}
    for idp,idn in ptn_pn:
        dim_idpn_mp.append( mp_product(dim_idp_mp[idp], dim_idn_mp[idn]) )
        mp_add( dim_mp, dim_idpn_mp[-1] )

    fp_wav = open(fn_wav,'rb')
    neig, mtotal = numpy.fromfile( fp_wav, dtype=edn+'i', count=2 )
    e_val        = numpy.fromfile( fp_wav, dtype=edn+'d', count=neig )
    jj           = numpy.fromfile( fp_wav, dtype=edn+'i', count=neig )
    print "neig=%3i, M=%3i/2" % (neig,mtotal)
    print "  i     JP       E(MeV)"
    for i in range(neig):
        print "%3i %3i/2%s %12.5f"%(i+1,jj[i],"- +"[iprty+1], e_val[i])
    print

    outputlist = [ x for x in range(neig) if not jj[x] in [mtotal, -1] ]
    neig2 = len(outputlist)
    mtotal2 = mtotal + 2

    if neig2 == 0: 
        print "\n*** NO wave function ***\n"
        sys.exit(0)

    fp_out = open(fn_out,'wb')

    print "wav output file: " + fn_out
    print "neig=%3i, M=%3i/2" % (neig2,mtotal2)
    print "  i     JP       E(MeV)"

    for j,i in enumerate(outputlist):
        print "%3i %3i/2%s %12.5f" % (j+1, jj[i], "- +"[iprty+1], e_val[i])

    fp_out.write( struct.pack( edn + 'ii', neig2, mtotal2 ) )
    fp_out.write( struct.pack( edn + '%id'%neig2, *[e_val[x] for x in outputlist] ) )
    fp_out.write( struct.pack( edn + '%ii'%neig2, *[jj[x] for x in outputlist] ) )

    # print "calculate wave functions"
    for istate in range(neig):
        print istate+1, "/", neig
        if not istate in outputlist:
            fp_wav.seek( kwf*dim_mp[mtotal,iprty], 1 )
            continue
        ndim_sum = 0
        for idpn in range(len(ptn_pn)):
            idp,idn = ptn_pn[idpn]
            occs = ptn_p[idp] + ptn_n[idn]
            if not dim_idpn_mp[idpn].has_key( (mtotal, iprty) ): continue
            ndim = dim_idpn_mp[idpn][mtotal, iprty]
            wfdata  = numpy.fromfile( fp_wav, dtype=edn+dt_wf, count=ndim )
            wfdata_mat = {}
            idim = 0
            for mtotp in range(min(dim_idp_m[idp]), max(dim_idp_m[idp])+2, 2):
                mtotn = mtotal - mtotp
                dimp = dim_idp_m[idp][mtotp]
                dimn = dim_idn_m[idn].get(mtotn, 0)
                if dimn==0: continue
                wfdata_mat[mtotp] = numpy.matrix( numpy.array( 
                    wfdata[idim:idim+dimp*dimn]).reshape((dimn,dimp) ) )
                idim += dimp*dimn
            wfdata_up = operate_jup(wfdata_mat, idp, idn, mtotal)
            for mtotp in range(min(dim_idp_m[idp]), max(dim_idp_m[idp])+2,2):
                mtotn = mtotal2 - mtotp
                dimp = dim_idp_m[idp][mtotp]
                dimn = dim_idn_m[idn].get(mtotn, 0)
                if dimn == 0: continue
                fp_out.write( struct.pack( 
                    edn + '%id'%(dimp*dimn), 
                    *numpy.array( wfdata_up[mtotp] ).reshape(-1)
                    / ((jj[istate] - mtotal)*(jj[istate] + mtotal + 2) / 4)**.5  ) )


    fp_wav.close()
    fp_out.close()
 
