#!/usr/bin/python
#
usage = \
""" 
 usage: overlap_ptns.py foo.snt bar1.ptn bar1.wav bar2.ptn bar2.wav 
 same snt, M, parity with different truncation scheme
"""

import sys, time, struct, itertools
import numpy as np
from count_dim import mp_add, mp_product, mps_product, set_dim_singlej

fb = '<' # little endian for intel CPU 
# fb = '>' # big endian for SPARC CPU

chunksize = 20 # TODO

def read_snt(filename):
    fp = open(filename,'r')

    arr = fp.readline().split()
    
    while arr[0][0] in ['!','#']:
        arr = fp.readline().split()
    n_jorb = [int(arr[0]),int(arr[1])]
    n_core = [int(arr[2]),int(arr[3])]

    norb, lorb, jorb, itorb = [], [], [], []

    for i in range(sum(n_jorb)):
        arr = fp.readline().split()
        norb.append(int(arr[1]))
        lorb.append(int(arr[2]))
        jorb.append(int(arr[3]))
        itorb.append(int(arr[4]))
    
    fp.close()
    return n_jorb,n_core,norb,lorb,jorb,itorb

            

class Partition:
    # n_ferm, iprty, ptn_p, ptn_n, ptn_pn
    def __init__(self, fn):

        fp = open(fn,'r')

        arr = fp.readline().split()
        while arr[0][0] in ['!','#']: arr = fp.readline().split()
        self.n_ferm = [int(arr[0]),int(arr[1])]
        self.iprty = int(arr[2])
        arr = fp.readline().split()

        while arr[0][0] in ['!','#']: arr = fp.readline().split()
        n_idp = int(arr[0])
        n_idn = int(arr[1])
        arr = fp.readline().split()


        while arr[0][0] in ['!','#']: arr = fp.readline().split()
        self.ptn_p = []
        for i in range(n_idp):
            self.ptn_p.append(tuple([int(x) for x in arr[1:]]))
            arr = fp.readline().split()

        while arr[0][0] in ['!','#']: arr = fp.readline().split()
        self.ptn_n = []
        for i in range(n_idn):
            self.ptn_n.append(tuple([int(x) for x in arr[1:]]))
            arr = fp.readline().split()
            
        while arr[0][0] in ['!','#']: arr = fp.readline().split()
        n_idpn = int(arr[0])
        self.ptn_pn = []
        for i in range(n_idpn):
            arr = fp.readline().split()
            self.ptn_pn.append((int(arr[0])-1, int(arr[1])-1))

        fp.close()

        return

    
    def set_dim(self, mtotal, n_jorb, jorb, lorb):

        dim_jnm = set_dim_singlej( jorb )

        # dimension for proton
        dim_idp_mp = []
        for ptn in self.ptn_p:
            mps = []
            for i,n in enumerate(ptn):
                p = (-1)**(lorb[i]*n) 
                j = jorb[i]
                mps.append( dict( ( (m, p), d ) for m,d in list(dim_jnm[j][n].items()) ) )
            dim_idp_mp.append(mps_product(mps))

        # dimension for neutron
        dim_idn_mp = []
        for ptn in self.ptn_n:
            mps = []
            for i,n in enumerate(ptn):
                p = (-1)**( lorb[ n_jorb[0]+i ] * n )
                j = jorb[ n_jorb[0]+i ]
                mps.append( dict( ( (m, p), d ) for m,d in list(dim_jnm[j][n].items()) ) )
            dim_idn_mp.append( mps_product(mps) )

        # product dimensions of proton and neutron
        self.dim_ptn = []
        for idp,idn in self.ptn_pn:
            dim_mp = {}
            mp_add( dim_mp, mp_product(dim_idp_mp[idp], dim_idn_mp[idn]) )
            self.dim_ptn.append( dim_mp[ mtotal, self.iprty ] )
        self.ndim = sum(self.dim_ptn)


        
if __name__ == '__main__':

    if len(sys.argv) != 6:
        print( usage )
        sys.exit(1)

    start_time = time.time()

    fn_snt, fn_ptn_l, fn_wav_l, fn_ptn_r, fn_wav_r = sys.argv[1:6]

    n_jorb, n_core, norb, lorb, jorb, itorb = read_snt( fn_snt )

    
    ptn_l = Partition( fn_ptn_l )
    ptn_r = Partition( fn_ptn_r )

    
    if ptn_l.n_ferm != ptn_r.n_ferm: sys.exit(1)
    if ptn_l.iprty  != ptn_r.iprty : sys.exit(1)

    print( "snt file : " + fn_snt )
    print( "ptn files: " + fn_ptn_l + " : " + fn_ptn_r )
    print( "wav files: " + fn_wav_l + " : " + fn_wav_r )
    Z = ptn_l.n_ferm[0] + n_core[0]
    N = ptn_l.n_ferm[1] + n_core[1]
    print( "Z=%3i, N=%3i" % ( Z, N ) )
    print( "parity=%2i" % ptn_l.iprty )
    sys.stdout.flush()

    fp_l = open( fn_wav_l, 'rb' )
    fp_r = open( fn_wav_r, 'rb' )

    def read_header(fp, iprty):
        neig, mtotal = struct.unpack( fb+'ii', fp.read(8) )
        evl  = struct.unpack( fb+'%id' % neig, fp.read(8*neig) )
        jj   = struct.unpack( fb+'%ii' % neig, fp.read(4*neig) )
        print( '-------------------------------------------------' )
        print( "neig=%3i, M=%3i/2" % (neig, mtotal) )
        print( "  i     JP       E(MeV)" )
        for i in range(neig):
            print( "%3i %3i/2%s %12.5f"%(i+1,jj[i],"- +"[iprty+1],evl[i]) )
        print()
        return neig, mtotal, evl, jj

    neig_l, mtotal_l, eval_l, jj_l = read_header(fp_l, ptn_l.iprty) 
    neig_r, mtotal_r, eval_r, jj_r = read_header(fp_r, ptn_r.iprty)



    ptn_l.set_dim( mtotal_l, n_jorb, jorb, lorb )
    ptn_r.set_dim( mtotal_r, n_jorb, jorb, lorb )

    print( "left  dim.", ptn_l.ndim )
    print( "right dim.", ptn_r.ndim )

    sys.stdout.flush()


    ptn_matrix = []
    for idl, dim_l in enumerate( ptn_l.dim_ptn  ):
        for idr, dim_r in enumerate( ptn_r.dim_ptn ):

            if dim_l != dim_r: continue
            
            if ptn_l.ptn_p[ ptn_l.ptn_pn[idl][0] ] == \
               ptn_r.ptn_p[ ptn_r.ptn_pn[idr][0] ] and \
               ptn_l.ptn_n[ ptn_l.ptn_pn[idl][1] ] == \
               ptn_r.ptn_n[ ptn_r.ptn_pn[idr][1] ] :
                ptn_matrix.append( ( dim_l, sum( ptn_l.dim_ptn[:idl]), sum( ptn_r.dim_ptn[:idr]) ) )
                break
    print( "partition matrix : ", ptn_matrix)
    print()

    sys.stdout.flush()
    

    olp = np.zeros( (neig_l, neig_r), np.float64)
    
    for dm, dsl, dsr in ptn_matrix:

        wv_l = np.empty( (dm, neig_l), np.float64 )
        for il in range(neig_l):
            fp_l.seek( 8 + 12*neig_l + 8*ptn_l.ndim*il + 8*dsl , 0 )
            wv_l[:,il]= np.asarray( struct.unpack( fb+str(dm)+'d', fp_l.read( 8*dm ) ) )
            
        wv_r = np.empty( (dm, neig_r), np.float64 )
        for ir in range(neig_r):
            fp_r.seek( 8 + 12*neig_r + 8*ptn_r.ndim*ir + 8*dsr , 0)
            wv_r[:,ir] = np.asarray( struct.unpack( fb+str(dm)+'d', fp_r.read( 8*dm ) ) )

        olp += np.dot( np.transpose(wv_l), wv_r )

    for il in range(neig_l):
        for ir in range(neig_r):
            print( "overlap:{:>5},{:>5},  {:>9.6f}".format( il+1, ir+1,  olp[il, ir] ) )
            
    # for il in range(neig_l):
    #     for ir in range(neig_r):
    #         for dm, dsl, dsr in ptn_matrix:

    #             fp_l.seek( 8 + 12*neig_l + 8 * ptn_l.ndim * il + 8 * dsl , 0 )
    #             wv_l = np.asarray( struct.unpack( fb+str(dm)+'d', fp_l.read( 8*dm ) ) )
                
    #             fp_r.seek( 8 + 12*neig_r + 8 * ptn_r.ndim * ir + 8 * dsr , 0)
    #             wv_r = np.asarray( struct.unpack( fb+str(dm)+'d', fp_r.read( 8*dm ) ) )
                
    #             olp[il, ir] += np.dot(  wv_l, wv_r )

    #         print( "overlap:{:>5},{:>5},  {:>9.6f}".format( il+1, ir+1,  olp[il, ir] ) )
    #         sys.stdout.flush()
            
            
    print()
    print( "elapsed time : {:12.3f} sec.".format( time.time() - start_time ) )
    print()

