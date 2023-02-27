#!/usr/bin/python
#
usage = \
""" 
 usage: convert_m2j.py foo.snt bar.ptn input.wav (output.jwv)
"""
#

from functools import reduce


# parameters 

chunk_save = 10000 #  chunk number of partitions to save at once
chunk_calc =   50 #  chunk number of partitions for parallel computation
chunk_calc =   10 #  chunk number of partitions for parallel computation
#chunk_save = 1 #  chunk number of partitions to save at once
#chunk_calc = 1 #  chunk number of partitions for parallel computation


fb = '<' # little endian for intel CPU 
# fb = '>' # big endian for SPARC CPU

nthreads = None
# nthreads = 8

fn_sb = 'seniority_basis.dat' # single-j j-scheme basis states by seniority_basis_exact.py
fb_sb = '<' # little endian



def orbit_coupling( n_jorb, norb, lorb, jorb, itorb ):
    # order of j coupling
    #
    # Ex. order=[(0,1),(2,3),(4,5)] for 4 orbits
    # 0+1 => 4, 2+3 => 5, 4+5 => 6
    # 
    type = 1
    order = []
    if type == 1:
        orbp = 0
        for i in range(1, n_jorb[0]):
            order.append((orbp,i))
            orbp = sum(n_jorb) + len(order) - 1
        orbn = n_jorb[0]
        for i in range(n_jorb[0]+1, sum(n_jorb)):
            order.append((orbn,i))
            orbn = sum(n_jorb) + len(order) - 1
        order.append( (orbp, orbn) )

    print() 
    print( "orbit coupling information" )
    for i in range(sum(n_jorb)):
        print( "{:2}: {} {}{}{:_>2}/2".format(
            i+1, "pn"[(itorb[i]+1)//2], norb[i],
            "spdfghijk"[lorb[i]], jorb[i]) )
    for i,(j,k) in enumerate(order):
        print(  "{:2}: {:2}+{:2}".format(sum(n_jorb)+i+1, j+1, k+1) )

    return order



# ---

import sys, struct, itertools, os, time
import multiprocessing
import numpy


# number of Young diagrams with n boxes (height <= h, width <= w)

young_cache = {}

def young(h, w, n):
    if h<0 or w<0 or n<0:
        return 0
    if (h,w,n) in young_cache:
        return young_cache[h, w, n]
    if n == 0:
        return 1
    young_cache[h, w, n] = young(h, w-1, n-h) + young(h-1, w, n)
    return young_cache[h, w, n]
   
# number of states of n fermions in j orbit with Jz=M is young(2j+1-n,n,(2j+1-n)*n/2)-M)

def calc_dim_vj(j):
    # return dim_vj[v,2J]
    dim_vj={}
    for v in range(j//2+2):
        for jj in range((j+1-v)*v,-1,-2):
            dim =   young(j+1-v,v,((j+1-v)*v-jj)/2) \
                  - young(j+1-v,v,((j+1-v)*v-jj)/2-1) \
                  - young(j+3-v,v-2,((j+3-v)*(v-2)-jj)/2) \
                  + young(j+3-v,v-2,((j+3-v)*(v-2)-jj)/2-1)
            if dim!=0:
                dim_vj[v,jj] = dim
    return dim_vj

dcgvalue = dict()

# functions of rotational group
# taken from rmcsm

def dbinomial(n, m):
    #
    #  binomial coefficient: n_C_m
    #  s: double precision
    #
    s = 1.0
    m1 = min(m, n-m)
    if m1 == 0: return s
    if n < 250:
        s1 = 1.0
        s2 = 1.0
        for i in range(1, m1+1):
            s1 = s1 * (n-i+1)
            s2 = s2 * (m1-i+1)
        s = s1 / s2
    else:
        for i in range(1, m1+1):
            s = (s * (n-i+1)) / (m1-i+1)
    return s

def triangle(j1, j2, j3):
    #
    #  triangle often used in calculation of 3j, 6j etc.
    #  delta
    #  = sqrt(((j1+j2-j3)/2)!((j1-j2+j3)/2)!((-j1+j2+j3)/2)!/(1+(j1+j2+j3)/2)!)
    #
    delta = 0.0
    info = 0
    js = j1 + j2 + j3
    k1 = j2 + j3 - j1
    k2 = j3 + j1 - j2
    k3 = j1 + j2 - j3
    if k1 < 0 or k2 < 0 or k3 <0 or js % 2 !=0:
        info = 1
        return (delta, info)
    js = js // 2
    k1 = k1 // 2
    delta = 1.0 / ( dbinomial(js, j3)**0.5 * dbinomial(j3, k1)**0.5
                    * (js+1.0)**0.5 )
    return (delta, info)

def dcg_calc(j1, m1, j2, m2, j3, m3):
    #
    #  Clebsch-Gordan coefficient 
    #
    #  dcg(j1, m1, j2, m2, j3, m3) 
    #  = ((j1)/2, (m1)/2, (j2)/2, (m2)/2 | (j3)/2, (m3)/2)
    #
    #  using the formula by Racah (1942)
    #
    if (j1, m1, j2, m2, j3, m3) in dcgvalue:
        return dcgvalue[(j1, m1, j2, m2, j3, m3)]
    s = 0
    jm1 = j1 - m1
    jm2 = j2 - m2
    jm3 = j3 - m3
    (delta, info) = triangle(j1, j2, j3)
    if info == 1: return s
    if m3 != m1+m2: return s
    jm1 = jm1 // 2
    jm2 = jm2 // 2
    jm3 = jm3 // 2
    js = (j1 + j2 + j3) // 2
    k1 = (j2 + j3 - j1) // 2
    k2 = (j3 + j1 - j2) // 2
    k3 = (j1 + j2 - j3) // 2
    tmp = ( (  dbinomial(j1, k2)/dbinomial(j1, jm1))**0.5
            * (dbinomial(j2, k3)/dbinomial(j2, jm2))**0.5 
            * (dbinomial(j3, k1)/dbinomial(j3, jm3))**0.5
            * ((j3+1.0))**0.5 * delta )
    izmin = max(0, jm1-k2, k3-jm2)
    izmax = min(k3, jm1, j2-jm2)
    isign = (-1)**izmin
    for iz in range(izmin, izmax+1):
        s = s + isign * dbinomial(k3,iz) \
            * dbinomial(k2,jm1-iz) * dbinomial(k1,j2-jm2-iz)
        isign = isign * (-1)
    s = s * tmp 
    dcgvalue[(j1, m1, j2, m2, j3, m3)] = s
    return s

dcg_matrix_cache = {}

def dcg_matrix(j1, j2, m3):
    if (j1,j2,m3) in dcg_matrix_cache:
        return dcg_matrix_cache[j1,j2,m3]
    size = (j1+j2-max(abs(j1-j2), abs(m3)))//2 + 1
    r = numpy.matrix( numpy.empty((size,)*2) )
    for i in range(size):
        j3 = max(abs(j1-j2), abs(m3)) + i*2
        for j in range(size):
            m1 = max(-j1, m3-j2)+j*2
            m2 = m3 - m1
            r[i,j] = dcg_calc(j1, m1, j2, m2, j3, m3)
    dcg_matrix_cache[j1, j2, m3] = r
    return r
   
def mm2jm(j1, j2, m1, m2):
    m = m1+m2
    if m > j1 - j2:
        j = m1 + j2
    else:
        j = j1 - m2
    return j,m
   
def jm2mm(j1,j2,j,m):
    if m > j1-j2:
        m1 = j - j2
        m2 = m - m1
    else:
        m2 = j1 - j
        m1 = m - m2
    return m1, m2



class ModelSpace:
    def __init__(self, fn_snt):
        # read snt file
        self.n_jorb, self.n_core, self.norb, self.lorb, self.jorb, self.itorb = read_snt( fn_snt )

    def read_ptn(self, fn_ptn):
        # read ptn file
        self.n_ferm, self.iprty, self.ptn_p, self.ptn_n, self.ptn_pn = read_ptn( fn_ptn )

    def read_singlej_basis(self, fn_sn_bs):
        # read seniority_basis.dat

        dim_jnm, mbit_jn, m_j_mbit, mbit_jnm, mbit_jnm_rev = {}, {}, {}, {}, {}

        for j in set(jorb):
            dim_jnm[j]  = [{} for x in range(j+2)]
            mbit_jn[j]  = [[] for x in range(j+2)]
            mbit_jnm[j] = [{} for x in range(j+2)]
            mbit_jnm_rev[j] = [{} for x in range(j+2)]
            m_j_mbit[j] = []
            for mbit in range(2**(j+1)):
                n, m = 0, 0
                for i in range(j+1):
                    n +=  (mbit&2**i)>>i
                    m += ((mbit&2**i)>>i) * (2*i-j)
                dim_jnm[j][n][m] = dim_jnm[j][n].get(m,0) + 1
                mbit_jn[j][n].append(mbit)
                if not m in mbit_jnm[j][n]:
                    mbit_jnm[j][n][m] = []
                    mbit_jnm_rev[j][n][m] = {}
                mbit_jnm_rev[j][n][m][mbit] = len(mbit_jnm[j][n][m])
                mbit_jnm[j][n][m].append(mbit)
                m_j_mbit[j].append(m)

        dim_idp_mp, dim_idp_m, mbit_idp_m, idx_idp_m_mbit \
            = zip(*map(init_mbit_p, ptn_p))
        dim_idn_mp, dim_idn_m, mbit_idn_m, idx_idn_m_mbit \
            = zip(*map(init_mbit_n, ptn_n))

        dim_idpn_mp = []
        dim_mp = {}
        for idp,idn in ptn_pn:
            dim_idpn_mp.append( mp_product(dim_idp_mp[idp], dim_idn_mp[idn] ) )
            mp_add( dim_mp, dim_idpn_mp[-1] )

        dim_j_vj = dict( (j, calc_dim_vj(j)) for j in set(jorb) )
        jva_jn   = dict( (j, [ [ (jj,v,alpha) for (v,jj),dim
                                in sorted( dim_j_vj[j].items(),
                                           key=lambda x:x[0][::-1])
                                 if (n+v)%2==0 and v<=min(n,j+1-n)
                                 for alpha in range(dim) ]
                               for n in range(j+2)] )
                          for j in set(jorb) )
        j_jn     = dict( (j, [sorted(list(set([x[0] for x in jva_jn[j][n]])))
                             for n in range(j+2)])
                         for j in set(jorb) )
        va_jnj   = dict( (j, [dict([(jj,[x[1:] for x in jva_jn[j][n] if x[0]==jj])
                                   for jj in j_jn[j][n]])
                             for n in range(j+2)])
                         for j in set(jorb) )
        dim_jnj  = dict( (j, [dict([(jj,len(va_jnj[j][n][jj]))
                                   for jj in j_jn[j][n]])
                             for n in range(j+2)])
                         for j in set(jorb) )
        jva_jnm  = dict( (j, [dict([(m,[x for x in jva_jn[j][n] if x[0]>=abs(m)])
                                   for m in range(-n*(j+1-n),n*(j+1-n)+2,2)])
                             for n in range(j+2)])
                         for j in set(jorb) )
        jva_jnm_rev = dict( (j, [dict([(m,dict([(y,x) for x,y
                                               in enumerate(jva_jnm[j][n][m])]))
                                      for m in range(-n*(j+1-n),n*(j+1-n)+2,2)])
                                for n in range(j+2)])
                            for j in set(jorb) )


        if not os.path.exists(fn_sb):
            bindir = os.path.dirname( __file__ )
            fn_sb = bindir + '/' + fn
        if not os.path.exists(fn_sb):
            raise "execute seniority_basis_exact.py to make seniority_basis.dat"

        jvabit_jnm = read_basis_all( fn_sb )
        
        def read_basis(f, j):
            jvabit_nm = {}
            for i in range(2**(j+1)):
                n, v, alpha, jj, m, num = struct.unpack( fb_sb+"6i", f.read(24) )
                if not (n,m) in jvabit_nm:
                    jvabit_nm[n,m] = numpy.zeros((len(self.mbit_jnm[j][n][m]),)*2)
                for k in range(num):
                    bit, value = struct.unpack( fb_sb+"id", f.read(12) )
                    jvabit_nm[n,m][ jva_jnm_rev[j][n][m][jj,v,alpha],
                                    mbit_jnm_rev[j][n][m][bit] ] = value
            return jvabit_nm

        fp = open(fn_sn_bs, 'rb')
        jvabit_jnm = {}
        is_maxj = False
        while True:
            j = fp.read(4)
            if not j: break
            j = struct.unpack(fb_sb+"i", j)[0]
            if j > max(self.jorb): break
            if j == max(self.jorb): is_maxj = True
            jvabit_jnm.update(
                dict([((j,)+x,y) for x,y in read_basis(f,j).items()]) )
        fp.close()
        if not is_maxj:
            raise 'maxj fail: increase maxj in seniority_basis_exact.py and run '
        
        self.jvabit_jnm = jvabit_jnm
        
        return 
        

   
def read_basis(f, j):
    jvabit_nm = {}
    for i in range(2**(j+1)):
        n, v, alpha, jj, m, num = struct.unpack( fb_sb+"6i", f.read(24) )
        if not (n,m) in jvabit_nm:
            jvabit_nm[n,m] = numpy.zeros((len(mbit_jnm[j][n][m]),)*2)
        for k in range(num):
            bit, value = struct.unpack( fb_sb+"id", f.read(12) )
            jvabit_nm[n,m][ jva_jnm_rev[j][n][m][jj,v,alpha],
                            mbit_jnm_rev[j][n][m][bit] ] = value
    return jvabit_nm


def read_basis_all(fn):
    f = open(fn, 'rb')
    jvabit_jnm = {}
    is_maxj = False
    while True:
        j = f.read(4)
        if not j: break
        j = struct.unpack(fb_sb+"i", j)[0]
        if j > max(jorb): break
        if j == max(jorb): is_maxj = True
        jvabit_jnm.update(
            dict([((j,)+x,y) for x,y in read_basis(f,j).items()]) )
    f.close()
    if not is_maxj:
        raise 'maxj fail: increase maxj in seniority_basis_exact.py and run '
    return jvabit_jnm
   
   
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

   
def read_ptn(filename):

    fp = open(filename,'r')
    
    arr = fp.readline().split()
    while arr[0][0] in ['!','#']: arr = fp.readline().split()
    n_ferm = [int(arr[0]),int(arr[1])]
    iprty = int(arr[2])
    arr = fp.readline().split()
    
    while arr[0][0] in ['!','#']: arr = fp.readline().split()
    n_idp = int(arr[0])
    n_idn = int(arr[1])
    arr = fp.readline().split()

    
    while arr[0][0] in ['!','#']: arr = fp.readline().split()
    ptn_p = []
    for i in range(n_idp):
        ptn_p.append(tuple([int(x) for x in arr[1:]]))
        arr = fp.readline().split()
        
    while arr[0][0] in ['!','#']: arr = fp.readline().split()
    ptn_n = []
    for i in range(n_idn):
        ptn_n.append(tuple([int(x) for x in arr[1:]]))
        arr = fp.readline().split()
        
    while arr[0][0] in ['!','#']: arr = fp.readline().split()
    n_idpn = int(arr[0])
    ptn_pn = []
    for i in range(n_idpn):
        arr = fp.readline().split()
        ptn_pn.append((int(arr[0])-1, int(arr[1])-1))
        
    fp.close()
    return n_ferm, iprty, ptn_p, ptn_n, ptn_pn


   
def mp_add(mp1,mp2):
    for (m,p),value in mp2.items():
        mp1[m,p] = mp1.get((m,p),0)+value

        
def mp_product(mp1,mp2):
    mp = dict()
    for (m1,p1),value1 in mp1.items():
        for (m2,p2),value2 in mp2.items():
            mp[m1+m2,p1*p2] = mp.get((m1+m2,p1*p2),0) + value1*value2
    return mp
   
def mps_product(mps):
    while True:
        if len(mps)==1: break
        mp1 = mps.pop(0)
        mp2 = mps.pop(0)
        mps.append(mp_product(mp1,mp2))
    return mps[0]
   
def init_mbit_p(ptn):
    mps = [dict([((m,(-1)**(lorb[i]*n)),dim_jnm[jorb[i]][n][m])
               for m in dim_jnm[jorb[i]][n]])
         for i,n in enumerate(ptn)]
    dim_mp = mps_product(mps)
    dim_m = dict([(m,x) for (m,p),x in dim_mp.items()])
    mbit_m = {}
    idx_m_mbit = {}
    for mbits in itertools.product( *[mbit_jn[x][y] for x,y
                                      in zip(jorb[:n_jorb[0]],ptn)][::-1] ):
        m = sum([m_j_mbit[x][y] for x,y in zip(jorb[:n_jorb[0]],mbits[::-1])])
        if not m in mbit_m:
            mbit_m[m] = []
            idx_m_mbit[m] = {}
        mbit_m[m].append(mbits[::-1])
        idx_m_mbit[m][mbits[::-1]] = len(mbit_m[m])-1
    return dim_mp,dim_m,mbit_m,idx_m_mbit
   
def init_mbit_n(ptn):
    mps = [dict([ ( (m,(-1)**(lorb[n_jorb[0]+i]*n)),
                    dim_jnm[jorb[n_jorb[0]+i]][n][m] )
                  for m in dim_jnm[jorb[n_jorb[0]+i]][n]] )
         for i,n in enumerate(ptn)]
    dim_mp = mps_product(mps)
    dim_m = dict([(m,x) for (m,p),x in dim_mp.items()])
    mbit_m = {}
    idx_m_mbit = {}
    for mbits in itertools.product(*[mbit_jn[x][y] for x,y
                                     in zip(jorb[n_jorb[0]:],ptn)][::-1]):
        m = sum( [ m_j_mbit[x][y] for x,y
                   in zip(jorb[n_jorb[0]:],mbits[::-1]) ] )
        if not m in mbit_m:
            mbit_m[m] = []
            idx_m_mbit[m] = {}
        mbit_m[m].append(mbits[::-1])
        idx_m_mbit[m][mbits[::-1]] = len(mbit_m[m])-1
    return dim_mp, dim_m, mbit_m, idx_m_mbit
   
def wrapper_convert_wf(x):
    return convert_wf_prl(x[0],**x[1])
   
def convert_wf_prl( idpn, wavdata, ptn_pn, jorb, occs, mtotal,
                    mbit_jnm, dim_idp_m, dim_idn_m,
                    mbit_idp_m, mbit_idn_m, m_j_mbit, mbit_jnm_rev,
                    n_jorb, jvabit_jnm, j_jn, va_jnj, dim_jnj, jva_jnm, jj):
    
    # print >> sys.stderr, idpn+1,"/",len(ptn_pn)
    # sys.stdout.write( '\r {:8d} / {:8d}'.format( idpn+1, len(ptn_pn) ) )
    # sys.stdout.flush()
   
    mslist = [ tuple(ms) for ms
               in itertools.product( *[ range(-n*(j+1-n), n*(j+1-n)+2, 2)
                                        for j,n in zip(jorb,occs) ] )
               if sum(ms) == mtotal ]
    wfdata = dict([ (ms, numpy.empty([ len(mbit_jnm[j][n][m]) for j,n,m
                                       in zip(jorb,occs,ms) ]) )
                    for ms in mslist ])
    pos = 0
    for mtotp in range(min(dim_idp_m),max(dim_idp_m)+2,2):
        mtotn = mtotal - mtotp
        dimp = dim_idp_m[mtotp]
        dimn = dim_idn_m.get(mtotn, 0)
        for idimn in range(dimn):
            for idimp in range(dimp):
                mbits = mbit_idp_m[mtotp][idimp] + \
                        mbit_idn_m[mtotn][idimn]
                ms = tuple([m_j_mbit[j][mbit] for j,mbit in zip(jorb,mbits)])
                wfdata[ms][
                    tuple([ mbit_jnm_rev[j][n][m][mbit]
                            for j,n,m,mbit
                            in zip(jorb,occs,ms,mbits)]) ] \
                            = struct.unpack(fb+'d', wavdata[pos:pos+8])[0]
                pos += 8
            
    for ms in wfdata:
        for i in range(sum(n_jorb)):
            wfdata[ms] = numpy.einsum(
                wfdata[ms],
                list(range(i)) + [sum(n_jorb)] + list(range(i+1, sum(n_jorb))),
                jvabit_jnm[jorb[i], occs[i], ms[i]],
                [i, sum(n_jorb)],
                range(sum(n_jorb)) )
         
    jslist = [ tuple(js) for js in
               itertools.product( *[j_jn[j][n] for j,n in zip(jorb,occs)] ) ]
    ms_js_rev = dict([ (js,dict([
        (tuple(ms),i) for i,ms
        in enumerate([ x for x in
                       itertools.product(*[range(-j,j+2,2) for j in js])
                       if sum(x)==mtotal])]))
                       for js in jslist ])
    vas_js_rev = dict([ (js,dict([
        (tuple(vas[0]+vas[1]),i) for i,vas
        in enumerate(sorted([list(zip(*x)) for x in
                             itertools.product(*[va_jnj[j][n][j2] for j,n,j2 in
                                                 zip(jorb,occs,js)])]))] ))
                       for js in jslist ])
    wfdata_new = dict([ (js, numpy.empty( (len(ms_js_rev[js]),
                                           reduce(lambda x,y:x*y,
                                                  [dim_jnj[j][n][j2]
                                                   for j,n,j2 in
                                                   zip(jorb,occs,js)]) ) ) )
                        for js in jslist])
    for ms in wfdata:
        for jvas,value in zip( itertools.product(*[jva_jnm[j][n][m]
                                                   for j,n,m in
                                                   zip(jorb,occs,ms)]),
                               wfdata[ms].flat ):
            js, vs, alphas = zip(*jvas)
            js = tuple(js)
            vas = tuple(vs+alphas)
            wfdata_new[js][ms_js_rev[js][ms], vas_js_rev[js][vas]] = value
    wfdata = wfdata_new
   
    r = b""
    for js in sorted(wfdata):
        for i in range(n_jorb[0]-1):
            if js[i+1]==0: continue
            if i==0 and js[0]==0: continue
            # initialize
            jcs=[]
            for i2 in range(i):
                jcs.append(abs((js[0] if i2==0 else jcs[i2-1])-js[i2+1]))
            msp = [ -js[i2] for i2 in range(i+2,n_jorb[0]) ]
            msn = [ -js[i2] for i2 in range(n_jorb[0],sum(n_jorb)) ]
            while True:
                # main
                m = mtotal - sum(msp + msn)
                j1 = (js[0] if i==0 else jcs[i-1])
                j2 = js[i+1]
                if j1 != 0 and abs(m) < j1+j2:
                    idx = []
                    for m1 in range(max(-j1,m-j2), min(j1,m+j2)+2, 2):
                        m2 = m - m1
                        ms = jcs+[m1,m2] + msp + msn
                        for i2 in range(i-1,-1,-1):
                            ms = ms[:i2] + \
                                 list(jm2mm(js[0] if i2==0 else ms[i2-1],
                                            js[i2+1], ms[i2], ms[i2+1])) \
                                 + ms[i2+2:]
                        idx.append(ms_js_rev[js][tuple(ms)])
                    wfdata[js][idx,:] = dcg_matrix(j1,j2,m)*wfdata[js][idx,:]
                #increment
                for i2 in range(sum(n_jorb)-1, n_jorb[0]-1, -1):
                    if j1 == 0: continue
                    if msn[i2-n_jorb[0]] == js[i2]: continue
                    msn[i2-n_jorb[0]] += 2
                    for i3 in range(i2+1, sum(n_jorb)):
                        msn[i3-n_jorb[0]] = -js[i3]
                    break
                else:
                    for i2 in range(n_jorb[0]-1, i+1, -1):
                        if j1 == 0: continue
                        if msp[i2-i-2] == js[i2]: continue
                        msp[i2-i-2] += 2
                        for i3 in range(i2+1,n_jorb[0]):
                            msp[i3-i-2] = -js[i3]
                        for i3 in range(n_jorb[0],sum(n_jorb)):
                            msn[i3-n_jorb[0]] = -js[i3]
                        break
                    else:
                        for i2 in range(i-1,-1,-1):
                            if jcs[i2]==( js[0] if i2==0 else
                                          jcs[i2-1])+js[i2+1]: continue
                            jcs[i2] += 2
                            for i3 in range(i2+1, i):
                                jcs[i3] = abs(jcs[i3-1] - js[i3+1])
                            for i3 in range(i+2, n_jorb[0]):
                                msp[i3-i-2] = -js[i3]
                            for i3 in range(n_jorb[0], sum(n_jorb)):
                                msn[i3-n_jorb[0]] = -js[i3]
                            break
                        else:
                            break
                     
        for i in range(n_jorb[0], sum(n_jorb)-1):
            if js[i+1] == 0: continue
            if i==n_jorb[0] and js[n_jorb[0]]==0: continue
            #initialize
            jcsp=[]
            for i2 in range(n_jorb[0]-1):
                jcsp.append(abs((js[0] if i2==0 else jcsp[i2-1])-js[i2+1]))
            jcsn=[]
            for i2 in range(n_jorb[0],i):
                jcsn.append(abs((js[n_jorb[0]]
                                 if i2==n_jorb[0]
                                 else jcsn[i2-n_jorb[0]-1])-js[i2+1]))
            mp = -jcsp[-1]
            msn = [-js[i2] for i2 in range(i+2,sum(n_jorb))]
            while True:
                #main
                m  = mtotal - sum([mp] + msn)
                j1 = (js[n_jorb[0]] if i==n_jorb[0] else jcsn[i-n_jorb[0]-1])
                j2 = js[i+1]
                if j1!=0 and abs(m)<j1+j2:
                    idx = []
                    for m1 in range(max(-j1,m-j2),min(j1,m+j2)+2,2):
                        m2 = m - m1
                        ms = jcsp + [mp] + jcsn + [m1,m2] + msn
                        for i2 in range(i-1, n_jorb[0]-1, -1):
                            ms=ms[:i2]+list(jm2mm(js[n_jorb[0]] if i2==n_jorb[0] else ms[i2-1],js[i2+1],ms[i2],ms[i2+1]))+ms[i2+2:]
                        for i2 in range(n_jorb[0]-2,-1,-1):
                            ms=ms[:i2]+list(jm2mm(js[0] if i2==0 else ms[i2-1],js[i2+1],ms[i2],ms[i2+1]))+ms[i2+2:]
                        idx.append(ms_js_rev[js][tuple(ms)])
                    wfdata[js][idx,:]=dcg_matrix(j1,j2,m)*wfdata[js][idx,:]
                #increment
                for i2 in range(sum(n_jorb)-1,i+1,-1):
                    if j1==0: continue
                    if msn[i2-i-2]==js[i2]: continue
                    msn[i2-i-2]+=2
                    for i3 in range(i2+1,sum(n_jorb)):
                        msn[i3-i-2]=-js[i3]
                    break
                else:
                    if j1!=0 and mp!=jcsp[-1]:
                        mp+=2
                        for i3 in range(i+2,sum(n_jorb)):
                            msn[i3-i-2]=-js[i3]
                    else:
                        for i2 in range(n_jorb[0]-2,-1,-1):
                            if j1==0: continue
                            if jcsp[i2]==(js[0] if i2==0 else jcsp[i2-1])+js[i2+1]: continue
                            jcsp[i2]+=2
                            for i3 in range(i2+1,n_jorb[0]-1):
                                jcsp[i3]=abs(jcsp[i3-1]-js[i3+1])
                            mp=-jcsp[-1]
                            for i3 in range(i+2,sum(n_jorb)):
                                msn[i3-i-2]=-js[i3]
                            break
                        else:
                            for i2 in range(i-1,n_jorb[0]-1,-1):
                                if jcsn[i2-n_jorb[0]]==(js[n_jorb[0]] if i2==n_jorb[0] else jcsn[i2-n_jorb[0]-1])+js[i2+1]: continue
                                jcsn[i2-n_jorb[0]]+=2
                                for i3 in range(i2+1,i):
                                    jcsn[i3-n_jorb[0]]=abs(jcsn[i3-n_jorb[0]-1]-js[i3+1])
                                for i3 in range(n_jorb[0]-1):
                                    jcsp[i3]=abs((js[0] if i3==0 else jcsp[i3-1])-js[i3+1])
                                mp=-jcsp[-1]
                                for i3 in range(i+2,sum(n_jorb)):
                                    msn[i3-i-2]=-js[i3]
                                break
                            else:
                                break
                        
        #initialize
        # jcsp = [ abs((js[0] if i2==0 else jcsp[i2-1])-js[i2+1]) for i2 in range(n_jorb[0]-1) ]
        jcsp=[]
        for i2 in range(n_jorb[0]-1):
            jcsp.append(abs((js[0] if i2==0 else jcsp[i2-1])-js[i2+1]))
        jcsn=[]
        for i2 in range(n_jorb[0],sum(n_jorb)-1):
            jcsn.append(abs((js[n_jorb[0]] if i2==n_jorb[0] else jcsn[i2-n_jorb[0]-1])-js[i2+1]))
        while True:
            #main
            m = mtotal
            j1 = jcsp[-1]
            j2 = jcsn[-1]
            if abs(j1-j2)<=jj<=j1+j2:
                idx = []
                for m1 in range(max(-j1,m-j2), min(j1,m+j2)+2, 2):
                    m2 = m - m1
                    ms = jcsp + [m1] + jcsn + [m2]
                    for i2 in range(sum(n_jorb)-2,n_jorb[0]-1,-1):
                        ms = ms[:i2] + list(jm2mm( js[n_jorb[0]] if i2==n_jorb[0]
                                                   else ms[i2-1],js[i2+1],ms[i2],ms[i2+1])) + ms[i2+2:]
                    for i2 in range(n_jorb[0]-2,-1,-1):
                        ms = ms[:i2]+list(jm2mm(js[0] if i2==0 else ms[i2-1],js[i2+1],ms[i2],ms[i2+1]))+ms[i2+2:]
                    idx.append(ms_js_rev[js][tuple(ms)])
                value = tuple((dcg_matrix(j1,j2,m)[(jj-max(abs(j1-j2),abs(m)))//2,:]*wfdata[js][idx,:]).flat)
                r += struct.pack( '<%id'%len(value), *value )
            #increment
            for i2 in range(sum(n_jorb)-2,n_jorb[0]-1,-1):
                if jcsn[i2-n_jorb[0]] == (js[n_jorb[0]] if i2==n_jorb[0] else jcsn[i2-n_jorb[0]-1]) + js[i2+1]:
                    continue
                jcsn[i2-n_jorb[0]]+=2
                for i3 in range(i2+1,sum(n_jorb)-1):
                    jcsn[i3-n_jorb[0]]=abs(jcsn[i3-n_jorb[0]-1]-js[i3+1])
                break
            else:
                for i2 in range(n_jorb[0]-2,-1,-1):
                    if jcsp[i2]==(js[0] if i2==0 else jcsp[i2-1])+js[i2+1]: continue
                    jcsp[i2]+=2
                    for i3 in range(i2+1,n_jorb[0]-1):
                        jcsp[i3]=abs(jcsp[i3-1]-js[i3+1])
                    for i3 in range(n_jorb[0],sum(n_jorb)-1):
                        jcsn[i3-n_jorb[0]]=abs((js[n_jorb[0]] if i3==n_jorb[0] else jcsn[i3-n_jorb[0]-1])-js[i3+1])
                    break
                else:
                    break
    return r
   

if __name__ == '__main__':

    if len(sys.argv) < 4:
        print( usage )
        sys.exit(1)

    start_time = time.time()

    fn_snt, fn_ptn, fn_wav = sys.argv[1:4]
    if len(sys.argv) == 4:
        fn_jwv = fn_wav[:-4] + '.jwv'
    else:
        fn_jwv = sys.argv[4]
        
    n_jorb, n_core, norb, lorb, jorb, itorb = read_snt( fn_snt )
    n_ferm, iprty, ptn_p, ptn_n, ptn_pn = read_ptn( fn_ptn )

    print( "snt file: " + fn_snt )
    print( "ptn file: " + fn_ptn )
    print( "wav file: " + fn_wav )
    print( "Z=%3i, N=%3i" % (n_ferm[0]+n_core[0],n_ferm[1]+n_core[1]) )
    print( "parity=%2i" % iprty )
    sys.stdout.flush()

    # sm = ModelSpace(fn_snt)
    # sm.read_ptn( fn_ptn )
    # sm.read_singlej_basis( fn_sb )
    
   
    dim_jnm, mbit_jn, m_j_mbit, mbit_jnm, mbit_jnm_rev = {}, {}, {}, {}, {}
    
#    for j in set(jorb):
    for j in range(1, max(jorb)+2, 2):
        dim_jnm[j]  = [{} for x in range(j+2)]
        mbit_jn[j]  = [[] for x in range(j+2)]
        mbit_jnm[j] = [{} for x in range(j+2)]
        mbit_jnm_rev[j] = [{} for x in range(j+2)]
        m_j_mbit[j] = []
        for mbit in range(2**(j+1)):
            n, m = 0, 0
            for i in range(j+1):
                n +=  (mbit&2**i)>>i
                m += ((mbit&2**i)>>i) * (2*i-j)
            dim_jnm[j][n][m] = dim_jnm[j][n].get(m,0) + 1
            mbit_jn[j][n].append(mbit)
            if not m in mbit_jnm[j][n]:
                mbit_jnm[j][n][m] = []
                mbit_jnm_rev[j][n][m] = {}
            mbit_jnm_rev[j][n][m][mbit] = len(mbit_jnm[j][n][m])
            mbit_jnm[j][n][m].append(mbit)
            m_j_mbit[j].append(m)
         
    dim_idp_mp, dim_idp_m, mbit_idp_m, idx_idp_m_mbit \
        = zip(*map(init_mbit_p, ptn_p))
    dim_idn_mp, dim_idn_m, mbit_idn_m, idx_idn_m_mbit \
        = zip(*map(init_mbit_n, ptn_n))
   
    dim_idpn_mp = []
    dim_mp = {}
    for idp,idn in ptn_pn:
        dim_idpn_mp.append( mp_product(dim_idp_mp[idp], dim_idn_mp[idn] ) )
        mp_add( dim_mp, dim_idpn_mp[-1] )
      
    dim_j_vj = dict( (j, calc_dim_vj(j)) for j in set(jorb) )
    jva_jn   = dict( (j, [ [ (jj,v,alpha) for (v,jj),dim
                            in sorted( dim_j_vj[j].items(),
                                       key=lambda x:x[0][::-1])
                             if (n+v)%2==0 and v<=min(n,j+1-n)
                             for alpha in range(dim) ]
                           for n in range(j+2)] )
                      for j in set(jorb) )
    j_jn     = dict( (j, [sorted(list(set([x[0] for x in jva_jn[j][n]])))
                         for n in range(j+2)])
                     for j in set(jorb) )
    va_jnj   = dict( (j, [dict([(jj,[x[1:] for x in jva_jn[j][n] if x[0]==jj])
                               for jj in j_jn[j][n]])
                         for n in range(j+2)])
                     for j in set(jorb) )
    dim_jnj  = dict( (j, [dict([(jj,len(va_jnj[j][n][jj]))
                               for jj in j_jn[j][n]])
                         for n in range(j+2)])
                     for j in set(jorb) )
    jva_jnm  = dict( (j, [dict([(m,[x for x in jva_jn[j][n] if x[0]>=abs(m)])
                               for m in range(-n*(j+1-n),n*(j+1-n)+2,2)])
                         for n in range(j+2)])
                     for j in set(jorb) )
    jva_jnm_rev = dict( (j, [dict([(m,dict([(y,x) for x,y
                                           in enumerate(jva_jnm[j][n][m])]))
                                  for m in range(-n*(j+1-n),n*(j+1-n)+2,2)])
                            for n in range(j+2)])
                        for j in set(jorb) )

    
    if not os.path.exists(fn_sb):
        bindir = os.path.dirname( __file__ )
        fn_sb = bindir + '/' + fn_sb
    if not os.path.exists(fn_sb):
        print( fn_sb )
        raise "execute seniority_basis_exact.py to make seniority_basis.dat"
    
    jvabit_jnm = read_basis_all( fn_sb )
  
    wavfile = open( fn_wav, 'rb' )
    neig, mtotal = struct.unpack(fb+'ii', wavfile.read(8))
    eval = struct.unpack(fb+'%id' % neig,wavfile.read(8*neig))
    jj   = struct.unpack(fb+'%ii' % neig,wavfile.read(4*neig))
    print( '-------------------------------------------------' )
    print( "neig=%3i, M=%3i/2" % (neig, mtotal) )
    print( "  i     JP       E(MeV)" )
    for i in range(neig):
        print( "%3i %3i/2%s %12.5f"%(i+1,jj[i],"- +"[iprty+1],eval[i]) )
    print()
    wavoutfile = open( fn_jwv, 'wb')
    outputlist = [x for x in range(neig) if not jj[x] in [-1]]
    neig2 = len(outputlist)
    mtotal2 = mtotal
    print( "wav output file: " + fn_jwv )
    print( "neig=%3i, M=%3i/2" % (neig2, mtotal2) )
    print( "  i     JP       E(MeV)" )
    for j,i in enumerate(outputlist):
        print( "%3i %3i/2%s %12.5f" % (j+1, jj[i], "- +"[iprty+1], eval[i]) )
    wavoutfile.write(struct.pack('<ii', neig2, mtotal2))
    wavoutfile.write(struct.pack('<%id'%neig2, *[eval[x] for x in outputlist]))
    wavoutfile.write(struct.pack('<%ii'%neig2, *[jj[x] for x in outputlist]))
   
    print( "convert wave functions" )
    sys.stdout.flush()

    p = multiprocessing.Pool(nthreads)

    def input_dict(x):
        return  { "wavdata":wavdata[x], "ptn_pn":ptn_pn,
                  "jorb":jorb, "occs":ptn_p[ptn_pn[x][0]]+ptn_n[ptn_pn[x][1]],
                  "mtotal":mtotal, "mbit_jnm":mbit_jnm,
                  "dim_idp_m":dim_idp_m[ptn_pn[x][0]],
                  "dim_idn_m":dim_idn_m[ptn_pn[x][1]],
                  "mbit_idp_m":mbit_idp_m[ptn_pn[x][0]],
                  "mbit_idn_m":mbit_idn_m[ptn_pn[x][1]],
                  "m_j_mbit":m_j_mbit,
                  "mbit_jnm_rev":mbit_jnm_rev,
                  "n_jorb":n_jorb, "jvabit_jnm":jvabit_jnm,
                  "j_jn":j_jn, "va_jnj":va_jnj, "dim_jnj":dim_jnj,
                  "jva_jnm":jva_jnm, "jj":jj[istate] }


    for istate in range(neig):
        
        print( '\n w.f. state {0:5d} / {1:5d}'.format( istate+1, neig ) )
        sys.stdout.flush()
        
        if not istate in outputlist:
            wavfile.seek( 8*dim_mp[mtotal, iprty], 1)
            continue
        for idpn in range(0, len(ptn_pn), chunk_save):
            sys.stdout.write( '\r {0:8d} / {1:8d}'.format(
                idpn+1, len(ptn_pn) ) )
            sys.stdout.flush()


            wavdata = {}
            for idpn2 in range(idpn, min(idpn+chunk_save, len(ptn_pn))):
                wavdata[idpn2] = wavfile.read(
                    8*dim_idpn_mp[idpn2].get((mtotal, iprty), 0) )
            # single process
            # wavoutdata = "".join( map(
            #     wrapper_convert_wf,[
            #         ( x, input_dict(x) )
            #         for x in range(idpn, min(idpn+chunk_save, len(ptn_pn))) ] ))
            wavoutdata = b"".join( p.map(
                wrapper_convert_wf,[
                    ( x, input_dict(x) )
                    for x in range(idpn, min(idpn+chunk_save, len(ptn_pn))) ],
                chunksize = chunk_calc ) )
            wavoutfile.write( wavoutdata )

        sys.stdout.write( '\r {0:8d} / {1:8d}'.format(len(ptn_pn), len(ptn_pn)) )
        sys.stdout.flush()
        

    t = time.time() - start_time
    th, t = divmod(t, 3600)
    tm, t = divmod(t, 60)
    print( '\nelapsed time  {0:02d}:{1:02d}:{2:05.2f} '.format( int(th), int(tm), t ) )
    print( 'finish.\n' )



