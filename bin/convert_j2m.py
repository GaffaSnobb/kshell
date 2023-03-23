#!/usr/bin/python
#
# usage: convert_m2j.py foo.snt bar.ptn input.wav output.jwav
#

# parameters

fn = 'seniority_basis.dat'

import sys, os, struct, itertools, numpy, multiprocessing
from functools import reduce

from convert_m2j import young, calc_dim_vj, dbinomial, triangle, dcg_calc, dcg_matrix, mm2jm, jm2mm, read_basis, read_basis_all, read_snt

def read_basis(f,j):
    jvabit_nm={}
    for i in range(2**(j+1)):
        n,v,alpha,jj,m,num=struct.unpack("<6i",f.read(24))
        if not (n,m) in jvabit_nm:
            jvabit_nm[n,m]=numpy.zeros((len(mbit_jnm[j][n][m]),)*2)
        for k in range(num):
            bit,value=struct.unpack("<id",f.read(12))
            jvabit_nm[n,m][jva_jnm_rev[j][n][m][jj,v,alpha],mbit_jnm_rev[j][n][m][bit]]=value
    return jvabit_nm
   
def read_basis_all(fn):
    f=open(fn,'rb')
    jvabit_jnm={}
    while True:
        j=f.read(4)
        if not j: break
        j=struct.unpack("<i",j)[0]
        if j>max(jorb): break
        jvabit_jnm.update(dict([((j,)+x,y) for x,y in read_basis(f,j).items()]))
    f.close()
    return jvabit_jnm
   
   
def read_snt(filename):
    inputfile=open(filename,'r')
    arr=inputfile.readline().split()
    while arr[0][0] in ['!','#']:
        arr=inputfile.readline().split()
    n_jorb=[int(arr[0]),int(arr[1])]
    n_core=[int(arr[2]),int(arr[3])]
    norb=[]
    lorb=[]
    jorb=[]
    itorb=[]
    for i in range(sum(n_jorb)):
        arr=inputfile.readline().split()
        norb.append(int(arr[1]))
        lorb.append(int(arr[2]))
        jorb.append(int(arr[3]))
        itorb.append(int(arr[4]))
    inputfile.close()
    return n_jorb,n_core,norb,lorb,jorb,itorb
   
def read_ptn(filename):
    inputfile=open(filename,'r')
    arr=inputfile.readline().split()
    while arr[0][0] in ['!','#']:
        arr=inputfile.readline().split()
    n_ferm=[int(arr[0]),int(arr[1])]
    iprty=int(arr[2])
    arr=inputfile.readline().split()
    while arr[0][0] in ['!','#']:
        arr=inputfile.readline().split()
    n_idp=int(arr[0])
    n_idn=int(arr[1])
    arr=inputfile.readline().split()
    while arr[0][0] in ['!','#']:
        arr=inputfile.readline().split()
    ptn_p=[]
    for i in range(n_idp):
        ptn_p.append(tuple([int(x) for x in arr[1:]]))
        arr=inputfile.readline().split()
    while arr[0][0] in ['!','#']:
        arr=inputfile.readline().split()
    ptn_n=[]
    for i in range(n_idn):
        ptn_n.append(tuple([int(x) for x in arr[1:]]))
        arr=inputfile.readline().split()
    while arr[0][0] in ['!','#']:
        arr=inputfile.readline().split()
    n_idpn=int(arr[0])
    ptn_pn=[]
    for i in range(n_idpn):
        arr=inputfile.readline().split()
        ptn_pn.append((int(arr[0])-1,int(arr[1])-1))
    inputfile.close()
    return n_ferm,iprty,ptn_p,ptn_n,ptn_pn
   
def mp_add(mp1,mp2):
    for (m,p),value in mp2.items():
        mp1[m,p]=mp1.get((m,p),0)+value
      
def mp_product(mp1,mp2):
    mp=dict()
    for (m1,p1),value1 in mp1.items():
        for (m2,p2),value2 in mp2.items():
            mp[m1+m2,p1*p2]=mp.get((m1+m2,p1*p2),0)+value1*value2
    return mp
   
def mps_product(mps):
    while 1:
        if len(mps)==1: break
        mp1=mps.pop(0)
        mp2=mps.pop(0)
        mps.append(mp_product(mp1,mp2))
    return mps[0]
   
def init_mbit_p(ptn):
    mps=[dict([((m,(-1)**(lorb[i]*n)),dim_jnm[jorb[i]][n][m])
               for m in dim_jnm[jorb[i]][n]])
         for i,n in enumerate(ptn)]
    dim_mp=mps_product(mps)
    dim_m=dict([(m,x) for (m,p),x in dim_mp.items()])
    mbit_m={}
    idx_m_mbit={}
    for mbits in itertools.product(*[mbit_jn[x][y] for x,y in zip(jorb[:n_jorb[0]],ptn)][::-1]):
        m=sum([m_j_mbit[x][y] for x,y in zip(jorb[:n_jorb[0]],mbits[::-1])])
        if not m in mbit_m:
            mbit_m[m]=[]
            idx_m_mbit[m]={}
        mbit_m[m].append(mbits[::-1])
        idx_m_mbit[m][mbits[::-1]]=len(mbit_m[m])-1
    return dim_mp,dim_m,mbit_m,idx_m_mbit
   
def init_mbit_n(ptn):
    mps=[dict([((m,(-1)**(lorb[n_jorb[0]+i]*n)),dim_jnm[jorb[n_jorb[0]+i]][n][m])
               for m in dim_jnm[jorb[n_jorb[0]+i]][n]])
         for i,n in enumerate(ptn)]
    dim_mp=mps_product(mps)
    dim_m=dict([(m,x) for (m,p),x in dim_mp.items()])
    mbit_m={}
    idx_m_mbit={}
    for mbits in itertools.product(*[mbit_jn[x][y] for x,y in zip(jorb[n_jorb[0]:],ptn)][::-1]):
        m=sum([m_j_mbit[x][y] for x,y in zip(jorb[n_jorb[0]:],mbits[::-1])])
        if not m in mbit_m:
            mbit_m[m]=[]
            idx_m_mbit[m]={}
        mbit_m[m].append(mbits[::-1])
        idx_m_mbit[m][mbits[::-1]]=len(mbit_m[m])-1
    return dim_mp,dim_m,mbit_m,idx_m_mbit
   
def wrapper_convert_wf(x):
    return convert_wf_prl(x[0],**x[1])
   
def convert_wf_prl(idpn,wavdata,ptn_pn,jorb,occs,mtotal2,mbit_jnm,dim_idp_m,dim_idn_m,mbit_idp_m,mbit_idn_m,m_j_mbit,mbit_jnm_rev,n_jorb,jvabit_jnm,j_jn,va_jnj,dim_jnj,jva_jnm,jj):

    sys.stdout.write( '\r {:8d} / {:8d}'.format( idpn+1, len(ptn_pn) ) )
    sys.stdout.flush()

   
    jslist=[tuple(js) for js in itertools.product(*[j_jn[j][n] for j,n in zip(jorb,occs)])]
    ms_js_rev=dict([(js,dict([(tuple(ms),i) for i,ms in enumerate([x for x in itertools.product(*[range(-j,j+2,2) for j in js]) if sum(x)==mtotal2])])) for js in jslist])
    vas_js_rev=dict([(js,dict([(tuple(vas[0]+vas[1]),i) 
                               for i,vas in enumerate(sorted([list(zip(*x)) for x in itertools.product(*[va_jnj[j][n][j2] for j,n,j2 in zip(jorb,occs,js)])]))]))
                     for js in jslist])
    wfdata=dict([(js,numpy.zeros((len(ms_js_rev[js]),reduce(lambda x,y:x*y,[dim_jnj[j][n][j2] for j,n,j2 in zip(jorb,occs,js)])))) for js in jslist])
    pos=0
    for js in sorted(wfdata):
        #initialize
        jcsp=[]
        for i2 in range(n_jorb[0]-1):
            jcsp.append(abs((js[0] if i2==0 else jcsp[i2-1])-js[i2+1]))
        jcsn=[]
        for i2 in range(n_jorb[0],sum(n_jorb)-1):
            jcsn.append(abs((js[n_jorb[0]] if i2==n_jorb[0] else jcsn[i2-n_jorb[0]-1])-js[i2+1]))
        while True:
            #main
            m=mtotal2
            j1=jcsp[-1]
            j2=jcsn[-1]
            if abs(j1-j2)<=jj<=j1+j2:
                idx=[]
                for m1 in range(max(-j1,m-j2),min(j1,m+j2)+2,2):
                    m2=m-m1
                    ms=jcsp+[m1]+jcsn+[m2]
                    for i2 in range(sum(n_jorb)-2,n_jorb[0]-1,-1):
                        ms=ms[:i2]+list(jm2mm(js[n_jorb[0]] if i2==n_jorb[0] else ms[i2-1],js[i2+1],ms[i2],ms[i2+1]))+ms[i2+2:]
                    for i2 in range(n_jorb[0]-2,-1,-1):
                        ms=ms[:i2]+list(jm2mm(js[0] if i2==0 else ms[i2-1],js[i2+1],ms[i2],ms[i2+1]))+ms[i2+2:]
                    idx.append(ms_js_rev[js][tuple(ms)])
                wfdata[js][idx,:] = (dcg_matrix(j1,j2,m).T[:,(jj-max(abs(j1-j2),abs(m)))//2]
                                     *numpy.matrix(struct.unpack('<%id'%wfdata[js].shape[1],wavdata[pos:pos+8*wfdata[js].shape[1]])))
                pos += 8*wfdata[js].shape[1]
            #increment
            for i2 in range(sum(n_jorb)-2,n_jorb[0]-1,-1):
                if jcsn[i2-n_jorb[0]]==(js[n_jorb[0]] if i2==n_jorb[0] else jcsn[i2-n_jorb[0]-1])+js[i2+1]: continue
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
               
        for i in range(sum(n_jorb)-2,n_jorb[0]-1,-1):
            if js[i+1]==0: continue
            if i==n_jorb[0] and js[n_jorb[0]]==0: continue
            #initialize
            jcsp=[]
            for i2 in range(n_jorb[0]-1):
                jcsp.append(abs((js[0] if i2==0 else jcsp[i2-1])-js[i2+1]))
            jcsn=[]
            for i2 in range(n_jorb[0],i):
                jcsn.append(abs((js[n_jorb[0]] if i2==n_jorb[0] else jcsn[i2-n_jorb[0]-1])-js[i2+1]))
            mp=-jcsp[-1]
            msn=[-js[i2] for i2 in range(i+2,sum(n_jorb))]
            while True:
                #main
                m=mtotal2-sum([mp]+msn)
                j1=(js[n_jorb[0]] if i==n_jorb[0] else jcsn[i-n_jorb[0]-1])
                j2=js[i+1]
                if j1!=0 and abs(m)<j1+j2:
                    idx=[]
                    for m1 in range(max(-j1,m-j2),min(j1,m+j2)+2,2):
                        m2=m-m1
                        ms=jcsp+[mp]+jcsn+[m1,m2]+msn
                        for i2 in range(i-1,n_jorb[0]-1,-1):
                            ms=ms[:i2]+list(jm2mm(js[n_jorb[0]] if i2==n_jorb[0] else ms[i2-1],js[i2+1],ms[i2],ms[i2+1]))+ms[i2+2:]
                        for i2 in range(n_jorb[0]-2,-1,-1):
                            ms=ms[:i2]+list(jm2mm(js[0] if i2==0 else ms[i2-1],js[i2+1],ms[i2],ms[i2+1]))+ms[i2+2:]
                        idx.append(ms_js_rev[js][tuple(ms)])
                    wfdata[js][idx,:]=dcg_matrix(j1,j2,m).T*wfdata[js][idx,:]
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
                        
        for i in range(n_jorb[0]-2,-1,-1):
            if js[i+1]==0: continue
            if i==0 and js[0]==0: continue
            #initialize
            jcs=[]
            for i2 in range(i):
                jcs.append(abs((js[0] if i2==0 else jcs[i2-1])-js[i2+1]))
            msp=[-js[i2] for i2 in range(i+2,n_jorb[0])]
            msn=[-js[i2] for i2 in range(n_jorb[0],sum(n_jorb))]
            while True:
                #main
                m=mtotal2-sum(msp+msn)
                j1=(js[0] if i==0 else jcs[i-1])
                j2=js[i+1]
                if j1!=0 and abs(m)<j1+j2:
                    idx=[]
                    for m1 in range(max(-j1,m-j2),min(j1,m+j2)+2,2):
                        m2=m-m1
                        ms=jcs+[m1,m2]+msp+msn
                        for i2 in range(i-1,-1,-1):
                            ms=ms[:i2]+list(jm2mm(js[0] if i2==0 else ms[i2-1],js[i2+1],ms[i2],ms[i2+1]))+ms[i2+2:]
                        idx.append(ms_js_rev[js][tuple(ms)])
                    wfdata[js][idx,:]=dcg_matrix(j1,j2,m).T*wfdata[js][idx,:]
                #increment
                for i2 in range(sum(n_jorb)-1,n_jorb[0]-1,-1):
                    if j1==0: continue
                    if msn[i2-n_jorb[0]]==js[i2]: continue
                    msn[i2-n_jorb[0]]+=2
                    for i3 in range(i2+1,sum(n_jorb)):
                        msn[i3-n_jorb[0]]=-js[i3]
                    break
                else:
                    for i2 in range(n_jorb[0]-1,i+1,-1):
                        if j1==0: continue
                        if msp[i2-i-2]==js[i2]: continue
                        msp[i2-i-2]+=2
                        for i3 in range(i2+1,n_jorb[0]):
                            msp[i3-i-2]=-js[i3]
                        for i3 in range(n_jorb[0],sum(n_jorb)):
                            msn[i3-n_jorb[0]]=-js[i3]
                        break
                    else:
                        for i2 in range(i-1,-1,-1):
                            if jcs[i2]==(js[0] if i2==0 else jcs[i2-1])+js[i2+1]: continue
                            jcs[i2]+=2
                            for i3 in range(i2+1,i):
                                jcs[i3]=abs(jcs[i3-1]-js[i3+1])
                            for i3 in range(i+2,n_jorb[0]):
                                msp[i3-i-2]=-js[i3]
                            for i3 in range(n_jorb[0],sum(n_jorb)):
                                msn[i3-n_jorb[0]]=-js[i3]
                            break
                        else:
                            break
                     
    mslist=[tuple(ms) for ms in itertools.product(*[range(-n*(j+1-n),n*(j+1-n)+2,2) for j,n in zip(jorb,occs)]) if sum(ms)==mtotal2]
    wfdata_new=dict([(ms,numpy.empty([len(mbit_jnm[j][n][m]) for j,n,m in zip(jorb,occs,ms)])) for ms in mslist])
    for ms in wfdata_new:
        for i,jvas in enumerate(itertools.product(*[jva_jnm[j][n][m] for j,n,m in zip(jorb,occs,ms)])):
            js,vs,alphas=zip(*jvas)
            js=tuple(js)
            vas=tuple(vs+alphas)
            wfdata_new[ms].flat[i]=wfdata[js][ms_js_rev[js][ms],vas_js_rev[js][vas]]
    wfdata=wfdata_new
   
    for ms in wfdata:
        for i in range(sum(n_jorb)):
            wfdata[ms] = numpy.einsum( wfdata[ms], list(range(i))+[sum(n_jorb)]+list(range(i+1,sum(n_jorb))),
                                       jvabit_jnm[jorb[i],occs[i],ms[i]],[sum(n_jorb),i], list(range(sum(n_jorb))))
         
    r = b""
    for mtotp in range(min(dim_idp_m),max(dim_idp_m)+2,2):
        mtotn=mtotal2-mtotp
        dimp=dim_idp_m[mtotp]
        dimn=dim_idn_m.get(mtotn,0)
        for idimn in range(dimn):
            for idimp in range(dimp):
                mbits=mbit_idp_m[mtotp][idimp]+mbit_idn_m[mtotn][idimn]
                ms=tuple([m_j_mbit[j][mbit] for j,mbit in zip(jorb,mbits)])
                r+=struct.pack('<d',wfdata[ms][tuple([mbit_jnm_rev[j][n][m][mbit] for j,n,m,mbit in zip(jorb,occs,ms,mbits)])])
            
    return r
   
if __name__ == '__main__':
    n_jorb,n_core,norb,lorb,jorb,itorb=read_snt(sys.argv[1])
    n_ferm,iprty,ptn_p,ptn_n,ptn_pn=read_ptn(sys.argv[2])
    print( "snt file: "+sys.argv[1] )
    print( "ptn file: "+sys.argv[2] )
    print( "wav file: "+sys.argv[3] )
    print( "Z=%3i, N=%3i"%(n_ferm[0]+n_core[0],n_ferm[1]+n_core[1]) )
    print( "parity=%2i"%iprty )
   
    dim_jnm={}
    mbit_jn={}
    m_j_mbit={}
    mbit_jnm={}
    mbit_jnm_rev={}
    for j in range(1,max(jorb)+2,2):
        dim_jnm[j] =[{} for x in range(j+2)]
        mbit_jn[j] =[[] for x in range(j+2)]
        mbit_jnm[j] =[{} for x in range(j+2)]
        mbit_jnm_rev[j] =[{} for x in range(j+2)]
        m_j_mbit[j]=[]
        for mbit in range(2**(j+1)):
            n=0
            m=0
            for i in range(j+1):
                n+=(mbit&2**i)>>i
                m+=((mbit&2**i)>>i)*(2*i-j)
            dim_jnm[j][n][m]=dim_jnm[j][n].get(m,0)+1
            mbit_jn[j][n].append(mbit)
            if not m in mbit_jnm[j][n]:
                mbit_jnm[j][n][m]=[]
                mbit_jnm_rev[j][n][m]={}
            mbit_jnm_rev[j][n][m][mbit]=len(mbit_jnm[j][n][m])
            mbit_jnm[j][n][m].append(mbit)
            m_j_mbit[j].append(m)
         
    dim_idp_mp,dim_idp_m,mbit_idp_m,idx_idp_m_mbit=zip(*map(init_mbit_p,ptn_p))
    dim_idn_mp,dim_idn_m,mbit_idn_m,idx_idn_m_mbit=zip(*map(init_mbit_n,ptn_n))
   
    dim_idpn_mp=[]
    dim_mp={}
    for idp,idn in ptn_pn:
        dim_idpn_mp.append(mp_product(dim_idp_mp[idp],dim_idn_mp[idn]))
        mp_add(dim_mp,dim_idpn_mp[-1])
      
    dim_j_vj=dict([(j,calc_dim_vj(j)) for j in range(1,max(jorb)+2,2)])
    jva_jn     =dict([(j,[[(jj,v,alpha) for (v,jj),dim in sorted(dim_j_vj[j].items(),key=lambda x:x[0][::-1])
                                        if (n+v)%2==0 and v<=min(n,j+1-n)
                                        for alpha in range(dim)]
                          for n in range(j+2)])
                      for j in range(1,max(jorb)+2,2)])
    j_jn       =dict([(j,[sorted(list(set([x[0] for x in jva_jn[j][n]])))
                          for n in range(j+2)])
                      for j in range(1,max(jorb)+2,2)])
    va_jnj     =dict([(j,[dict([(jj,[x[1:] for x in jva_jn[j][n] if x[0]==jj])
                                for jj in j_jn[j][n]])
                          for n in range(j+2)])
                      for j in range(1,max(jorb)+2,2)])
    dim_jnj    =dict([(j,[dict([(jj,len(va_jnj[j][n][jj]))
                                for jj in j_jn[j][n]])
                          for n in range(j+2)])
                      for j in range(1,max(jorb)+2,2)])
    jva_jnm    =dict([(j,[dict([(m,[x for x in jva_jn[j][n] if x[0]>=abs(m)])
                                for m in range(-n*(j+1-n),n*(j+1-n)+2,2)])
                          for n in range(j+2)])
                      for j in range(1,max(jorb)+2,2)])
    jva_jnm_rev=dict([(j,[dict([(m,dict([(y,x) for x,y in enumerate(jva_jnm[j][n][m])]))
                                for m in range(-n*(j+1-n),n*(j+1-n)+2,2)])
                          for n in range(j+2)])
                      for j in range(1,max(jorb)+2,2)])
   
    if not os.path.exists(fn):
        bindir = os.path.dirname( __file__ )
        fn = bindir + '/' + fn
    if not os.path.exists(fn):
        raise "execute seniority_basis_exact.py to make seniority_basis.dat"
    jvabit_jnm = read_basis_all( fn )
   
    wavfile=open(sys.argv[3],'rb')
    neig,mtotal=struct.unpack('<ii',wavfile.read(8))
    eval=struct.unpack('<%id'%neig,wavfile.read(8*neig))
    jj=struct.unpack('<%ii'%neig,wavfile.read(4*neig))
    print( "neig=%3i, M=%3i/2"%(neig,mtotal) )
    print( "  i     JP       E(MeV)" )
    for i in range(neig):
        print( "%3i %3i/2%s %12.5f"%(i+1,jj[i],"- +"[iprty+1],eval[i]) )
    print()
    wavoutfile = open(sys.argv[4],'wb')
    mtotal2=mtotal
    outputlist=[x for x in range(neig) if (not jj[x] in [-1]) and jj[x]>=mtotal2]
    neig2=len(outputlist)
    print( "wav output file: " + sys.argv[4] )
    print( "neig=%3i, M=%3i/2" % (neig2, mtotal2) )
    print( "  i     JP       E(MeV)" )
    for j,i in enumerate(outputlist):
        print( "%3i %3i/2%s %12.5f"%(j+1,jj[i],"- +"[iprty+1],eval[i]) )
    wavoutfile.write(struct.pack( '<ii', neig2, mtotal2 ))
    wavoutfile.write(struct.pack( '<%id'%neig2, *[eval[x] for x in outputlist] ))
    wavoutfile.write(struct.pack( '<%ii'%neig2, *[jj[x]   for x in outputlist] ))
   
    p = multiprocessing.Pool(6)
    print(  "convert wave functions" )
    sys.stdout.flush()
    
    for istate in range(neig):
        print(  '\n w.f. state {0:5d} / {1:5d}'.format( istate+1, neig ) )
        sys.stdout.flush()

        if not istate in outputlist:
            wavfile.seek(8*(dim_mp[jj[istate],iprty]-dim_mp.get((jj[istate]+2,iprty),0)),1)
            continue
        chunksize=100
        for idpn in range(0,len(ptn_pn),chunksize):
            wavdata={}
            for idpn2 in range(idpn,min(idpn+chunksize,len(ptn_pn))):
                wavdata[idpn2]=wavfile.read(8*(dim_idpn_mp[idpn2].get((jj[istate],iprty),0)-dim_idpn_mp[idpn2].get((jj[istate]+2,iprty),0)))
            wavoutdata = b"".join( p.map(wrapper_convert_wf,[(x,
             {"wavdata":wavdata[x], "ptn_pn":ptn_pn, "jorb":jorb,
              "occs":ptn_p[ptn_pn[x][0]]+ptn_n[ptn_pn[x][1]],
              "mtotal2":mtotal2, "mbit_jnm":mbit_jnm,
              "dim_idp_m" : dim_idp_m[ptn_pn[x][0]],
              "dim_idn_m" : dim_idn_m[ptn_pn[x][1]],
              "mbit_idp_m" : mbit_idp_m[ptn_pn[x][0]],
              "mbit_idn_m" : mbit_idn_m[ptn_pn[x][1]],
              "m_j_mbit":m_j_mbit, "mbit_jnm_rev":mbit_jnm_rev, "n_jorb":n_jorb, "jvabit_jnm":jvabit_jnm,
              "j_jn":j_jn, "va_jnj":va_jnj, "dim_jnj":dim_jnj, "jva_jnm":jva_jnm, "jj":jj[istate]}
            ) for x in range(idpn,min(idpn+chunksize,len(ptn_pn)))], chunksize=1 ))
            wavoutfile.write(wavoutdata)

    print( '\n finish.' )
            
