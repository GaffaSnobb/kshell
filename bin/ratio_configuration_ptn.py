#!/usr/bin/env python3
#
# show probabilities of particle occupations of KSHELL wave function 
# 
# usage: ratio_configuration_ptn.py foo.snt bar.ptn baz.wav
#
#  original by Y. Tsunoda 2017/04/20 
#  revised  by N. Shimizu 2017/06/20
#

kwf = 8 # bit for wave function
# kwf = 4 # bit for wave function

edn = '<' # little endian, for intel  
# edn = '>' # big endian, for SPARC

thd = 0.03


import sys, math, struct
from count_dim import read_snt, read_ptn, \
    mp_add, mp_product, mps_product, set_dim_singlej, jj2str, jjp2str


def main( fn_snt, fn_ptn, fn_wav, orbs):

    n_jorb, n_core, norb, lorb, jorb, itorb = read_snt(fn_snt)
    n_ferm, iprty, ptn_p, ptn_n, ptn_pn    = read_ptn(fn_ptn)

    print( "snt file: " + fn_snt )
    print( "ptn file: " + fn_ptn )
    print( "wav file: " + fn_wav )
    print()
    print( "Z=%3i, N=%3i" % (n_ferm[0]+n_core[0], n_ferm[1]+n_core[1]) )
    print( "parity=%2i"   % iprty )

    dim_jnm = set_dim_singlej(jorb)
    # dim_jnm = {}
    # for j in set(jorb):
    #     dim_jnm[j] = [ {} for x in range(j+2) ]
    #     for mbit in range(2**(j+1)):
    #         n, m = 0, 0
    #         for i in range(j+1):
    #             n += (mbit&2**i)>>i
    #             m += ((mbit&2**i)>>i)*(2*i-j)
    #         dim_jnm[j][n][m]=dim_jnm[j][n].get(m, 0)+1

    dim_idp_mp = []
    for ptn in ptn_p:
        mps = []
        for i,n in enumerate(ptn):
            mps.append( dict( ( (m, (-1)**(lorb[i]*n)), dim_jnm[jorb[i]][n][m] )
                              for m in dim_jnm[jorb[i]][n] ) )
        dim_idp_mp.append(mps_product(mps))

    dim_idn_mp = []
    for ptn in ptn_n:
        mps = []
        for i,n in enumerate(ptn):
            mps.append( dict( ( (m, (-1)**(lorb[n_jorb[0]+i]*n) ), 
                                dim_jnm[jorb[n_jorb[0]+i]][n][m] ) 
                              for m in dim_jnm[jorb[n_jorb[0]+i]][n] ) )
        dim_idn_mp.append(mps_product(mps))

    dim_idpn_mp = []
    dim_mp = {}
    for idp,idn in ptn_pn:
        dim_idpn_mp.append(mp_product(dim_idp_mp[idp], dim_idn_mp[idn]))
        mp_add(dim_mp, dim_idpn_mp[-1])

    fp_wav = open(sys.argv[3],'rb')
    neig, mtotal = struct.unpack( edn + 'ii', fp_wav.read(8) )
    e_val        = struct.unpack( edn + '%id'%neig, fp_wav.read(8*neig) )
    jj           = struct.unpack( edn + '%ii'%neig, fp_wav.read(4*neig) )

    print( "n_eigen=%3i, M= %s" % (neig, jj2str(mtotal)) )
    print( " n_eig JP       E(MeV)" )
    for i in range(neig):
        print( "%3i %3i/2%s %12.5f" % ( i+1, jj[i], "- +"[iprty+1], e_val[i] ) )
    print()
    print( "#  occ. for orbits    n_eigen  => " )
    print( "#  " +  " "*max(3*len(jorb)-9, 0) + "J = " \
           + " %8s"*neig % tuple(jjp2str(i, iprty) for i in jj) )
    print( "#" + "%3i"*len(jorb) % tuple(jorb) + \
           "%5i    "*neig % tuple(range(1,neig+1)) )

    ratio = [ [] for x in ptn_pn ]
    for istate in range(neig):
        for idpn in range(len(ptn_pn)):
            ndim = dim_idpn_mp[idpn].get((mtotal,iprty), 0)
            if kwf==8:
                wfdata = struct.unpack( edn + '%id'%ndim,
                                        fp_wav.read(8*ndim) )
            else:
                wfdata = struct.unpack( edn + '%if'%ndim,
                                        fp_wav.read(4*ndim) )
            olp = 0.
            for x in wfdata:
                olp += x*x
            ratio[idpn].append(olp)

    sumr = [ 0.0 for i in range(neig) ]

    for idpn in range(len(ptn_pn)):
        idp, idn = ptn_pn[idpn] 
        occs = ptn_p[idp] + ptn_n[idn]
        print( ( " " + "%3i"*len(occs) + "  " + "%9.6f"*neig ) \
               % ( occs + tuple(ratio[idpn]) ) )

        for i in range(neig): 
            sumr[i] += ratio[idpn][i]

    print()
    print( "# total " + " "*max(3*len(occs)-5, 0) + "%9.6f"*neig % tuple(sumr) )
    print()

    print( "### major components ###" )

    tratio = list(map(list, zip(*ratio)))

    for i in range(neig):
        print( '\n# ------------------------------------------------------#' )
        print( "state  # %5i  Jp= %-6s   E= %10.3f MeV" % \
               ( i+1, jjp2str(jj[i], iprty), e_val[i] ) )
        rn = dict( (r,n) for n,r in enumerate(tratio[i]) if r > thd )
        for r in sorted( rn.keys(), reverse=True ):
            idpn = rn[ r ] 
            idp, idn = ptn_pn[ idpn ] 
            occs = ptn_p[idp] + ptn_n[idn]
            print( ( " " + "%3i"*len(occs) + "  %9.6f" ) % ( occs + (r,) ) )
            

        if not orbs: continue

        oorbs = [ j-1 for j in orbs]
        rt = {}
        for idpn in range(len(ptn_pn)):
            idp, idn = ptn_pn[idpn] 
            occs = ptn_p[idp] + ptn_n[idn]
            nl = tuple( occs[j] for j in oorbs )
            rt[nl] = rt.get(nl, 0.) + ratio[idpn][i]

        print() 
        print( ' orbits : ' + ' %3i'*len(orbs) % tuple(orbs) )
        keys = sorted(rt.keys()) 
        for k in keys:
            print( ( ' ' + '%3i'*len(k)+'  %9.6f' ) % ( k + (rt[k],) ) )

        if len(orbs) != 2: continue
        xlist = sorted( list( set( x for x,y in keys ) ) )
        ylist = sorted( list( set( y for x,y in keys ) ) )
        print() 
        print(  '%7i'*len(ylist) % tuple(ylist) + '      sum' )
        vx = [0.,] * len(ylist)
        for x in xlist:
            vy = tuple( rt.get((x,y), 0.) for y in ylist )
            print( '%3i' % x + ' %6.4f'*len(ylist) % vy  \
                   + '   %6.4f' % sum(vy) )
            vx = [ vxx + vyy for vxx, vyy in zip(vx, vy) ]
        print( 'sum ' + '%6.4f '*len(vx) % tuple(vx) )
        print()
        sys.stdout.flush()
            

    

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print( '\nusage: ratio_configuration_ptn.py foo.snt bar.ptn baz.wav\n' )
        sys.exit(1)

    fn_snt, fn_ptn, fn_wav = sys.argv[1:4]

    orbs = [ int(i) for i in sys.argv[4:] ]
    
    main(fn_snt, fn_ptn, fn_wav, orbs)
 
