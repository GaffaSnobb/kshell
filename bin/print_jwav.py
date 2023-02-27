#!/usr/bin/python
#
# 
#
usage_message = \
'''
  usage: print_jwav.py foo.snt bar.ptn baz.jwav
'''

# threshould to show components
thd = 0.01
thd = 0.0001
# thd = 0.0

import sys, struct, itertools

from convert_m2j import read_snt, read_ptn, orbit_coupling, calc_dim_vj


if __name__ == "__main__":

    if len(sys.argv) < 4: 
        print( usage_message )
        sys.exit()

    n_jorb, n_core, norb, lorb, jorb, itorb = read_snt(sys.argv[1])
    n_ferm, iprty, ptn_p, ptn_n, ptn_pn = read_ptn(sys.argv[2])
    print( "snt file: " + sys.argv[1] )
    print( "ptn file: " + sys.argv[2] )
    print( "wav file: " + sys.argv[3] )
    print( "Z= %3i, N= %3i" % (n_ferm[0]+n_core[0], n_ferm[1]+n_core[1]) )
    print( "parity= %2i" % iprty )

    idx_j_n = {}
    for j in range(1, max(jorb)+2, 2):
        dim_vj = calc_dim_vj(j)
        idx_j_n[j] = [[] for n in range(j+2)]
        for n in range(j+2):
            for v,jj in dim_vj:
                # for v,jj in sorted(dim_vj, key=lambda vjj:(vjj[1],vjj[0])):
                if (n+v)%2==1 or v>min(n,j+1-n): continue
                for alpha in range(dim_vj[v,jj]):
                    idx_j_n[j][n].append( (jj, v, alpha) )

    wavfile = open(sys.argv[3],'rb')
    neig, mtotal = struct.unpack( '<ii', wavfile.read(8) )
    e_val        = struct.unpack( '<%id'%neig, wavfile.read(8*neig) )
    jj           = struct.unpack( '<%ii'%neig, wavfile.read(4*neig) )

    print( "neig=%3i, M=%3i/2" % (neig, mtotal) )
    print( "  i     JP       E(MeV)" )
    for i in range(neig):
        print( "%3i %3i/2%s %12.5f" % (i+1, jj[i], "- +"[iprty+1], e_val[i]) )

    order = orbit_coupling( n_jorb, norb, lorb, jorb, itorb )

    print()

    for istate in range(neig):
        print( "----------------------------------------" )
        print()
        print( "  i     JP       E(MeV)" )
        print( "%3i %3i/2%s %12.5f" % \
            (istate+1, jj[istate], "- +"[iprty+1], e_val[istate]) )
        print()
        sys.stdout.flush()

        rto_vs = {}
        value_idx = []

        for idp,idn in ptn_pn:

            occs = ptn_p[idp] + ptn_n[idn]

            idxs = [ tuple(zip(*x))
                     for x in itertools.product(
                             *[ idx_j_n[jorb[i]][occs[i]] 
                                for i in range(sum(n_jorb)) ] ) ]

            for orb1,orb2 in order:
                idxs_new = []
                for idx in idxs:
                    j1, j2 = idx[0][orb1], idx[0][orb2]
                    for j3 in range( abs(j1-j2), j1+j2+2, 2):
                        if (orb1,orb2) == order[-1] and j3 != jj[istate]:
                            continue
                        idxs_new.append( (idx[0]+(j3,),)+idx[1:] )
                idxs = idxs_new

            # ( (J1,J2,..,J12,..,Jc), (v1, v2, ...,vj), (a1,a2, ...aj) )
            idxs = sorted(idxs)

            for idx in idxs:
                value = struct.unpack( '<d', wavfile.read(8) )[0]
                vs = sum( idx[1] )
                rto_vs[vs]  = rto_vs.get(vs, 0.) + value**2
                
                if value**2 < thd: continue

                value_idx.append( (value, idx, occs) )

        value_idx = sorted( value_idx, key=lambda x : x[0]**2,
                            reverse=True )
        for value, idx, occs in value_idx:

            print( "index: " + ("%2i " * sum(n_jorb)) % tuple( 
                range(1, sum(n_jorb)+1) ) )
            print( "nl   : " + ("%2s " * sum(n_jorb)) % tuple(
                str(x) + "spdfghijk"[y] for x,y in zip(norb, lorb) ) )
            print( "2j   : " + ("%2i " * sum(n_jorb)) % tuple(jorb) )
            print( "n    : " + ("%2i " * sum(n_jorb)) % occs )
            print( "2J   : " + ("%2i " * sum(n_jorb)) % idx[0][:sum(n_jorb)] )
            print( "Jcoup:    " + ( "%2i "*(n_jorb[0]-1) +
                                    "   " + "%2i "*n_jorb[1]) \
                   % idx[0][sum(n_jorb):] )
            print( "v    : " + ("%2i " * sum(n_jorb)) % idx[1] )
            print( "alpha: " + ("%2i " * sum(n_jorb)) % tuple(
                x+1 for x in idx[2] ) )
            print( "value: %10.7f" % value )
            print( "v**2 : %10.7f" % value**2 )
            print()


        print( "*** seniority ratio *** " )
        print( "  v :  ratio " )
        for v in sorted( rto_vs.keys() ):
            print( "%3i : %10.7f" % ( v, rto_vs[v] ) )
        print()

