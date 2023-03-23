#!/usr/bin/python3
#
# ./espe.py foobar.snt #proton #neutron
#     e.g.   ./espe.py w.snt 4 0-8 
#
#  monopole int. and  effective SPE by snt file
#
#    by N. Shimizu 2014/05/24
# 

import sys
import numpy as np
from kshell_ui import element


# order of filling occupation (n, l, j)
order_orb = ( (0, 0, 1),  # s-shell
              (0, 1, 3),  (0, 1, 1),  # p-shell
              (0, 2, 5),  (1, 0, 1),  (0, 2, 3),  # sd-shell
              (0, 3, 7),  (1, 1, 3),  (1, 1, 1),  (0, 3, 5), #pf-shell
              (0, 4, 9),  # g9/2
              (0, 4, 7),  (1, 2, 5),  (1, 2, 3), 
              (2, 0, 1),  (0, 5, 11), # 50-82
              (0, 5, 9),  (1, 3, 7),  (1, 3, 5),  (2, 1, 3), 
              (2, 1, 1),  (0, 6, 13), # 82-126
              (1, 4, 9),  (2, 2, 5),  (0, 6, 11), (1, 4, 7), 
              (3, 0, 1),  (2, 2, 3),  (0, 7, 15)  # 126-184
          )


def orb2char(n,l,j,tz):
    lorb2c = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 
              'm', 'n', 'o']
    tz2c = { -1:'p', 1:'n' }
    return "%c%d%c%02d/2" % (tz2c[tz], n, lorb2c[l], j)

def read_comment_skip(fp):
    while True:
        arr = fp.readline().split()
        if not arr: return None
        for i in range(len(arr)): 
            if arr[i][0]=="!" or arr[i][0]=="#": 
                arr = arr[:i]
                break
        if not arr: continue
        try:
            return [int(i) for i in arr]
        except ValueError:
            try:
                return [int(i) for i in arr[:-1]]+[float(arr[-1])]
            except ValueError:
                return arr


class SMInt:
    def __init__(self, fn_snt):
        fp = open(fn_snt, 'r')
        self.n_jorb, self.n_core = [0,0], [0,0]
        self.n_jorb[0], self.n_jorb[1], self.n_core[0], self.n_core[1] \
            = read_comment_skip(fp)
        self.n_jorb_pn = sum(self.n_jorb)
        self.norb, self.lorb, self.jorb, self.itorb = [], [], [], []
        for i in range(sum(self.n_jorb)):
            arr = read_comment_skip(fp)
            if i+1!=arr[0]: raise "read error"
            self.norb.append(arr[1])
            self.lorb.append(arr[2])
            self.jorb.append(arr[3])
            self.itorb.append(arr[4])
            if (i< self.n_jorb[0] and arr[4]!=-1) or \
               (i>=self.n_jorb[0] and arr[4]!= 1): 
                raise "read error"

        self.lchar = [ orb2char(n,l,j,tz) for n,l,j,tz 
                  in zip(self.norb, self.lorb, self.jorb, self.itorb) ]

        vob = {}
        nline = read_comment_skip(fp)
        for i in range(nline[0]):
            arr = read_comment_skip(fp)
            i, j = int(arr[0]), int(arr[1])
            i, j = i-1, j-1
            vob[(i,j)] = float(arr[2])
            if i!=j: vob[(j,i)] = vob[(i,j)]


        vtb = {}
        nline = read_comment_skip(fp)
        if nline[1] == 0: 
            print( "# no mass dependece" )
            self.fmd_mass = 0
        elif nline[1] == 1:
            self.fmd_mass, self.fmd_power = float(nline[2]), float(nline[3])
            print( '# mass dependence (A/', self.fmd_mass, ')^', self.fmd_power )
        else:
            raise "mass dependence not supported"


        for i in range(nline[0]):
            arr = read_comment_skip(fp)
            i, j, k, l, J  = [ int(arr[i]) for i in range(5) ]
            i, j, k, l = i-1, j-1, k-1, l-1
            v = float(arr[5])
            if (i,j,k,l,J) in vtb:
                if abs( v - vtb[i,j,k,l,J] )>1.e-3:
                        print( 'WARNING duplicate TBME', i+1,j+1,k+1,l+1,J,v,vtb[(i,j,k,l,J)] )
            vtb[(i,j,k,l,J)] = v
            sij = (-1) **( (self.jorb[i]+self.jorb[j])/2-J + 1)
            skl = (-1) **( (self.jorb[k]+self.jorb[l])/2-J + 1)
            if i!=j:          vtb[(j,i,k,l,J)] = v*sij
            if k!=l:          vtb[(i,j,l,k,J)] = v*skl
            if i!=j and k!=l: vtb[(j,i,l,k,J)] = v*sij*skl
            if (i,j)!=(k,l):
                vtb[(k,l,i,j,J)] = v
                if i!=j:          vtb[(k,l,j,i,J)] = v*sij
                if k!=l:          vtb[(l,k,i,j,J)] = v*skl
                if i!=j and k!=l: vtb[(l,k,j,i,J)] = v*sij*skl
        fp.close()

        nj = sum(self.n_jorb)
        self.spe = np.zeros( (nj) )
        for i,j in vob.keys(): 
            if i != j: continue
            self.spe[i] = vob[(i,i)]

        self.sort_orbits()

        self.vm = np.zeros( (nj, nj) )
        # non-diagonal
        for i in range(nj):
            for j in range(nj):
                Jmin = abs(self.jorb[i]-self.jorb[j]) // 2
                Jmax = (self.jorb[i]+self.jorb[j]) // 2
                skip = 1
                if i==j: skip = 2
                for J in range(Jmin, Jmax+1, skip):
                    if not (i,j,i,j,J) in vtb.keys(): 
                        print( "# Warning not found TBME entry",i+1,j+1,i+1,j+1,J )
                v = sum([ (2.*J+1.)*vtb.get((i,j,i,j,J), 0.0)
                          for J in range(Jmin, Jmax+1, skip) ]) 
                d = sum([ 2.*J+1.
                          for J in range(Jmin, Jmax+1, skip) ]) 
                self.vm[i,j] = v/d

        for t1,t2,ct in ((-1,-1,'p-p'), (1,1,'n-n'), (-1,1,'p-n')):
            print( '# ---- monopole int. type ', ct )
            print( '#           ', end='' )
            for j in range(nj): 
                if self.itorb[j]==t2: print( self.lchar[j]+' ', end='' )
            print()
            for i in range(nj):
                if self.itorb[i]!=t1: continue
                print( '#%2d %s' % (i+1, self.lchar[i]), end='' )
                for j in range(nj):
                    if self.itorb[j]!=t2: continue
                    print( '%8.3f' % (self.vm[i,j],), end='' )
                print()

    def fmd(self, mass):
        # mass dependence
        if self.fmd_mass == 0: return 1.0
        return (mass*1.0/self.fmd_mass)**self.fmd_power

    def energy_occ(self, occ):
        f = self.fmd( sum(occ) + sum(self.n_core) )
        e = sum( np.asarray(occ) * self.spe )
        for i,ni in enumerate(occ):
            for j,nj in enumerate(occ):
                if i<j: 
                    continue
                elif i==j:
                    e += occ[i]*(occ[i]-1)*0.5 * self.vm[i,i] * f
                else:
                    e += occ[i]*occ[j] * self.vm[i,j] * f
        return e

    def sort_orbits(self):
        nj = sum(self.n_jorb)
        sped = np.asarray( [ ( order_orb.index( 
            (self.norb[i], self.lorb[i], self.jorb[i]) ), 
                               self.itorb[i], i ) for i in range(nj) ] )
        # occupy order ... s.p.e.
        # sped = np.asarray( [ (self.spe[i], self.itorb[i], i) 
        #                      for i in range(nj) ] )
        sped = sped[ sped[:,0].argsort() ]
        self.orbsort = [ [int(sped[i,2]) for i in range(nj) if sped[i,1]==-1], 
                         [int(sped[i,2]) for i in range(nj) if sped[i,1]== 1] ]
        print( '# filling order proton  orbit ', \
            [ i+1 for i in self.orbsort[0] ] )
        print( '# filling order proton  orbit ', \
            [ self.lchar[i] for i in self.orbsort[0] ] )
        print( '# filling order neutron orbit ', \
            [ i+1 for i in self.orbsort[1] ] )
        print( '# filling order neutron orbit ', \
            [ self.lchar[i] for i in self.orbsort[1] ] )
        



def cal_espe(sm, nf):
    mass = sum(sm.n_core) + sum(nf)
    nj = sum(sm.n_jorb)

    # occupation
    nfull = np.asarray(sm.jorb)+1
    nocc = np.zeros( (nj,), dtype=int )
    for t in range(2):
        nr = nf[t]
        for i in sm.orbsort[t]:
            if nr > nfull[i]:
                nocc[i] = nfull[i]
                nr -= nfull[i]
            else:
                nocc[i] = nr
                nr -= nr
                break
        if nr!=0: raise 'particle number error'

    # print "# nocc",nocc
    e_core = sm.energy_occ(nocc)
    # print "# core energy ",e_core
    
    espe = []
    for i in range(nj):
        occ = nocc.copy()
        if nfull[i] == nocc[i]:
            occ[i] -= 1
            e = e_core - sm.energy_occ(occ)
        else:
            occ[i] += 1
            e = sm.energy_occ(occ) - e_core
        espe.append(e)

    # espe = []
    # for i in range(nj):
    #     e = sm.spe[i]
    #     if nfull[i] == nocc[i]:
    #         e += sm.vm[i,i]*sm.fmd(mass) * (nocc[i]-1)
    #     else:
    #         e += sm.vm[i,i]*sm.fmd(mass) * nocc[i]
    #     for j in range(nj):
    #         if j==i: continue
    #         e += sm.vm[i,j]*sm.fmd(mass) * nocc[j] 
    #         # print " %3d %10.5f # %s" % (i+1, e, sm.lchar[i])
    #     espe.append(e)
    
    return espe



def read_log_occ():


    return occ
        
    


if __name__ == "__main__":
    if len(sys.argv)==1: 
        print( 'usage: espe.py hoge.snt #proton #neutron ' )
        print( '  # of valence particles with range' )
        print( '  e.g.   ./espe.py w.snt 4 0-8  ' )
        sys.exit(1)

    fn_snt = sys.argv[1]
    
    if len(sys.argv)>=3:
        nf1 = [ int(i) for i in sys.argv[2].split('-') ]
    else:
        nf1 = [ 0, ]
    if len(sys.argv)>=4:
        nf2 = [ int(i) for i in sys.argv[3].split('-') ]
    else: 
        nf2 = [ 0, ]

    if len(nf1)==1: nf1 = [ nf1[0], nf1[0] ]
    if len(nf2)==1: nf2 = [ nf2[0], nf2[0] ]
    nflist = [ (i,j) for i in range(nf1[0], nf1[1]+1) 
               for j in range(nf2[0], nf2[1]+1)  ]
    
    sm = SMInt(fn_snt)

    print( "# Z  N   ", end='' )
    for i in range(sm.n_jorb_pn):
        print(  sm.lchar[i] + "  ", end='' )
    print()
    for nf in nflist:
        espe = cal_espe(sm, nf)
        print( "%3d %3d" % (nf[0]+sm.n_core[0], nf[1]+sm.n_core[1]), end='' )
        for e in espe:
            print( " %8.3f" % e, end='' )
        print()

