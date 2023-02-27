#!/usr/bin/env python
#
# ./add_snt.py foo.snt bar.snt output.snt
# 
#   output = foo + bar
#
#    by N. Shimizu 2019/07/09
# 

import sys, os

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
            i, j = int(arr[0])-1, int(arr[1])-1
            if i>j: i,j = j,i
            vob[(i,j)] = float(arr[2])
        self.vob = vob


        nline = read_comment_skip(fp)
        if nline[1] == 0: 
            # print "# no mass dependece"
            self.fmd_mass = 0
        elif nline[1] == 1:
            self.fmd_mass, self.fmd_power = float(nline[2]), float(nline[3])
            # print '# mass dependence (A/', self.fmd_mass, ')^', self.fmd_power
        else:
            raise "mass dependence not supported"
        
        vtb = {}
        for i in range(nline[0]):
            arr = read_comment_skip(fp)
            i, j, k, l = [ int(arr[i])-1 for i in range(4) ]
            J = int(  arr[4])
            v = float(arr[5])

            if i > j: 
                v *= (-1) **( (self.jorb[i]+self.jorb[j])/2-J + 1)
                i, j = j, i
            if k > l:
                v *= (-1) **( (self.jorb[k]+self.jorb[l])/2-J + 1)
                k, l = l, k
            if i > k or (i == k and j > l):
                i, j, k, l = k, l, i, j

            if (i,j,k,l,J) in vtb:
                if abs( v - vtb[(i,j,k,l,J)] )>1.e-4:
                        print 'WARNING duplicate TBME', i,j,k,l,J,v,vtb[(i,j,k,l,J)]
                        continue
            vtb[(i,j,k,l,J)] = v
        fp.close()

        self.tbme = vtb


    def fmd(self, mass):
        # mass dependence
        if self.fmd_mass == 0: return 1.0
        return (mass*1.0/self.fmd_mass)**self.fmd_power

        

    def add(self, smr):
        for i in range(len(self.norb)):
            if    self.norb[i] != smr.norb[i] or self.lorb[i] != smr.lorb[i] \
               or self.jorb[i] != smr.jorb[i] or self.itorb[i] != smr.itorb[i]:
                raise ' inconsistent model space '
            pass
        
        for ij in set( self.vob.keys() + smr.vob.keys() ):
            self.vob[ij] = self.vob.get( ij, 0. ) + smr.vob.get( ij, 0. )

        for ijklJ in set( self.tbme.keys() + smr.tbme.keys() ):
            self.tbme[ijklJ] = self.tbme.get( ijklJ, 0. ) + smr.tbme.get( ijklJ, 0. ) 

        return



    def print_snt(self, comment=None):
        out = '!\n'
        if comment: out += '! ' + comment
        out += '!\n'
        out += '! model space \n'
        out += ' %3d %3d   %3d %3d\n' % tuple(self.n_jorb + self.n_core)
        for i, (n, l, j, tz) in enumerate( zip( self.norb, self.lorb, 
                                                self.jorb, self.itorb ) ):
            out += "%5d   %3d %3d %3d %3d  !  %2d = %s \n" \
                % ( i+1, n, l, j, tz, i+1, orb2char( n, l, j, tz ) )
        # one-body int.
        out += "! interaction\n"
        out += "! num, method=,  hbar_omega\n"
        out += "!  i  j     <i|H(1b)|j>\n"
        out += " %3d %3d\n" % (len(self.vob), 0)
        for i,j in sorted(self.vob):
            out += "%3d %3d % 15.8f\n" % (i+1, j+1, self.vob[(i,j)])
        # two-body int.
        out += "! TBME\n"
        if self.fmd_mass == 0:
            out += " %5d  %3d\n" % (len(self.tbme), 0)
        else:
            out += " %5d  %3d %3d % 15.8f\n" % (len(self.tbme), 1, self.fmd_mass, self.fmd_power )
            pass
        idx = sorted( self.tbme.keys() )
        for Tz in (-2, 0, 2): 
            for i, j, k, l, J in idx:
                v = self.tbme[ (i,j,k,l,J) ]
                if self.itorb[i] + self.itorb[j] != Tz: continue
                out += "%3d %3d %3d %3d  %3d   % 15.8f\n" % (i+1, j+1, k+1, l+1, J, v)

        return out

    def vob_get(self, i, j):
        return self.vob.get0( i-1, j-1 )

    def vob_get0(self, i, j):
        return self.vob.get( (i,j), 0.0 )

    def tbme_get(self, i, j, k, l, J):
        i, j, k, l = i-1, j-1, k-1, l-1
        return self.tbme_get0( i, j, k, l, J )

    def tbme_get0(self, i, j, k, l, J):
        v = 1.
        if i > j: 
            v *= (-1) **( (self.jorb[i]+self.jorb[j])/2-J + 1)
            i, j = j, i
        if k > l:
            v *= (-1) **( (self.jorb[k]+self.jorb[l])/2-J + 1)
            k, l = l, k
        if i > k or (i == k and j > l):
            i, j, k, l = k, l, i, j
        return v * self.tbme.get( (i,j,k,l,J), 0.0 )
        


if __name__ == "__main__":
    if len(sys.argv) < 3: 
        print 'usage: add_snt.py foo.snt bar.snt (foobar.snt)'
        sys.exit(1)

    sm  = SMInt( sys.argv[1] )
    smr = SMInt( sys.argv[2] )
    sm.add(smr)

    comment = ' add_snt.py ' + sys.argv[1] + ' ' + sys.argv[2] + '\n'
        
    out = sm.print_snt( comment )
    
    fp = sys.stdout if len(sys.argv) < 4 else open(sys.argv[3], 'w')
    fp.write(out) 

    
