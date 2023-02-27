#!/usr/bin/env python
#
# generate "seniority_basis.dat", which contains j-scheme basis of single orbit
#
#  usage: ./seniority_basis_exact.py
#


# parameters 

maxj = 13 # maximum j is 13/2 for the model space
# maxj = 15
# maxj = 9

fb = '<' # little endian for intel CPU 
# fb = '>' # big endian for SPARC CPU

fname = 'seniority_basis.dat'



# save_memory = True
save_memory = False

if not save_memory:
    isqrt_cache={}


import sympy, fractions, sys, struct
from convert_m2j import calc_dim_vj
   
def isqrt(n):
    if not save_memory:
        if n in isqrt_cache:
            return isqrt_cache[n]
    x = n
    if x<=0: return 0
    a,b = divmod(x.bit_length(), 2)
    m = 2**(a+b)
    while True:
        y = (m+x/m)>>1
        if y >= m:
            if not save_memory:
                isqrt_cache[n] = m
            return m
        m = y
      
def fraction_sqrt(x):
    if x<0:
        return -fractions.Fraction(isqrt(-x.numerator), isqrt(x.denominator))
    else:
        return  fractions.Fraction(isqrt( x.numerator), isqrt(x.denominator))
      
class Sqrt_fraction:
    # sgn(fraction)*sqrt(abs(fraction))
    def __init__(self, fraction):
        self.fraction = fractions.Fraction(fraction)
        
    def __str__(self):
        if self.fraction<0:
            return "-sqrt(%s)" % str(-self.fraction)
        else:
            return  "sqrt(%s)" % str( self.fraction)
        
    def __lt__(self,other):
        if other:
            return NotImplemented
        else:
            return self.fraction < 0
        
    def __nonzero__(self):
        return bool(self.fraction)
    
    def __add__(self, other):
        if not self.fraction:
            return other
        tmp = 1 + fraction_sqrt(other.fraction/self.fraction)
        if tmp<0:
            tmp = Sqrt_fraction(-tmp**2)
        else:
            tmp = Sqrt_fraction( tmp**2)
        return self * tmp
    
    def __sub__(self, other):
        return self+-other
    
    def __mul__(self, other):
        return Sqrt_fraction(self.fraction*other.fraction)
    
    def __div__(self,other):
        return Sqrt_fraction(self.fraction/other.fraction)
    
    def __radd__(self,other):
        return self+other
    
    def __rsub__(self,other):
        return -(self-other)
    
    def __rmul__(self,other):
        return self*other
    
    def __rdiv__(self,other):
        return Sqrt_fraction(other.fraction/self.fraction)
    
    def __neg__(self):
        return Sqrt_fraction(-self.fraction)
    
    def sqrt(self):
        return Sqrt_fraction(fraction_sqrt(self.fraction))
    
    def __float__(self):
        if self.fraction<0:
            return -(-self.fraction)**.5
        else:
            return self.fraction**.5
         
sf0 = Sqrt_fraction(0)
sf1 = Sqrt_fraction(1)

def clebsch_gordan_2j(j1, j2, j3, m1, m2, m3):
    import sympy.physics.wigner
    s2 = sympy.S(2)
    tmp = sympy.physics.wigner.clebsch_gordan(j1/s2, j2/s2, j3/s2, m1/s2, m2/s2, m3/s2)
    if tmp < 0:
        return Sqrt_fraction( "-"+str(tmp**2) )
    else:
        return Sqrt_fraction(     str(tmp**2) )
      
   
# bitwf
# Ex. {0b0011:2,0b0101:3} represents 2|0011>+3|0101> or 2 c+_0010 c+_0001 + 3 c+_0100 c+_0001

def strip_bitwf(c):
    # remove mbit with 0 coefficient in |c> 
    for key in list(c):
        if not c[key]: del c[key]
    return c 

def bitwf_add(a, b):
    c = dict(a)
    for key in b:
        c[key] = c.get(key, sf0) + b[key]
    return strip_bitwf(c)
   
def bitwf_times_s(a, b):
    # a:scalar, b:bitwf
    # return a * |b>
    if a == 0: return {}
    return { key: a*b[key] for key in b }
   
bit_sign_cache = {}

def bit_sign(a, b):
    # sign for bit operation, a*|M_b>
    if (a,b) in bit_sign_cache:
        return bit_sign_cache[a,b]
    r = 1
    c = a
    while c:
        r *= (-1)**bin( (~c^(c-1)) & b ).count("1")
        c = c & (c-1)  # remove the last '1' bit 
    bit_sign_cache[a,b] = r
    return r
   
def bitwf_times_c(a, b):
    # a:creation operators, b:bitwf
    # return |c> = a*|b>
    c = {}
    for key1 in a:
        for key2 in b:
            if key1 & key2: continue
            key = key1 | key2
            c[key] = c.get(key, sf0) \
                     + Sqrt_fraction(bit_sign(key1, key2)) * a[key1] * b[key2]
    return strip_bitwf(c)
   
def bitwf_times_a(a, b):
    # a: creation operators, b: bitwf
    # return dagger(a)*|b>
    c={}
    for key1 in a:
        s = 1 if bin(key1).count("1")%4 in [0, 3] else -1
        for key2 in b:
            if key1 & ~key2: continue
            key = ~key1 & key2
            c[key] = c.get(key, sf0) \
                     + Sqrt_fraction( bit_sign(key1, key2) * s) * a[key1] * b[key2]
    return strip_bitwf(c)
   
def bitwf_dot(a, b):
    # return <a | b>
    r = sf0
    for key in a:
        if not (key in b): continue
        r += a[key] * b[key]
    return r
   
def bitwf_jdown(a, j):
    # return [J-,a]
    c = {}
    for key in a:
        for i in range(j):
            if (key & (3<<i))>>i == 2:
                c[key-(1<<i)] = c.get(key-(1<<i), sf0) \
                                + Sqrt_fraction((j-i)*(i+1)) * a[key]
    return strip_bitwf(c)
   
def bitwf_jup(a, j):
    # return [J+,a]
    c={}
    for key in a:
        for i in range(j):
            if (key & (3<<i))>>i == 1:
                c[key+(1<<i)] = c.get(key+(1<<i), sf0) \
                                + Sqrt_fraction((j-i)*(i+1)) * a[key]
    return strip_bitwf(c)
   
def bitwf_jz(a,j):
    # return [Jz,a]
    c = {}
    for key in a:
        m = sum( i*2-j for i in range(j+1) if key & (1<<i) )
        x = fractions.Fraction(m, 2)**2
        c[key] = Sqrt_fraction( -x if m<0 else x ) * a[key]
    return strip_bitwf(c)
   
def seniority_up(j):
    # return S+
    return dict( ((1<<x)+(1<<(j-x)), Sqrt_fraction((-1)**x))
                 for x in range((j+1)/2) )
   
def bitwf_szero(a, j):
    # return [S0, a]
    c = {}
    for key in a:
        s = bin(key).count("1") - (j+1)/2
        x = fractions.Fraction(s,2)**2
        c[key] = Sqrt_fraction( -x if s<0 else x ) * a[key]
    return strip_bitwf(c)
   
def str_bitwf(a, j):
    r = ""
    for key in sorted(a):
        r += ( "{:0>"+str(j+1)+"b}: {}\n" ).format( key, str(a[key]) )
    return r
   
def bitwf_ortho(a, b):
    # orthogonalize |a> with b=(b1,b2,...)
    if not isinstance(b, list):
        # |a> - |b> * <b|a>
        return bitwf_add( a, bitwf_times_s(-bitwf_dot(a,b), b) )
    for bitwf in b:
        a = bitwf_ortho(a, bitwf)
    return a
      
def gs_bitwf(bitwfs):
    # Gram Schmidt orthogonalization
    bitwfs = [x for x in bitwfs if x!={}]
    if not bitwfs: return []

    # pick up v(i) and fix its phase
    bit = min([x for bitwf in bitwfs for x in bitwf])
    for i in range(len(bitwfs)):
        if bit in bitwfs[i]: break
    a = bitwfs[i]
    if a[bit] < 0:
        a = bitwf_times_s(-sf1, a)

    # |b> = |b> - b(bit)/a(bit) * |a>  ?
    b = bitwfs[:i] + bitwfs[i+1:]
    for i in range(len(b)):
        if bit in b[i]:
            b[i] = bitwf_add(b[i], bitwf_times_s(-b[i][bit]/a[bit], a))
            
    b = gs_bitwf(b) # recursive call
    a = bitwf_ortho(a, b)
    a = bitwf_times_s(sf1/bitwf_dot(a,a).sqrt(), a)  # normalize
    return [a]+b
   
def pair_creation(jpair, j):
    # pair creation operator with J=jpair in j/2 orbit
    op = {}
    for m1 in range(-j, j+2, 2):
        for m2 in range( max(m1+2,-m1-jpair*2), min(j,-m1+jpair*2)+2, 2 ):
            op[ (1<<((m1+j)/2)) + (1<<((m2+j)/2)) ] \
                = Sqrt_fraction(2) * \
                clebsch_gordan_2j(j, j, jpair*2, m2, m1, m1+m2)
    return op
   
def bitwf_ortho_j9(bitwfs, bitwfv2):
    # choose solvable state as alpha=1
    # only for j=9/2, v=4, J=4,6
    op = pair_creation(2, 9)
    bitwf = bitwf_times_a(op, bitwfv2) 
    a = bitwf_add(bitwf_times_s( bitwf_dot(bitwf, bitwf_times_a(op, bitwfs[1])), bitwfs[0]),
                  bitwf_times_s(-bitwf_dot(bitwf, bitwf_times_a(op, bitwfs[0])), bitwfs[1]))
    a = gs_bitwf([a])[0]
    b = gs_bitwf([bitwf_ortho(bitwfs[0],a)])[0]
    return [a,b]
   
def calc_seniority_basis(j):
    # return bitwf_vj[n,v,alpha,2J,2M], dim_vj[v,2J]
    dim_vj = calc_dim_vj(j)
    bitwf_vj = {}
    sup = seniority_up(j)
    for v in range((j+1)/2+1):
        if v==0:
            bitwf_vj[0,0,0,0,0] = {0:sf1}
        else:
            for j1 in [key[1] for key in sorted(dim_vj) if key[0]==v]:
                print>>sys.stderr, "v=", v, "2j=", j1, "dim=", dim_vj[v,j1]
                bitwfs = []
                for j2,alpha in sorted([
                        (key[1], x)
                        for key in dim_vj if key[0]==v-1 and abs(j1-j)<=key[1]<=j1+j
                        for x in range(dim_vj[key]) ]):
                    bitwf = {}
                    for m in range(max(j1-j2,-j), min(j1+j2,j)+2, 2):
                        if save_memory:
                            if not (v-1,v-1,alpha,j2,j1-m) in bitwf_vj:
                                bitwf_vj[v-1,v-1,alpha,j2,j1-m] \
                                    = bitwf_times_s(
                                        sf1/Sqrt_fraction((j2+j1-m+2)*(j2-j1+m)/4),
                                        bitwf_jdown( bitwf_vj[v-1, v-1, alpha, j2, j1-m+2], j ) )
                        bitwf = bitwf_add( bitwf,bitwf_times_s(
                            clebsch_gordan_2j( j, j2, j1, m, j1-m, j1 ),
                            bitwf_times_c( {1<<((j+m)/2):sf1},
                                           bitwf_vj[v-1,v-1,alpha,j2,j1-m]) ))
                    bitwf = bitwf_ortho( bitwf,
                                         [ bitwf_vj[key] for key in bitwf_vj
                                           if key[0]==v and key[1]==v-2 and key[3]==key[4]==j1 ] )
                    if bitwf:
                        bitwfs.append( bitwf )
                        if len(bitwfs) == dim_vj[v,j1]:
                            bitwfs = gs_bitwf(bitwfs)
                        if len(bitwfs) == dim_vj[v,j1]:
                            break
                if j==9 and v==4 and j1 in [8,12]:
                    bitwfs = bitwf_ortho_j9(bitwfs, bitwf_vj[v,v-2,0,j1,j1])
                for i in range(dim_vj[v,j1]):
                    bitwf_vj[v,v,i,j1,j1] = bitwfs[i]
        if save_memory:
            for key in bitwf_vj.keys():
                if (key[0]!=key[1] and key[0]<=v) or (key[3]!=key[4] and key[0]<v):
                    del bitwf_vj[key]
        for j1 in [key[1] for key in sorted(dim_vj) if key[0]==v]:
            for alpha in range(dim_vj[v,j1]):
                print>>sys.stderr, "v=",v,"2j=",j1,"dim=",alpha+1,"/",dim_vj[v,j1]
                if save_memory:
                    if v+2<=(j+1)/2:
                        bitwf_vj[v+2,v,alpha,j1,j1] \
                            = bitwf_times_s( sf1/Sqrt_fraction((j+1-v*2)/2),
                                             bitwf_times_c(sup,bitwf_vj[v,v,alpha,j1,j1]) )
                else:
                    for m1 in range(j1,-j1-2,-2):
                        if m1 != j1:
                            bitwf_vj[v,v,alpha,j1,m1] \
                                = bitwf_times_s( sf1/Sqrt_fraction((j1+m1+2)*(j1-m1)/4),
                                                 bitwf_jdown(bitwf_vj[v,v,alpha,j1,m1+2],j) )
                        for n in range(v+2,j-v+3,2):
                            bitwf_vj[n,v,alpha,j1,m1] \
                                = bitwf_times_s( sf1/Sqrt_fraction((n-v)*(j-n-v+3)/4),
                                                 bitwf_times_c(sup, bitwf_vj[n-2,v,alpha,j1,m1]) )
    return bitwf_vj, dim_vj
   
def check_seniority_basis(bitwf_vj, dim_vj, j):
    sqrt_1_over_4 = Sqrt_fraction(fractions.Fraction(1,4)),
    sup = seniority_up(j)
    for key in sorted(bitwf_vj, key=lambda x:(x[0],x[1],x[3],x[2],-x[4])):
        jjwf = bitwf_add( bitwf_times_s(
            sqrt_1_over_4, 
            bitwf_add( bitwf_jdown(bitwf_jup(bitwf_vj[key], j),j),
                       bitwf_jup(bitwf_jdown(bitwf_vj[key], j),j) ) ),
                          bitwf_jz(bitwf_jz(bitwf_vj[key],j), j) )
        sswf = bitwf_add( bitwf_times_s(
            sqrt_1_over_4, 
            bitwf_add( bitwf_times_a(sup, bitwf_times_c(sup,bitwf_vj[key])),
                       bitwf_times_c(sup, bitwf_times_a(sup,bitwf_vj[key])) ) ), 
                          bitwf_szero( bitwf_szero(bitwf_vj[key],j), j) )
        jj = bitwf_dot(jjwf, bitwf_vj[key])
        ss = bitwf_dot(sswf, bitwf_vj[key])
        print "n=%2i v=%1i alpha=%2i/%2i j=%4s m=%5s: jj=%s ss=%s" % \
            ( key[0], key[1], key[2]+1, dim_vj[key[1],key[3]],
              "%2i"%(key[3]/2) if key[3]%2==0 else "%2i/2"%key[3],
              "%3i"%(key[4]/2) if key[4]%2==0 else "%3i/2"%key[4],
              str(jj.sqrt().fraction), str(ss.sqrt().fraction) )
        for alpha in range(key[2]):
            print "overlap with alpha=%2i: %s" % \
                ( alpha + 1,
                  bitwf_dot(bitwf_vj[key], bitwf_vj[key[:2]+(alpha,)+key[3:]]) )
        if j==9 and key[1]==4 and key[2]==0 and key[3] in [8,12]:
            for jpair in range(0, 10, 2):
                op = pair_creation(jpair, 9)
                if not save_memory:
                    print "off-diagonal TBME of J=%i with v=2    : %s" % \
                        ( jpair, 
                          bitwf_dot( bitwf_times_a(op, bitwf_vj[key[:1]+(2,)+key[2:]]),
                                     bitwf_times_a(op, bitwf_vj[key]) ) )
                print "off-diagonal TBME of J=%i with alpha=2: %s" % \
                    ( jpair, 
                      bitwf_dot( bitwf_times_a(op, bitwf_vj[key[:2]+(1,)+key[3:]]),
                                 bitwf_times_a(op, bitwf_vj[key]) ) )


if __name__ == "__main__":

    outputfile = open( "seniority_basis.dat", 'wb' )
    for j in range(1, maxj+2, 2):
        bitwf_vj, dim_vj = calc_seniority_basis(j)
        outputfile.write( struct.pack(fb+"i", j) )
        # print "j=%i/2"%j
        # for v,jj in sorted(dim_vj):
        #   print "v=%i j=%4s: %i"%(v,"%2i"%(jj/2) if jj%2==0 else "%2i/2"%jj,dim_vj[v,jj])
        # print
        # check_seniority_basis(bitwf_vj,dim_vj,j)
        if save_memory:
            sup = seniority_up(j)
            for n in range(j+2):
                for v in range(min(n,j+1-n),-1,-2)[::-1]:
                    for jj in sorted([x[1] for x in dim_vj if x[0]==v]):
                        for alpha in range(dim_vj[v,jj]):
                            for m in range(jj, -jj-2, -2):
                                if m == jj:
                                    bitwf = bitwf_vj[v, v, alpha, jj, jj]
                                    for i in range(v+2, n+2, 2):
                                        bitwf = bitwf_times_s( sf1/Sqrt_fraction((i-v)*(j-i-v+3)/4),
                                                               bitwf_times_c(sup, bitwf) )
                                else:
                                    bitwf = bitwf_times_s( sf1/Sqrt_fraction((jj+m+2)*(jj-m)/4),
                                                           bitwf_jdown(bitwf, j) )
                                print>>sys.stderr, "n=%2i v=%1i alpha=%2i/%2i j=%4s m=%5s:" % \
                                    ( n, v, alpha+1, dim_vj[v,jj],
                                      "%2i"%(jj/2) if jj%2==0 else "%2i/2"%jj,
                                      "%3i"%( m/2) if m %2==0 else "%3i/2"%m )
                                outputfile.write(struct.pack(
                                    fb+"6i", n, v, alpha, jj, m, len(bitwf) ) )
                                outputfile.write( struct.pack(
                                    fb+"id"*len(bitwf),
                                    *[ [bit,float(bitwf[bit])][x]
                                       for bit in sorted(bitwf)
                                       for x in [0,1] ] ) )
    #       print "n=%2i v=%1i alpha=%2i/%2i j=%4s m=%5s:"%(n,v,alpha+1,dim_vj[v,jj],
    #                                                       "%2i"%(jj/2) if jj%2==0 else "%2i/2"%jj,
    #                                                       "%3i"%(m/2) if m%2==0 else "%3i/2"%m)
    #       print str_bitwf(bitwf,j)
        else:
            for key in sorted( bitwf_vj, key=lambda x:(x[0],x[1],x[3],x[2],-x[4]) ):
                outputfile.write( struct.pack(
                    fb+"6i", key[0], key[1], key[2], key[3], key[4], len(bitwf_vj[key]) ))
                outputfile.write( struct.pack(
                    fb+"id"*len(bitwf_vj[key]),
                    *[ [bit, float(bitwf_vj[key][bit])][x]
                       for bit in sorted(bitwf_vj[key])
                       for x in [0,1] ] ))
    #   print "n=%2i v=%1i alpha=%2i/%2i j=%4s m=%5s:"%(key[0],key[1],key[2]+1,dim_vj[key[1],key[3]],
    #                                                   "%2i"%(key[3]/2) if key[3]%2==0 else "%2i/2"%key[3],
    #                                                   "%3i"%(key[4]/2) if key[4]%2==0 else "%3i/2"%key[4])
    #   print str_bitwf(bitwf_vj[key],j)
    outputfile.close()
