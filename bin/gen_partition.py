#!/usr/bin/env python
# ./gen_partition.py hoge.snt hoge.ptn
# generate partition file from snt file
# usage: gen_partiton.py hoge.snt #proton #neutron parity
#

import sys, operator, random, os.path
from functools import reduce

output_ans = ""

def raw_input_save(comment: str = "") -> str:
    """
    Fetch user input from 'input' and append it to 'output_ans'.

    Parameters
    ----------
    comment : str
        Comment displayed on screen when user is prompted for input.

    Returns
    -------
    ans : str
        Input from user.
    """
    # if c is None: r = input()
    # else: r = input(c)
    ans = input(comment)
    global output_ans
    output_ans += ans + '\n'
    return ans

def read_comment_skip(fp):
    """
    NOTE: This can probably be replaced by 'read_comment_skip' in 
    'kshell_ui'.
    """
    while True:
        arr = fp.readline().split()
        if not arr: return None
        for i in range(len(arr)): 
            if (arr[i][0] == "!") or (arr[i][0] == "#"): 
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

def orb2char(n, l, j, tz):
    lorb2c = [
        's', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o'
    ]
    tz2c = { -1:'p', 1:'n' }
    return "%c_%d%c%d/2" % (tz2c[tz], n, lorb2c[l], j)

class ModelSpace:
    def __init__(self,
        valence_p_n: tuple,
        norb: list,
        lorb: list,
        jorb: list,
        itorb: list,
    ):
        """
        Parameters
        ----------
        valence_p_n : tuple
            Number of valence protons and neutrons. For example, if the
            interaction is USDA and the nucleus is 22Ne, then this
            will be (2, 2).

        norb : list
            "Principal quantum number" of each orbital. Not really the
            principal quantum number as in atomic physics, but rather
            a counter which counts the number of times an orbital of a
            certain total angular momentum has occurred. For USDA we
            have the following model space:

                1     0   2   3  -1  !   1 = p 0d_ 3/2
                2     0   2   5  -1  !   2 = p 0d_ 5/2
                3     1   0   1  -1  !   3 = p 1s_ 1/2
                4     0   2   3   1  !   4 = n 0d_ 3/2
                5     0   2   5   1  !   5 = n 0d_ 5/2
                6     1   0   1   1  !   6 = n 1s_ 1/2

            where `norb` is [0, 0, 1, 0, 0, 1]. 0d_ 3/2 is the first
            d_ 3/2 orbital, while 1s_ 1/2 is the second s_ 1/2 orbital.
            This is because USDA has an inert 18O core which contains an
            s_ 1/2 orbital.

        lorb : list
            Orbital angular momentum of each orbital. For USDA we have
            `lorb` equal to [2, 2, 0, 2, 2, 0].

        jorb : list
            Total angular momentum of each orbital. For USDA we have
            `jorb` equal to [3, 5, 1, 3, 5, 1]. Note that the values are
            multiplied by 2 to avoid having to deal with floats.

        itorb : list
            Isospin of each orbital. For USDA we have `itorb` equal to
            [-1, -1, -1, 1, 1, 1]. -1 corresponds to protons, while 1
            corresponds to neutrons.
        """
        self.valence_p_n: tuple[int, int] = valence_p_n
        self.norb: list[int] = norb
        self.lorb: list[int] = lorb
        self.jorb: list[int] = jorb
        self.itorb: list[int] = itorb
        self.parityorb: list[int] = [(-1)**l for l in lorb]

        self.norb_pn = [
            [n for n, t in zip(norb, itorb) if t == -1],
            [n for n, t in zip(norb, itorb) if t == 1],
        ]
        self.lorb_pn = [
            [l for l, t in zip(lorb, itorb) if t == -1],
            [l for l, t in zip(lorb, itorb) if t == 1],
        ]
        self.jorb_pn = [    # Make separate jorb lists for protons and neutrons.
            [j for j, t in zip(jorb, itorb) if t == -1],
            [j for j, t in zip(jorb, itorb) if t == 1],
        ]
        self.iporb_pn = [
            [p for p, t in zip(self.parityorb, itorb) if t == -1],
            [p for p, t in zip(self.parityorb, itorb) if t == 1],
        ]
        self.phtrunc_t = []
        self.phtrunc_orb = []
        self.phtrunc_mask_pn = []
        self.phtrunc_maxt_pn = []
        self.phtrunc_mint_pn = []

        self.minhw, self.maxhw = 0, 0
        self.hworb_pn = [
            [0 for t in self.itorb if t == -1 ],
            [0 for t in self.itorb if t == 1],
        ]
        self.maxhw_pn, self.minhw_pn = (0, 0), (0, 0)
        self.is_called_hw = False
        self.is_monopole_trunc = False
        self.min_eocc = 1.e10

    def set_hw_for_phtrunc(self, phorb, mmhw):
        # phorb ... orbit list
        # mmhw ... min. and max occupation
        # can be called  once
        hworb = [0]*len(self.norb)
        for i in phorb: hworb[i] += 1
        self.hworb_pn = [ [ p for p,t in zip(hworb, self.itorb) 
                            if t == -1], 
                          [ p for p,t in zip(hworb, self.itorb) 
                            if t ==  1] ]
        self.set_hw_truncation(mmhw, is_hw_exct=False)

    def set_hw_truncation(self,
            hw_trunctation: list[int],
            is_hw_exct: bool = True,
        ):
        """
        Parameters
        ----------
        hw_trunctation : list[int]
            Minimum and maximum number of nucleons to be excited across
            the major shell gaps.
        """
        if self.is_called_hw:
            msg = "'set_hw_truncation' can only be called once but it was called twice."
            raise RuntimeError(msg)

        self.is_called_hw = True            
        self.minhw, self.maxhw = hw_trunctation
        
        if is_hw_exct:
            self.hworb_pn = [
                [2*n + l for n, l, t in zip(self.norb, self.lorb, self.itorb) if t == tz ] for tz in (-1, 1)
            ]
        
        lowest_pn, highest_pn = self.cal_hw_low_high_pn(self.valence_p_n)
        
        if is_hw_exct:
            self.minhw += sum(lowest_pn)
            self.maxhw += sum(lowest_pn)
            print("lowest hw, maxhw ", self.minhw, self.maxhw)
        
        self.maxhw_pn = (
            min(self.maxhw - lowest_pn[1], highest_pn[0]), min(self.maxhw - lowest_pn[0], highest_pn[1])
        )
        self.minhw_pn = (
            max(self.minhw - highest_pn[1], lowest_pn[0]), max(self.minhw - highest_pn[0], lowest_pn[1])
        )

    def set_ph_truncation(self, orb_list, t_list):
        self.phtrunc_t = t_list
        self.phtrunc_orb = orb_list

        for orb, t in zip(self.phtrunc_orb, self.phtrunc_t):
            mask = [0,]*len(self.norb)
            for i in orb: mask[i] = 1
            self.phtrunc_mask_pn.append(
                [ [mask[i] for i,t in enumerate(self.itorb) if t==-1], 
                  [mask[i] for i,t in enumerate(self.itorb) if t== 1] ] )

        lowest_pn, highest_pn = self.cal_phtrunc_t_low_high_pn(self.valence_p_n)

        self.phtrunc_maxt_pn = [
            ( min( pht[1]-lpn[1], hpn[0] ) , min( pht[1]-lpn[0], hpn[1] ) ) 
            for i, (pht, lpn, hpn) in enumerate( zip(
                    self.phtrunc_t, lowest_pn, highest_pn ) ) ]

        self.phtrunc_mint_pn = [
            ( max( pht[0]-hpn[1], lpn[0] ) , max( pht[0]-hpn[0], lpn[1] ) ) 
            for i, (pht, lpn, hpn) in enumerate( zip(
                    self.phtrunc_t, lowest_pn, highest_pn ) ) ]

        # for i in range(len(self.phtrunc_t)):
        #     self.phtrunc_maxt_pn.append(
        #         ( min( self.phtrunc_t[i][1] - lowest_pn[i][1], 
        #                highest_pn[i][0] ) , 
        #           min( self.phtrunc_t[i][1] - lowest_pn[i][0], 
        #                highest_pn[i][1] ) ) )
        #     self.phtrunc_mint_pn.append(
        #         ( max( self.phtrunc_t[i][0] - highest_pn[i][1], 
        #                lowest_pn[i][0] ) , 
        #           max( self.phtrunc_t[i][0] - highest_pn[i][0], 
        #                lowest_pn[i][1] ) ) )

    def set_monopole_truncation(self, model_space_filename, thd_energy):
        self.is_monopole_trunc = True
        from espe import SMInt
        self.SMInt = SMInt(model_space_filename)
        self.monopole_e_thd = thd_energy

    def gen_ptn_pn(self):
        """
        Variables
        ---------
        orb_hw : list[list[int]]
            `orb_hw` is a list of lists. Each sublist contains the
            occupation number for each orbital in the same major shell.

        self.hworb_pn : list[list[int]]
            `self.hworb_pn` is possibly a list of which major shell each
            orbital belongs to. For sdpf-mu the list is
            
                self.hworb_pn = [
                    [2, 2, 2, 3, 3, 3, 3], [2, 2, 2, 3, 3, 3, 3]
                ]
            
            which might mean that the 1d5/2, 1d3/2 and 2s1/2 orbitals
            belong to the same major shell (number 2), while 1f7/2,
            1f5/2, 2p3/2 and 2p1/2 belong to the same major shell
            (number 3). There is one sublist for protons and one for
            neutrons.

        self.jorb_pn : list[list[int]]
            Each sublist contains the total angular momentum for each
            orbital. One sublist for protons and one for neutrons. For
            sdpf-mu the list is

                self.jorb_pn = [
                    [5, 3, 1, 7, 5, 3, 1], [5, 3, 1, 7, 5, 3, 1]
                ]

        orb_nhw : list[int]
            `orb_nhw` is a list with the maximum occupation number for
            each major shell. For sdpf-mu the list is

                orb_nhw = [12, 20] 

            meaning that the sd shell can have a maximum of 12
            protons and 12 neutrons, while the pf shell can have a
            maximum of 20 protons and 20 neutrons.
        """
        # nocc_orb_pn = [ [ j+1 for j,t in zip(self.jorb, self.itorb) 
        #                   if t==-1], 
        #                 [ j+1 for j,t in zip(self.jorb, self.itorb) 
        #                   if t== 1] ]
        self.ptn_pn = [[], []]

        def gen_hw_nocc(orb_hw: list[list[int]], hwnocc):
            """
            Parameters
            ----------
            orb_hw : list[list[int]]
                Max. occupation number for each orbital. There is one
                sublist for each major shell. Hence, USDA will only have
                one sublist while sdpf-mu has two.
            """
            if self.valence_p_n == 0:
                yield (0,)*sum([len(i) for i in orb_hw])
                return
            if len(orb_hw) == 1:
                for i in self.gen_nocc(orb_hw[0], hwnocc[0]):
                    yield i
                return
            for i in self.gen_nocc(orb_hw[0], hwnocc[0]):
                for j in gen_hw_nocc(orb_hw[1:], hwnocc[1:]):
                    yield i + j

        print(f"{self.hworb_pn = }")
        print(f"{self.jorb_pn = }")
        for tz in range(2):
            hw_list, orb_hw = [], []
            hw0 = -0.1 # initialized, not integer
            ihw = 0
            for hw, j in zip(self.hworb_pn[tz], self.jorb_pn[tz]):
                if hw != hw0: 
                    hw_list.append(hw)
                    ihw += 1
                    hw0 = hw
                    orb_hw.append([j+1])
                else:
                    orb_hw[-1].append(j+1)
            orb_nhw = [ sum(arr) for arr in orb_hw ]
            hw_nocc = []
            print(f"{orb_hw = }")
            print(f"{orb_nhw = }")
            print(f"{hw_list = }")
            
            for arr in self.gen_nocc(orb_nhw, self.valence_p_n[tz]):
                nhw = sum( hw*n for hw,n in zip(hw_list, arr))
                if nhw > self.maxhw_pn[tz]: continue
                if nhw < self.minhw_pn[tz]: continue
                hw_nocc.append( arr )

            def check_trunc_pn( arr ):
                for i in range(len(self.phtrunc_t)):
                    if not (
                        self.phtrunc_mint_pn[i][tz] <= sum( [ m*n for m, n in zip(self.phtrunc_mask_pn[i][tz], arr)] ) <= self.phtrunc_maxt_pn[i][tz]
                    ):
                        return False
                return True

            for hwnocc in hw_nocc:
                for arr in gen_hw_nocc(orb_hw, hwnocc):
                    if check_trunc_pn( arr ):
                        self.ptn_pn[tz].append( arr )
            self.ptn_pn[tz].sort()

    def ptn_combined(self, parity):
        # parity
        self.ptn_pn_parity = [
            [ reduce(operator.mul, 
                     [ p**n for p,n in zip(self.iporb_pn[tz], arr) ])
              for arr in self.ptn_pn[tz] ]
            for tz in range(2) ]
        # hw 
        self.ptn_pn_hw = [ 
            [ sum( self.hworb_pn[tz][i]*arr[i] for i in range(len(arr)))
              for arr in self.ptn_pn[tz] ]
            for tz in range(2) ]
        # p-h truncation
        ptn_pn_phtrunc_t = []
        for orb, t in zip(self.phtrunc_orb, self.phtrunc_t):
            mask = [0,]*len(self.norb)
            for i in orb: mask[i] += 1
            maskpn = [ [mask[i] for i,t in enumerate(self.itorb) if t==-1], 
                       [mask[i] for i,t in enumerate(self.itorb) if t== 1] ]
            ptn_pn_phtrunc_t.append( 
                [ [ sum( n*m for n,m in zip(arr, maskpn[tz]) ) 
                    for arr in self.ptn_pn[tz] ] 
                  for tz in range(2) ] )
            # for i in orb: mask[i] = 1
            # maskpn = [ [mask[i] for i,t in enumerate(self.itorb) if t==-1], 
            #            [mask[i] for i,t in enumerate(self.itorb) if t== 1] ]
            # occorb = [ [i for i,m in enumerate(maskpn[tz]) if m==1]
            #            for tz in range(2)]
            # ptn_pn_phtrunc_t.append( 
            #     [ [ sum( arr[i] for i in occorb[tz] ) 
            #         for arr in self.ptn_pn[tz] ] 
            #       for tz in range(2) ] )
            
        def check_trunc(i_p, i_n):
            # parity
            if self.ptn_pn_parity[0][i_p] * self.ptn_pn_parity[1][i_n] \
               != parity: return False
            # hw excitation
            hw = self.ptn_pn_hw[0][i_p] + self.ptn_pn_hw[1][i_n]
            if not self.minhw <= hw <= self.maxhw: return False
            # ph trunc 
            for tpn, t in zip(ptn_pn_phtrunc_t, self.phtrunc_t):
                n = tpn[0][i_p] + tpn[1][i_n]
                if not t[0] <= n <= t[1]: 
                    return False
            return True

        # monopole trunc
        if self.is_monopole_trunc:
            ptn_list = [ (i, j)
                         for i in range(len(self.ptn_pn[0]))
                         for j in range(len(self.ptn_pn[1]))
                         if check_trunc(i, j) ]
            elist = [ self.SMInt.energy_occ(
                self.ptn_pn[0][i] + self.ptn_pn[1][j] )
                      for (i,j) in ptn_list ]
            self.min_eocc = min(elist)
            self.monopole_e_thd += self.min_eocc
            self.ptn_list = []
            for (i,j), e in zip(ptn_list, elist):
                nocc = self.ptn_pn[0][i] + self.ptn_pn[1][j]
                if e > self.monopole_e_thd:
                    print('SKIP partition', nocc, ' : %10.5f' % e)
                else:
                    self.ptn_list.append( (i,j) )
                    print('PASS partition', nocc, ' : %10.5f' % e)
            return
                
        
            
        self.ptn_list = [ (i, j)
                          for i in range(len(self.ptn_pn[0]))
                          for j in range(len(self.ptn_pn[1]))
                          if check_trunc(i, j) ]

    def strip_ptn_pn(self):
        is_ptn_pn = [ [False,]*len(self.ptn_pn[0]),  
                      [False,]*len(self.ptn_pn[1]) ]
        for pp,nn in self.ptn_list:
            is_ptn_pn[0][pp] = True
            is_ptn_pn[1][nn] = True
        map_orb = [ [-1,]*len(is_ptn_pn[0]), [-1,]*len(is_ptn_pn[1]) ]
        for tz in range(2):
            j = 0
            for i,f in enumerate(is_ptn_pn[tz]):
                if f: 
                    map_orb[tz][i] = j
                    j += 1
            ptn_pn = []
            for i,f in enumerate(map_orb[tz]):
                if f != -1:
                    ptn_pn.append( self.ptn_pn[tz][i] )
            self.ptn_pn[tz] = ptn_pn
        ptn_list = []
        for i,j in self.ptn_list:
            ni, nj = map_orb[0][i], map_orb[1][j]
            if ni == -1 or nj == -1: continue
            ptn_list.append( (ni, nj) )
        self.ptn_list = ptn_list
                                  
    def write_ptn_pn(self, fp, parity, model_space_filename):
        # output partition of proton and neutron separately
        fp.write( "# partition file of %s  Z=%d  N=%d  parity=%+d\n" 
                  % (model_space_filename, self.valence_p_n[0], self.valence_p_n[1], parity ) )
        fp.write( " %d %d %d\n" % (self.valence_p_n[0], self.valence_p_n[1], parity) )
        fp.write( "# num. of  proton partition, neutron partition\n" )
        fp.write( " %d %d\n" % (len(self.ptn_pn[0]), len(self.ptn_pn[1]) ))
        for tz in range(2):
            if tz==0: fp.write( "# proton partition\n" )
            if tz==1: fp.write( "# neutron partition\n" )
            for i, arr in enumerate(self.ptn_pn[tz]):
                fp.write( " %5d   " % (i+1,) )
                for a in arr:
                    fp.write( " %2d" % (a,) )
                fp.write("\n")

    def write_ptn_combined(self, fp):
        fp.write( "# partition of proton and neutron\n" )
        out = ""
        nline = len(self.ptn_list)
        # random.shuffle(ptn_list)  # shuffle order of p-n partitions
        out=""
        for i,j in self.ptn_list:
            out += "%5d %5d\n" % (i+1, j+1)
        fp.write( "%d\n" % (nline,) )
        fp.write( out )
        if len(self.ptn_list)==0: 
            sys.stdout.write( "\n *** WARNING NO PARTITION *** \n" )

    def cal_hw_low_high_pn(self, valence_p_n: tuple[int, int]):
        """
        total hw excitation of the lowest and highest configuration

        Parameters
        ----------
        valence_p_n : tuple[int, int]
            The number of valence protons and neutrons.
        """
        nhw = [ [], [] ]
        for tz in range(2):
            for i in range(len(self.jorb_pn[tz])):
                nhw[tz] += [ self.hworb_pn[tz][i], ]*(self.jorb_pn[tz][i] + 1)
        for tz in range(2): nhw[tz].sort()
        return ( sum(nhw[0][:valence_p_n[0]]), sum(nhw[1][:valence_p_n[1]]) ), \
            ( sum(nhw[0][-valence_p_n[0]:]), sum(nhw[1][-valence_p_n[1]:]) )

    def cal_phtrunc_t_low_high_pn(self, valence_p_n):
        lowest_pn = []
        highest_pn = []
        for mask_pn in self.phtrunc_mask_pn:
            nhw = [ [], [] ]
            for tz in range(2):
                for i in range(len(self.jorb_pn[tz])):
                    nhw[tz] += [ mask_pn[tz][i], ]*(self.jorb_pn[tz][i]+1) 
            for tz in range(2): nhw[tz].sort()
            lowest_pn.append((sum(nhw[0][:valence_p_n[0]]),sum(nhw[1][:valence_p_n[1]])))
            highest_pn.append((sum(nhw[0][-valence_p_n[0]:]),sum(nhw[1][-valence_p_n[1]:])))
        return lowest_pn, highest_pn

    def gen_nocc(self, nlist, valence_p_n):
        if valence_p_n == 0: 
            yield (0,)*len(nlist)
            return
        if len(nlist)==1:
            yield (valence_p_n,)
            return
        ns, nrest = nlist[0], nlist[1:]
        # for i in range(max(0, valence_p_n-sum(nrest)), min(ns, valence_p_n)+1): 
        for i in range(min(ns, valence_p_n), max(0, valence_p_n - sum(nrest)) -1, -1): 
            for j in self.gen_nocc(nrest, valence_p_n - i):
                yield (i,) + j

def main(
    model_space_filename: str,
    partition_filename: str,
    valence_p_n: tuple,
    parity: int
    ):
    """
    Parameters
    ----------
    model_space_filename:
        The filename of the model space (.snt) file.

    partition_filename:
        The filename of the partition (.ptn) file.

    valence_p_n:
        Tuple containing the number of valence protons and neutrons.
        Example: (#p, #n).
    
    parity:
        The parity of the current partition.
    """
    
    try:
        fp = open(model_space_filename, 'r')
    except FileNotFoundError:
        print(f"File '{model_space_filename=}' not found")
        sys.exit()
    
    n_jorb, n_core = [0, 0], [0, 0]
    n_jorb[0], n_jorb[1], n_core[0], n_core[1] = read_comment_skip(fp)
    norb, lorb, jorb, itorb = [], [], [], []
    
    for i in range(sum(n_jorb)):
        arr = read_comment_skip(fp)
        if i+1!=arr[0]: raise "read error"
        norb.append(arr[1])
        lorb.append(arr[2])
        jorb.append(arr[3])
        itorb.append(arr[4])
        if (i < n_jorb[0] and arr[4] != -1) or (i >= n_jorb[0] and arr[4] != 1): 
            raise "ERROR to read snt: proton orbit should come first"
    spe = [0.]*sum(n_jorb)
    nline = read_comment_skip(fp)
    for i in range(nline[0]):
        arr = read_comment_skip(fp)
        if arr[0] != arr[1]: continue
        spe[int(arr[0])-1] = float(arr[2])
    # print spe

    fp.close()

    class_ms = ModelSpace(
        valence_p_n, norb, lorb, jorb, itorb
    )

    # parity check
    prty_list = [ set(ip) for ip in class_ms.iporb_pn ]
    for i in range(2): 
        if valence_p_n[i] % 2 == 0 and prty_list[i] == set([-1]): 
            prty_list[i] = set( [1] )
        if valence_p_n[i] == 0: prty_list[i] = set( [1] )

    if parity == 1:
        if not (1 in prty_list[0] and 1 in prty_list[1] ) \
                and not (-1 in prty_list[0] and -1 in prty_list[1] ):
            print("No states in  positive parity")
            return
            # sys.exit()
    elif parity == -1:
        if not (1 in prty_list[0] and -1 in prty_list[1] ) \
                and not (-1 in prty_list[0] and 1 in prty_list[1] ):
            print("No states in negative parity")
            return
            # sys.exit()
    else:
        raise "illegal input"

    fpout = open(partition_filename, 'w')

    print(" truncation scheme ?\n" \
        + "      0 : No truncation (default) \n" \
        + "      1 : particle-hole truncation for orbit(s) \n" \
        + "      2 : hw truncation \n" \
        + "      3 : Both (1) and (2) \n")

    ans = raw_input_save()
    ans = ans.rstrip()
    if not ans: ans = 0
    truncation_mode = int(ans)

    if not 0 <= truncation_mode <= 4: raise 'input out of range'

    if (truncation_mode == 2) or (truncation_mode == 3):
        """
        Prompt user for hw truncation.
        """
        hw_truncation = raw_input_save( " (min. and) max hw for excitation : " )
        hw_truncation = hw_truncation.replace(",", " ").split()
        hw_truncation = [int(i) for i in hw_truncation]
        if len(hw_truncation) == 1: hw_truncation = [0, hw_truncation[0]]
        class_ms.set_hw_truncation(hw_truncation)
        fpout.write("# hw trucnation,  min hw = "+str(hw_truncation[0]) 
                    +" ,   max hw = "+str(hw_truncation[1])+"\n")

    if (truncation_mode == 1) or (truncation_mode == 3):
        print("   #    n   l   j   tz    spe ")
        for i in range(len(norb)):
            """
            Print the properties of the available valence orbitals.

            Example:
            #    n   l   j   tz    spe
            1    0   2   3  -1     1.980     p_0d3/2
            2    0   2   5  -1    -3.944     p_0d5/2
            3    1   0   1  -1    -3.061     p_1s1/2
            4    0   2   3   1     1.980     n_0d3/2
            5    0   2   5   1    -3.944     n_0d5/2
            6    1   0   1   1    -3.061     n_1s1/2
            """
            n, l, j, tz = norb[i], lorb[i], jorb[i], itorb[i],
            # print(" %3d  %3d %3d %3d %3d %9.3f     %s" \
            #     % ( i+1, n, l, j, tz, spe[i], orb2char(n, l, j, tz) ))
            msg = f" {i + 1:3d}  {n:3d} {l:3d} {j:3d} {tz:3d} {spe[i]:9.3f}"
            msg += f"     {orb2char(n, l, j, tz)}"
            print(msg)
        
        print(' specify # of orbit(s) and min., max. occupation numbers ' \
            + 'for restriction')
        orb_list, t_list = [], []
        while True:
            ans = raw_input_save(
                "\n # of orbit(s) for restriction?  (<CR> to quit): ")
            ans = ans.replace(',', ' ').split()
            if not ans: break
            orb_list.append( [int(i)-1 for i in ans] )
            ans = raw_input_save(
                ' min., max. restricted occupation numbers' 
                + 'for the orbit(s) (or max only) : ')
            ans = ans.replace(',', ' ').split()
            if len(ans) == 1: ans = [0,] + ans
            if len(ans) != 2: raise 'read error'
            t_list.append( [int(i) for i in ans] )

        if truncation_mode==1 and len(orb_list)>0:
            # class_ms.set_ph_truncation(orb_list[:], t_list[:])
            class_ms.set_hw_for_phtrunc(orb_list[0], t_list[0])
            if len(orb_list)>1:
               class_ms.set_ph_truncation(orb_list[1:], t_list[1:])
        else: # truncation_mode == 3
            class_ms.set_ph_truncation(orb_list, t_list)
        fpout.write("# particle-hole truncation orbit(s) : min., max.\n")
        for orb,t in zip(orb_list, t_list):
            fpout.write("#      " + str([i+1 for i in orb]) + " :  " + \
                            str(t[0]) + " " + str(t[1]) + "\n")
    if truncation_mode == 4:
        ans = raw_input_save(
            " monopole trunc, threashold energy (relative to min): " )
        thd = float(ans)
        fpout.write( "# monopole-based partition truncation, thd= %10.5f\n"
                     %  thd)
        class_ms.set_monopole_truncation(model_space_filename, thd)

    sys.stdout.write( "generating partition file ..." )
    sys.stdout.flush()

    class_ms.gen_ptn_pn()

    sys.stdout.write( "..." )
    if class_ms.is_monopole_trunc:     sys.stdout.write( "\n" )
    sys.stdout.flush()

    class_ms.ptn_combined(parity)

    sys.stdout.write( "..." )
    sys.stdout.flush()

    class_ms.strip_ptn_pn()

    sys.stdout.write( "..." )
    sys.stdout.flush()

    class_ms.write_ptn_pn(fpout, parity, model_space_filename)
    class_ms.write_ptn_combined(fpout)

    sys.stdout.write( " done.\n" )
    if class_ms.is_monopole_trunc:
        print('\nminimum energy for partition %10.5f, threashold %10.5f\n'\
            % (class_ms.min_eocc, class_ms.monopole_e_thd))

    fpout.close()

    ret = None
    if 'orb_list' in locals():
        # orbit list in the first truncation
        ret = [i+1 for i in orb_list[0]] 
    return ret

if __name__ == "__main__":

    if len(sys.argv)<3: 
        print('usage: gen_partiton.py hoge.snt ' \
            + 'output.ptn #proton #neutron parity')
        sys.exit(1)

    if os.path.exists(sys.argv[2]): raise "partition file exists"

    model_space_filename, fn_out = sys.argv[1], sys.argv[2]
    valence_p_n = (int(sys.argv[3]), int(sys.argv[4]))
    parity = 1
    if len(sys.argv) > 5: 
        if   sys.argv[5] == "+": parity =  1
        elif sys.argv[5] == "-": parity = -1
        else: parity = int(sys.argv[5])

    if not parity in (1, -1): raise "parity error"

    main(model_space_filename, fn_out, valence_p_n, parity)
