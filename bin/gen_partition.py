#!/usr/bin/env python
# ./gen_partition.py hoge.snt hoge.ptn
# generate partition file from snt file
# usage: gen_partiton.py hoge.snt #proton #neutron parity
#

import sys, operator, os.path
from functools import reduce
from io import TextIOWrapper

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
        n_valence_pn: tuple,
        norb: list,
        lorb: list,
        jorb: list,
        itorb: list,
    ):
        """
        Parameters
        ----------
        n_valence_pn : tuple
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
        self.n_valence_pn: tuple[int, int] = n_valence_pn
        self.norb: list[int] = norb
        self.lorb: list[int] = lorb
        self.jorb: list[int] = jorb
        self.itorb: list[int] = itorb
        self.parityorb: list[int] = [(-1)**l for l in lorb]

        self.norb_pn: list[list[int]] = [
            [n for n, t in zip(norb, itorb) if t == -1],
            [n for n, t in zip(norb, itorb) if t == 1],
        ]
        self.lorb_pn: list[list[int]] = [
            [l for l, t in zip(lorb, itorb) if t == -1],
            [l for l, t in zip(lorb, itorb) if t == 1],
        ]
        self.jorb_pn: list[list[int]] = [    # Make separate jorb lists for protons and neutrons.
            [j for j, t in zip(jorb, itorb) if t == -1],
            [j for j, t in zip(jorb, itorb) if t == 1],
        ]
        self.parityorb_pn: list[list[int]] = [
            [p for p, t in zip(self.parityorb, itorb) if t == -1],
            [p for p, t in zip(self.parityorb, itorb) if t == 1],
        ]
        self.phtrunc_t = []
        self.phtrunc_orb = []
        self.phtrunc_mask_pn = []
        self.phtrunc_maxt_pn = []
        self.phtrunc_mint_pn = []

        self.minhw, self.maxhw = 0, 0   # Used for limiting the possible proton-neutron configuration combinations based on the hw truncation.
        self.hworb_pn: list[list[int]] = [
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
            hw_truncation: list[int],
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
        self.minhw, self.maxhw = hw_truncation
        
        if is_hw_exct:
            self.hworb_pn = [
                [2*n + l for n, l, t in zip(self.norb, self.lorb, self.itorb) if t == tz ] for tz in (-1, 1)
            ]
        # if is_hw_exct:
            # self.hworb_pn = []
            # for tz in (-1, 1):
            #     tmp_list = []
            #     for n, l, t in zip(self.norb, self.lorb, self.itorb):
            #         if t == tz:
            #             tmp_list.append(2*n + l)
            #     self.hworb_pn.append(tmp_list)

        lowest_pn, highest_pn = self.cal_hw_low_high_pn(self.n_valence_pn)
        
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

        lowest_pn, highest_pn = self.cal_phtrunc_t_low_high_pn(self.n_valence_pn)

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

    def set_monopole_truncation(self, filename_interaction, thd_energy):
        self.is_monopole_trunc = True
        from espe import SMInt
        self.SMInt = SMInt(filename_interaction)
        self.monopole_e_thd = thd_energy

    def gen_ptn_pn(self):
        """
        Variables
        ---------
        orb_hw : list[list[int]]
            `orb_hw` is a list of lists. Each sublist contains the
            occupation number for each orbital in the same major shell.
            For `sdpf-mu` it is

                orb_hw = [[6, 4, 2], [8, 6, 4, 2]]

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

        def gen_hw_nocc(orb_hw: list[list[int]], hwnocc, n_valence: int | None = None):
            """
            Parameters
            ----------
            orb_hw : list[list[int]]
                Max. occupation number for each orbital. There is one
                sublist for each major shell. Hence, USDA will only have
                one sublist while sdpf-mu has two.
            """
            # if self.n_valence_pn == 0:
            if n_valence == 0:
                """
                Some combinations of interaction and nucleus might have
                0 valence protons or neutrons. Eg. sn100pn with any Sn
                isotope.
                """
                yield (0,)*sum([len(i) for i in orb_hw])
                return
            if len(orb_hw) == 1:
                """
                If the model space has only one major shell or if the
                recursive calls of this function has reached the final
                major shell.
                """
                for i in self.gen_nocc(orb_hw[0], hwnocc[0]):
                    yield i
                return
            for i in self.gen_nocc(orb_hw[0], hwnocc[0]):
                for j in gen_hw_nocc(orb_hw[1:], hwnocc[1:]):
                    yield i + j

        for tz in range(2):
            hw_list: list[int] = []
            orb_hw: list[list[int]] = []
            hw0 = -0.1 # initialized, not integer
            ihw = 0 # I dont think this is used for anything.
            for hw, j in zip(self.hworb_pn[tz], self.jorb_pn[tz]):
                """
                `hw` is an index which specifies the major shell of each
                orbital. Each orbital has an hw index. For USDA hw is 0
                for all orbitals while for sdpf-mu hw is 2 and 3 for the
                different orbitals. `j` is the total angular momentum of
                the accompanying orbital.
                """
                if hw != hw0:
                    """
                    The fist time this value of `hw` shows up. Append
                    the `hw` index to the `hw_list` and create a sublist
                    in `orb_hw` for the current major shell.
                    """
                    hw_list.append(hw)
                    ihw += 1
                    hw0 = hw
                    orb_hw.append([j+1])    # Append max occupation number for the current orbital.
                else:
                    """
                    A sublist for `hw` already exists, so append the
                    max occupation number for the current orbital to
                    that sublist.
                    """
                    orb_hw[-1].append(j+1)

            orb_nhw = [ sum(arr) for arr in orb_hw ]    # Sum the max occupation number for each major shell.
            hw_nocc = []
            for arr in self.gen_nocc(orb_nhw, self.n_valence_pn[tz]):
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
                """
                For 44Sc sdpf-mu with 2hw truncation:

                    hw_nocc = [(12, 1), (11, 2), (10, 3)]
                    hw_nocc = [(12, 3), (11, 4), (10, 5)]

                For 44Sc sdpf-mu with no truncation:

                    hw_nocc = [(13,)]
                    hw_nocc = [(15,)]

                for 20Ne USDA:
                
                    hw_nocc = [(2,)]
                    hw_nocc = [(2,)]
                """
                for arr in gen_hw_nocc(orb_hw, hwnocc, n_valence=self.n_valence_pn[tz]):
                    if check_trunc_pn( arr ):
                        self.ptn_pn[tz].append( arr )
            self.ptn_pn[tz].sort()

    def ptn_combined(self, parity):
        # parity
        # self.ptn_pn_parity = [
        #     [reduce(operator.mul, [p**n for p, n in zip(self.parityorb_pn[tz], arr)] + [1]) for arr in self.ptn_pn[tz]] for tz in range(2)
        # ]
        self.ptn_pn_parity: list[list[int]] = []    # List of two lists, one for protons and one for neutrons, each containing the parity of each configuration.

        for isospin in range(2):
            """
            Isospin is 0 for protons and 1 for neutrons.
            """
            configuration_parity: list[int] = []    # To store the parity of each configuration.

            for configuration in self.ptn_pn[isospin]:
                """
                Calculate the parity of each configuration.
                `configuration` is a list of occupation numbers for each 
                orbital.
                """
                occupation_parity: list[int] = []

                for p, occupation in zip(self.parityorb_pn[isospin], configuration):
                    """
                    Take the parity of each orbital to the power of the
                    number of particles in that orbital, giving the
                    parity of each orbital occupation.
                    """
                    occupation_parity.append(p**occupation)
                
                occupation_parity.append(1)  # In case the list is empty?
                configuration_parity.append(reduce(operator.mul, occupation_parity))    # This is where the parity of the configuration is calculated.
            
            self.ptn_pn_parity.append(configuration_parity)

        # hw 
        # self.ptn_pn_hw = [ [ sum( self.hworb_pn[tz][i]*arr[i] for i in range(len(arr))) for arr in self.ptn_pn[tz] ] for tz in range(2) ]

        self.ptn_pn_hw: list[list[int]] = [] # List of two lists, one for protons and one for neutrons

        for isospin in range(2):
            tmp_list = []
            
            for configuration in self.ptn_pn[isospin]:
                sum_elements = 0
                
                # for i in range(len(configuration)):
                #     sum_elements += self.hworb_pn[isospin][i]*configuration[i]

                for hw_idx, occupation in zip(self.hworb_pn[isospin], configuration):
                    sum_elements += hw_idx*occupation
                
                tmp_list.append(sum_elements)
            
            self.ptn_pn_hw.append(tmp_list)

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
            
        def check_trunc(p_idx: int, n_idx: int) -> bool:
            """
            Check if the combination of proton configuration `p_idx` and
            neutron configuration `n_idx` satisfies the truncation
            conditions.

            Parameters
            ----------
            p_idx : int
                Index of the proton configuration.

            n_idx : int
                Index of the neutron configuration.

            Returns
            -------
            bool
                True if the truncation conditions are satisfied, False
                otherwise.
            """
            if self.ptn_pn_parity[0][p_idx]*self.ptn_pn_parity[1][n_idx] != parity:
                """
                Parity of the proton configuration times the parity of
                the neutron configuration must equal the parity of the
                partition file.
                """
                return False

            # hw excitation
            hw = self.ptn_pn_hw[0][p_idx] + self.ptn_pn_hw[1][n_idx]
            if not self.minhw <= hw <= self.maxhw:
                """
                
                """
                return False
            
            # ph trunc
            for tpn, t in zip(ptn_pn_phtrunc_t, self.phtrunc_t):
                n = tpn[0][p_idx] + tpn[1][n_idx]
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
                
        # self.ptn_list = [ (i, j) for i in range(len(self.ptn_pn[0])) for j in range(len(self.ptn_pn[1])) if check_trunc(i, j) ]

        self.ptn_list: list[tuple[int, int]] = []
        for p_idx in range(len(self.ptn_pn[0])):
            """
            Loop over all proton and neutron indices and check if the
            truncation conditions are satisfied. If so, add the indices
            to the list of allowed proton-neutron configurations.
            """
            for n_idx in range(len(self.ptn_pn[1])):
                
                if check_trunc(p_idx, n_idx):
                    self.ptn_list.append((p_idx, n_idx))

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
                                  
    def write_ptn_pn(self,
        fp: TextIOWrapper,
        parity: int,
        filename_interaction: str,
    ):
        """
        Write the proton and neutron partitions to the partition file.

        Parameters
        ----------
        fp : TextIOWrapper
            The file pointer of the partition file.

        parity : int
            The parity of the partition.

        filename_interaction : str
            The filename of the interaction file.
        """
        # fp.write( "# partition file of %s  Z=%d  N=%d  parity=%+d\n" % (filename_interaction, self.n_valence_pn[0], self.n_valence_pn[1], parity ) )
        # fp.write( " %d %d %d\n" % (self.n_valence_pn[0], self.n_valence_pn[1], parity) )
        # fp.write( " %d %d\n" % (len(self.ptn_pn[0]), len(self.ptn_pn[1]) ))
        fp.write(f"# partition file of {filename_interaction}  Z={self.n_valence_pn[0]}  N={self.n_valence_pn[1]}  parity={parity}\n")
        fp.write(f" {self.n_valence_pn[0]} {self.n_valence_pn[1]} {parity}\n")
        fp.write( "# num. of  proton partition, neutron partition\n" )
        fp.write(f" {len(self.ptn_pn[0])} {len(self.ptn_pn[1])}\n")
        
        for tz in range(2):
            """
            Isospin (tz) 0 is proton and 1 is neutron
            """
            if tz == 0: fp.write( "# proton partition\n" )
            if tz == 1: fp.write( "# neutron partition\n" )
            
            for i, arr in enumerate(self.ptn_pn[tz]):
                """
                `i` is the index of the partition. `arr` is the
                occupation numbers for a specific configuration.
                """
                # fp.write( " %5d   " % (i+1,) )
                fp.write(f" {i+1:5d}   ")
                
                for a in arr:
                    # fp.write( " %2d" % (a,) )
                    fp.write(f" {a:2d}")
                
                fp.write("\n")

    def write_ptn_combined(self, fp: TextIOWrapper):
        """
        Write the combinations of proton and neutron partitions to the
        partition file.

        Parameters
        ----------
        fp : TextIOWrapper
            The file pointer of the partition file.
        """
        fp.write( "# partition of proton and neutron\n" )
        out = ""
        nline = len(self.ptn_list)
        
        if nline == 0:
            sys.stdout.write( "\n *** WARNING NO PARTITION *** \n" )
        
        for i, j in self.ptn_list:
            # out += "%5d %5d\n" % (i+1, j+1)
            out += f"{i+1:5d} {j+1:5d}\n"
        
        fp.write( "%d\n" % (nline,) )
        fp.write( out )

    def cal_hw_low_high_pn(self, n_valence_pn: tuple[int, int]):
        """
        total hw excitation of the lowest and highest configuration

        Parameters
        ----------
        n_valence_pn : tuple[int, int]
            The number of valence protons and neutrons.
        """
        nhw = [ [], [] ]
        
        for isospin in range(2):
            for i in range(len(self.jorb_pn[isospin])):
                """
                NOTE: `self.jorb_pn` contains the total angular momentum
                of the orbitals. `self.jorb_pn[0]` is the j list for
                protons and `self.jorb_pn[1]` is the j list for
                neutrons.

                `self.hworb_pn` contains the harmonic oscillator quanta
                of the orbitals.
                """
                # nhw[isospin] += [ self.hworb_pn[isospin][i], ]*(self.jorb_pn[isospin][i] + 1)
                nhw[isospin].extend([self.hworb_pn[isospin][i]]*(self.jorb_pn[isospin][i] + 1))
        
        for isospin in range(2): nhw[isospin].sort()
        print(f"{nhw = }")
        print(f"{sum(nhw[0][:n_valence_pn[0]]) = }")
        return ( sum(nhw[0][:n_valence_pn[0]]), sum(nhw[1][:n_valence_pn[1]]) ), ( sum(nhw[0][-n_valence_pn[0]:]), sum(nhw[1][-n_valence_pn[1]:]) )

    def cal_phtrunc_t_low_high_pn(self, n_valence_pn):
        lowest_pn = []
        highest_pn = []
        for mask_pn in self.phtrunc_mask_pn:
            nhw = [ [], [] ]
            for tz in range(2):
                for i in range(len(self.jorb_pn[tz])):
                    nhw[tz] += [ mask_pn[tz][i], ]*(self.jorb_pn[tz][i]+1) 
            for tz in range(2): nhw[tz].sort()
            lowest_pn.append((sum(nhw[0][:n_valence_pn[0]]),sum(nhw[1][:n_valence_pn[1]])))
            highest_pn.append((sum(nhw[0][-n_valence_pn[0]:]),sum(nhw[1][-n_valence_pn[1]:])))
        return lowest_pn, highest_pn

    def gen_nocc(self, nlist, n_valence_pn: int):
        """
        Parameters
        ----------
        n_valence_pn : int
            The number of valence protons or neutrons.
        Returns
        -------
        For 44Sc sdpf-mu, returns a generator with the following:

            [(12, 1), (11, 2), (10, 3), (9, 4), (8, 5), (7, 6), (6, 7),
            (5, 8), (4, 9), (3, 10), (2, 11), (1, 12), (0, 13)]

            [(12, 3), (11, 4), (10, 5), (9, 6), (8, 7), (7, 8), (6, 9),
            (5, 10), (4, 11), (3, 12), (2, 13), (1, 14), (0, 15)]
        """
        if n_valence_pn == 0: 
            yield (0,)*len(nlist)
            return
        if len(nlist) == 1:
            yield (n_valence_pn,)
            return
        ns, nrest = nlist[0], nlist[1:]
        # for i in range(max(0, n_valence_pn-sum(nrest)), min(ns, n_valence_pn)+1): 
        for i in range(min(ns, n_valence_pn), max(0, n_valence_pn - sum(nrest)) -1, -1): 
            for j in self.gen_nocc(nrest, n_valence_pn - i):
                yield (i,) + j

def main(
    filename_interaction: str,
    filename_partition: str,
    n_valence_pn: tuple,
    parity: int
    ):
    """
    Parameters
    ----------
    filename_interaction:
        The filename of the model space (.snt) file.

    filename_partition:
        The filename of the partition (.ptn) file.

    n_valence_pn:
        Tuple containing the number of valence protons and neutrons.
        Example: (#p, #n).
    
    parity:
        The parity of the current partition.
    """
    
    try:
        fp = open(filename_interaction, 'r')
    except FileNotFoundError:
        print(f"File '{filename_interaction=}' not found")
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
        n_valence_pn, norb, lorb, jorb, itorb
    )
    # print(f"{n_valence_pn = }")
    # print(f"{norb = }")
    # print(f"{lorb = }")
    # print(f"{jorb = }")
    # print(f"{itorb = }")
    # return

    # parity check
    prty_list = [ set(ip) for ip in class_ms.parityorb_pn ]
    for i in range(2): 
        if n_valence_pn[i] % 2 == 0 and prty_list[i] == set([-1]): 
            prty_list[i] = set( [1] )
        if n_valence_pn[i] == 0: prty_list[i] = set( [1] )

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

    fpout = open(filename_partition, 'w')

    print(" truncation scheme ?\n" \
        + "      0 : No truncation (default) \n" \
        + "      1 : particle-hole truncation for orbit(s) \n" \
        + "      2 : hw truncation \n" \
        + "      3 : Both (1) and (2) \n" \
        + "      4 : Monopole-based partition truncation ")

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

        if (truncation_mode == 1) and len(orb_list) > 0:
            # class_ms.set_ph_truncation(orb_list[:], t_list[:])
            class_ms.set_hw_for_phtrunc(orb_list[0], t_list[0])
            if len(orb_list) > 1:
               class_ms.set_ph_truncation(orb_list[1:], t_list[1:])
        
        else: # truncation_mode == 3
            class_ms.set_ph_truncation(orb_list, t_list)
        
        fpout.write("# particle-hole truncation orbit(s) : min., max.\n")

        for orb, t in zip(orb_list, t_list):
            fpout.write("#      " + str([i+1 for i in orb]) + " :  " + \
                            str(t[0]) + " " + str(t[1]) + "\n")
    if truncation_mode == 4:
        ans = raw_input_save(
            " monopole trunc, threshold energy (relative to min): "
        )
        thd = float(ans)
        fpout.write("# monopole-based partition truncation, thd= %10.5f\n" %  thd)
        class_ms.set_monopole_truncation(filename_interaction, thd)

    sys.stdout.write( "generating partition file ..." )
    sys.stdout.flush()

    class_ms.gen_ptn_pn()

    sys.stdout.write( "..." )
    if class_ms.is_monopole_trunc: sys.stdout.write( "\n" )
    sys.stdout.flush()

    class_ms.ptn_combined(parity)

    sys.stdout.write( "..." )
    sys.stdout.flush()

    class_ms.strip_ptn_pn()

    sys.stdout.write( "..." )
    sys.stdout.flush()

    class_ms.write_ptn_pn(fpout, parity, filename_interaction)
    class_ms.write_ptn_combined(fpout)

    sys.stdout.write( " done.\n" )
    if class_ms.is_monopole_trunc:
        print('\nminimum energy for partition %10.5f, threshold %10.5f\n'\
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

    filename_interaction, fn_out = sys.argv[1], sys.argv[2]
    n_valence_pn = (int(sys.argv[3]), int(sys.argv[4]))
    parity = 1
    if len(sys.argv) > 5: 
        if   sys.argv[5] == "+": parity =  1
        elif sys.argv[5] == "-": parity = -1
        else: parity = int(sys.argv[5])

    if not parity in (1, -1): raise "parity error"

    main(filename_interaction, fn_out, n_valence_pn, parity)
