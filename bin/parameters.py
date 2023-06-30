GS_FREE_PROTON = 5.585
GS_FREE_NEUTRON = -3.826

recommended_quenching_factors = {
    "GCLSTsdpfsdgix5pn.snt": f"0.75*GS_FREE = {round(0.75*GS_FREE_PROTON, 3), round(0.75*GS_FREE_NEUTRON, 3)}",
    "gs8.snt": f"0.75*GS_FREE = {round(0.75*GS_FREE_PROTON, 3), round(0.75*GS_FREE_NEUTRON, 3)}",
    "jun45.snt": f"0.7*GS_FREE = {round(0.7*GS_FREE_PROTON, 3), round(0.7*GS_FREE_NEUTRON, 3)}",
    "gxpf1a.snt": f"0.9*GS_FREE = {round(0.9*GS_FREE_PROTON, 3), round(0.9*GS_FREE_NEUTRON, 3)}",
    "gxpf1.snt": f"0.9*GS_FREE = {round(0.9*GS_FREE_PROTON, 3), round(0.9*GS_FREE_NEUTRON, 3)}",
    "sdpf-mu.snt": f"0.9*GS_FREE = {round(0.9*GS_FREE_PROTON, 3), round(0.9*GS_FREE_NEUTRON, 3)}"
}

element = [
    'NA', 
    'H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne', 
    'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca',
    'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
]