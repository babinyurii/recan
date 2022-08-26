# -*- coding: utf-8 -*-


from math import log, sqrt

"""
U	Uracil (RNA)	U
W	Weak	A/T
S	Strong	C/G
M	Amino	A/C
K	Keto	G/T
R	Purine	A/G
Y	Pyrimidine	C/T
B	Not A	C/G/T
D	Not C	A/G/T
H	Not G	A/C/T
V	Not T	A/C/G
N	Any	A/C/G/T
"""
DEGENERATE_NUCS = ("U", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N")

def estimate_nucleotide_frequencies(seq):
    
    seq = seq.replace("-","").upper()
    
    A = seq.count("A")
    C = seq.count("C")
    G = seq.count("G")
    T = seq.count("T")
    
    length = float(len(seq))
    
    return [ x / length for x in [A, C ,G, T] ]
    
    
def p_distance(seq1, seq2):
    """calculates pairwise distance between two sequences
    distance = num of different nucleotides / total nucleotides
    """
    
    different_nucs = 0
    nuc_pairs = []
    
    for x in zip(seq1, seq2):
        if '-' not in x and x[0] not in DEGENERATE_NUCS and x[1] not in DEGENERATE_NUCS: # skip gaps
            nuc_pairs.append(x)
            
    for (nuc_1, nuc_2) in nuc_pairs:
        if nuc_1 != nuc_2:
            different_nucs += 1
    
    total_nucs = len(nuc_pairs)
    
    try:
        distance = float(different_nucs / total_nucs)
        
        return distance
    
    except ValueError:
        print("""the reasons for the ValueError maybe: 1. too many gaps in some
              region of the alignment. 2. too short window span""")
    
        
    
    
def jc_distance(seq1, seq2):
    
    """ 
    Jukes-Cantor 
    jc_distance = - b log(1 - p_dist / b)
    ------------
    b is a constant. b = 3/4 for nucleotides and 19/20 for proteins.
    p_dist = pairwise distance, which is uncorrected distance between seq1 and seq2
    """
    b = 0.75
    p_dist = p_distance(seq1,seq2)
    
    try: 
        distance = - b * log(1 - p_dist / b)
        
        return distance
    
    except ValueError: 
        print("""the reasons for the ValueError maybe: 1. too many gaps in some
              region of the alignment. 2. too short window span""")
        

       
def k2p_distance(seq1,seq2):
    """
    Kimura 2-Parameter distance 
    k2p_distance = - 0.5 log( (1 - 2p - q)*sqrt(1 - 2q) )
    where:
    p : transition frequency
    q : transversion frequency
    """
    nuc_pairs = []

    # collect nuc pairs without gaps
    for x in zip(seq1,seq2):
        if '-' not in x and x[0] not in DEGENERATE_NUCS and x[1] not in DEGENERATE_NUCS: 
            nuc_pairs.append(x)
        
    ts_count = 0
    tv_count = 0
    total_nucs = len(nuc_pairs)
    
    transitions = [ "AG", "GA", "CT", "TC"]
    transversions = [ "AC", "CA", "AT", "TA",
                      "GC", "CG", "GT", "TG" ]

    for (nuc_1, nuc_2) in nuc_pairs:
        if nuc_1 + nuc_2 in transitions: 
            ts_count += 1 
        elif nuc_1 + nuc_2 in transversions: 
            tv_count += 1
    
    ts_freq = float(ts_count) / total_nucs
    tv_freq = float(tv_count) / total_nucs
    
    try: 
        distance = -0.5 * log((1 - 2 * ts_freq - tv_freq) * sqrt( 1 - 2 * tv_freq ))
        
        return distance
    
    except ValueError: 
        print("""the reasons for the ValueError maybe: 1. too many gaps in some
              region of the alignment. 2. too short window span""")
        
     
def tamura_distance(seq1,seq2):
    """
    Tamura distance = -C log( 1 - P/C - Q ) - 0.5( 1 - C )log( 1 - 2Q )
    where:
    P = transition frequency
    Q = transversion frequency
    C = GC1 + GC2 - 2 * GC1 * GC2
    GC1 = GC-content of sequence 1
    GC2 = GC-coontent of sequence 2
    """
    
    nuc_pairs = []
    
    #collect ungapped pairs
    for x in zip(seq1,seq2):
        if '-' not in x and x[0] not in DEGENERATE_NUCS and x[1] not in DEGENERATE_NUCS: 
            nuc_pairs.append(x)
        
    ts_count = 0
    tv_count = 0
    total_nucs = len(nuc_pairs)
    
    transitions = [ "AG", "GA", "CT", "TC"]
    transversions = [ "AC", "CA", "AT", "TA",
                      "GC", "CG", "GT", "TG" ]

    for (nuc_1, nuc_2) in nuc_pairs:
        if nuc_1 + nuc_2 in transitions: 
            ts_count += 1 
        elif nuc_1 + nuc_2 in transversions: 
            tv_count += 1
    
    ts_freq = float(ts_count) / total_nucs # p
    tv_freq = float(tv_count) / total_nucs # q
    
    gc1 = sum(estimate_nucleotide_frequencies(seq1)[1:3])
    gc2 = sum(estimate_nucleotide_frequencies(seq2)[1:3])
    
    c = gc1 + gc2 - 2 * gc1 * gc2

    try: 
        distance = -c * log( 1 - ts_freq / c - tv_freq) - 0.5 * ( 1 - c ) * log ( 1 - 2 * tv_freq )
        
        return distance
    
    except ValueError: 
        print("""the reasons for the ValueError maybe: 1. too many gaps in some
              region of the alignment. 2. too short window span""")
        
    
# TODO tajima nei isn't included still
def tn_distance(seq1, seq2):
    """ 
    Tajima-Nei distance = -b log(1 - p / b)
    where:
    b = 0.5 * [ 1 - Sum i from A to T(Gi^2+p^2/h) ]
    h = Sum i from A to G( Sum j from C to T (Xij^2/2*Gi*Gj))
    p = p-distance, i.e. uncorrected distance between seq1 and seq2
    Xij = frequency of pair (i,j) in seq1 and seq2, with gaps removed
    Gi = frequency of base i over seq1 and seq2 """
    from math import log

    ns = ['A','C','G','T']
    G = estimate_nucleotide_frequencies(seq1 + seq2)
    p = p_distance(seq1,seq2)
    pairs = []
    h = 0

    #collect ungapped pairs
    for x in zip(seq1,seq2):
        if '-' not in x: pairs.append(x)
       
    #pair frequencies are calculated for AC, AG, AT, CG, CT, GT (and reverse order)
    for i in range(len(ns)-1):
        for j in range(i+1,len(ns)):
            if i != j: 
                paircount = pairs.count( (ns[i], ns[j]) ) + pairs.count( (ns[j], ns[i]) )
            Xij_sq = (float(paircount)/len(pairs))**2
            GiGj = G[i]*G[j]
            h += 0.5*Xij_sq/GiGj  #h value used to calculate b
    
    b = 0.5*(1-sum([x**2 for x in G])+p**2/h)
    
    try: 
        d = -b * log(1 - p/b)
        return d 
    
    except ValueError: 
        print("""the reasons for the ValueError maybe: 1. too many gaps in some
              region of the alignment. 2. too short window span""")
        
       
def calc_pairwise_distance(seq1, seq2, dist_method):
    
    if dist_method == "pdist":
        distance = p_distance(seq1, seq2)
        
    elif dist_method == "jcd":
        distance = jc_distance(seq1, seq2)
        
    elif dist_method == "k2p":
        distance = k2p_distance(seq1, seq2)
        
    elif dist_method == "td":
        distance = tamura_distance(seq1, seq2)
    
    #elif dist_method == "tnd":
    #    distance = tn_distance(seq1, seq2)
           
    return 1 - distance




