# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 12:26:58 2022

@author: yuriy
"""

def _estimate_nucleotide_frequencies(seq):
    seq = seq.replace('-','').upper()
    A = seq.count('A')
    C = seq.count('C')
    G = seq.count('G')
    T = seq.count('T')
    length = float(len(seq))
    return [ x/length for x in [A,C,G,T] ]
    
    
def _p_distance(seq1, seq2):
    """calculates 1 - pairwise distance between two sequences"""
    p = 0
    pairs = []
    for x in zip(seq1, seq2):
        if '-' not in x:
            pairs.append(x)
    for (x, y) in pairs:
        if x != y:
            p += 1
    length = len(pairs)
    #assert length > 0, "AssertionError: perhaps your alignment contains only or too many gaps"
    try:
        dist = float(1 - p / length)  # '1 - p' to take plot 'upside down'
        return dist
    except ZeroDivisionError as e:
        #print(e, ": perhaps your alignment contains only gaps")
        pass
    
    
def _jc_distance(p_distance):
    """ 
    distance = -b log(1 - p / b)
    where:
    b = 3/4
    and p = p-distance, i.e. uncorrected distance between seq1 and seq2
    """
    from math import log
    b = 0.75
    #p = _p_distance(seq1,seq2)
    try: 
        d = -b * log( 1 - p_distance / b )
    except ValueError: 
        print ("Tried to take log of a negative number")
        return None
    return 1 - d


       
def _k2p_distance(seq1,seq2):
    """
    Kimura 2-Parameter distance = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
    where:
    p = transition frequency
    q = transversion frequency
    """
    from math import log, sqrt
    pairs = []

    #collect ungapped pairs
    for x in zip(seq1,seq2):
        if '-' not in x: 
            pairs.append(x)
        
    ts_count=0
    tv_count=0
    length = len(pairs)
    
    transitions = [ "AG", "GA", "CT", "TC"]
    transversions = [ "AC", "CA", "AT", "TA",
                      "GC", "CG", "GT", "TG" ]

    for (x,y) in pairs:
        if x+y in transitions: 
            ts_count += 1 
        elif x+y in transversions: 
            tv_count += 1
    
    p = float(ts_count) / length
    q = float(tv_count) / length
    
    try: 
        d = -0.5 * log((1 - 2 * p - q) * sqrt( 1 - 2 * q ))
    except ValueError: 
        print ("Tried to take log of a negative number")
        return None
    return 1 - d 

def _tn_distance(seq1, seq2):
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
    G = _estimate_nucleotide_frequencies(seq1 + seq2)
    p = _p_distance(seq1,seq2)
    pairs = []
    h = 0

    #collect ungapped pairs
    for x in zip(seq1,seq2):
        if '-' not in x: pairs.append(x)
       
    #pair frequencies are calculated for AC, AG, AT, CG, CT, GT (and reverse order)
    for i in range(len(ns)-1):
        for j in range(i+1,len(ns)):
            if i != j: paircount = pairs.count( (ns[i], ns[j]) ) + pairs.count( (ns[j], ns[i]) )
            Xij_sq = (float(paircount)/len(pairs))**2
            GiGj = G[i]*G[j]
            h += 0.5*Xij_sq/GiGj  #h value used to calculate b
    
    b = 0.5*(1-sum([x**2 for x in G])+p**2/h)
    try: d = -b * log(1 - p/b)
    except ValueError: 
        print ("Tried to take log of a negative number")
        return None
    return d      

def _tamura_distance(seq1,seq2):
    """
    Tamura distance = -C log( 1 - P/C - Q ) - 0.5( 1 - C )log( 1 - 2Q )
    where:
    P = transition frequency
    Q = transversion frequency
    C = GC1 + GC2 - 2 * GC1 * GC2
    GC1 = GC-content of sequence 1
    GC2 = GC-coontent of sequence 2
    """
    from math import log
    pairs = []
    
    #collect ungapped pairs
    for x in zip(seq1,seq2):
        if '-' not in x: pairs.append(x)
        
    ts_count=0
    tv_count=0
    length = len(pairs)
    
    transitions = [ "AG", "GA", "CT", "TC"]
    transversions = [ "AC", "CA", "AT", "TA",
                      "GC", "CG", "GT", "TG" ]

    for (x,y) in pairs:
        if x+y in transitions: ts_count += 1 
        elif x+y in transversions: tv_count += 1
    
    p = float(ts_count) / length
    q = float(tv_count) / length
    gc1 = sum(_estimate_nucleotide_frequencies(seq1)[1:3])
    gc2 = sum(_estimate_nucleotide_frequencies(seq2)[1:3])
    c = gc1 + gc2 - 2 * gc1 * gc2

    try: d = -c * log( 1 - p/c - q) - 0.5 * ( 1 - c ) * log ( 1 - 2*q )
    except ValueError: 
        print ("Tried to take log of a negative number")
        return None
    return 1 - d


def calc_pairwise_distance(seq1, seq2, dist_method):
    
    if dist_method == "pdist":
        distance = _p_distance(seq1, seq2)
        
    elif dist_method == "jcd":
        p_distance = _p_distance(seq1, seq2)
        distance = _jc_distance(p_distance)
        
    elif dist_method == "k2p":
        distance = _k2p_distance(seq1, seq2)
        
    elif dist_method == "tnd":
        distance = _tn_distance(seq1, seq2)
       
    elif dist_method == "td":
        distance = _tamura_distance(seq1, seq2)
       
        
    
    
    
    
    return distance