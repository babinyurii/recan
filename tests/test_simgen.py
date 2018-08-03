# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 11:07:59 2018

@author: babin
"""

from recan.simgen import _pdistance, simgen



def test_pdistance():
    """
    if sequences are totally equal, return 1
    if totally different, return 0
    """
    assert _pdistance("A", "A") == 1
    assert _pdistance("A", "T") == 0
    assert _pdistance("AAAAATTTTT", "GGGGGTTTTT") == 0.5
    assert _pdistance("AAAAA-----", "TTTTT-----") == 0
    assert _pdistance("----------", "----------") == None
    assert _pdistance("A-A-A-A-A-", "-A-A-A-A-A") == None


def test_simgen():
    a, b = simgen("C:\\recan\\data\\equal.fasta", pot_rec=0, draw=False).values()
    assert a, b == (1, 1)
    a, b = simgen("C:\\recan\\data\\empty.fasta", pot_rec=0, draw=False).values()
    assert a, b == (None, None)
    a, b = simgen("C:\\recan\\data\\half.fasta", pot_rec=0, draw=False).values()
    assert a, b == (0.5, 0.5)
    a, b = simgen("C:\\recan\\data\\half_and_equal.fasta", pot_rec=0, draw=False).values()
    assert a, b == (0.5, 1)

    

    