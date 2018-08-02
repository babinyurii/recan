# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 17:59:26 2018
@author: babin
"""
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import pandas as pd
import plotly.graph_objs as go
import plotly.plotly as py
from Bio import AlignIO
init_notebook_mode(connected=True)


def _pdistance(seq1, seq2):
    """calculates pairwise distance between two sequences"""
    p = 0
    pairs = []
    for x in zip(seq1, seq2):
        if '-' not in x:
            pairs.append(x)
    for (x, y) in pairs:
        if x != y:
            p += 1
    length = len(pairs)
    return float(1 - p / length)  # '1 - p' to take plot 'upside down'


def _draw_simplot(distance_data, tick_container):
    """draws similarity plot"""

    data = []
    for key in distance_data.keys():
        trace = go.Scatter(y=distance_data[key], x=tick_container, name=key)
        data.append(trace)

    layout = go.Layout(
        title="similarity plot",
        xaxis=dict(
            title="nucleotide position"),
        yaxis=dict(
            title="sequence identity"),
        legend=dict(x=-0.1, y=1.5))

    fig = go.Figure(data=data, layout=layout)
    iplot(fig)


def simgen(align, pot_rec, window=500, shift=100, region=False, return_data=False):
    """slices the alignment, collects the distance data

    Parameters:
    -----------
    align: the location of the fasta alignment
    pot_rec: int
        the number of the sequence under study, starts with 0,
        like the 'x' dimension in the numpy array
        window: int
        sliding window size
    shift: int
        the step window slides downstream the alignment
    region: a tuple or a list of two integers
        the region of the alignment to analyze. the start
        and the end nucleotide positions
    return_data: bool, optional
        return the data in pandas DataFrame
        """

    if region:
        align = AlignIO.read(align, "fasta")
        align = align[:, region[0] : region[1]]
        left_border = region[0]   # border for the first tick
    else:
        align = AlignIO.read(align, "fasta")
        left_border = 1  # border for the first tick

    distance_data = {}
    parents = list(range(0, len(align)))
    parents.remove(pot_rec)
    align_length = len(align[0, :])

    # creating tick labels for the plot
    tick_container = []
    tick_container.append(left_border)
    while tick_container[-1] < align_length:
        tick_container.append(tick_container[-1] + shift)

    for par in parents:
        dist_container = []
        start = 0
        finish = shift

        while start < align_length:
            # potential recombinant slice
            seq1 = align[pot_rec, start:finish].seq
            seq2 = align[par, start:finish].seq  # a parent's slice
            dist_container.append(_pdistance(seq1, seq2))
            start += shift
            finish = start + window

        distance_data[align[par].id] = dist_container

    if return_data:
        # [1:] to map data to index
        data = pd.DataFrame(data=distance_data, index=tick_container[1:])
        return data
    else:
        _draw_simplot(distance_data, tick_container)


simgen("C:\\for_recan_lsdv_1.fasta", pot_rec=2)