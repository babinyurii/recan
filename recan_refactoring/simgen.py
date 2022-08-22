# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 12:56:47 2022

@author: yuriy
"""
import plotly.graph_objs as go
from plotly.offline import init_notebook_mode, iplot
from rolling_window import RollingWindowOnAlignment


class Simgen():
    
    def __init__(self, in_file):
        
        self.alignment_roll_window = RollingWindowOnAlignment(in_file)
        self.align_sliced = None
        self.ticks_for_x_axis = None
        self.distance_data = {}
        self.pot_rec_index = None
        self.pot_rec_id = None
        
        
    def _plot_distance_by_plotly(self):
                
        """draws similarity plot using plotly"""
        init_notebook_mode()
        data = []
        for key in self.distance_data.keys():
            trace = go.Scatter(y=self.distance_data[key], x=self.ticks_for_x_axis, name=key)
            data.append(trace)
        
        layout = go.Layout(
            title = self.pot_rec_id,
            xaxis=dict(
                title="nucleotide position"),
            yaxis=dict(
                title="sequence identity"),
            legend=dict(x=-0.1, y=1.5, orientation="h"))
            #legend=dict(x=-0.1, y=1.5))
            
        fig = go.Figure(data=data, layout=layout)
        iplot(fig)
        
    
    def _pdistance(self, seq1, seq2):
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
        
    def _get_ticks_for_x_axis(self):
        
        ticks = []
        for start_stop_nuc_index in self.align_sliced.keys():
            #ticks.append(key[0])
            #ticks.append(key[1])
            ticks.append((start_stop_nuc_index[0] + start_stop_nuc_index[1]) / 2)
        
        return ticks
    
    def _get_pot_rec_id(self):
        
        for slice in self.align_sliced.values():
            slice_0 = slice
            break
        
        self.pot_rec_id = slice_0[self.pot_rec_index].id
        
    
    
    def _prepare_distance_data(self):
        
        for slice in self.align_sliced.values():
            slice_0 = slice
            break
        
        for seq in slice_0:
            if seq.id == self.pot_rec_id:
                continue
            self.distance_data[seq.id] = []
            
    
    def _calculate_distance(self):
        
        #print(self.align_sliced)
        
        for window_borders, alignment_slice in self.align_sliced.items():
            for seq in alignment_slice:
                if seq.id == self.pot_rec_id:
                    continue
            
                distance = self._pdistance(alignment_slice[self.pot_rec_index], seq.seq)
                self.distance_data[seq.id].append(distance)
        
            
        
    
    
    def simgen(self, pot_rec, window, shift, region=False, dist="pdist"):
        
        
        
        self.align_sliced = self.alignment_roll_window.roll_window_along_alignment(window_len=window, 
                                                                          window_step=shift)
        
        self.pot_rec_index = pot_rec
        self._get_pot_rec_id()
        
        self._get_ticks_for_x_axis()
        
        self._prepare_distance_data()
        
        self._calculate_distance()
        
        self._plot_distance_by_plotly()
        
        
        
        
        
        
       
    

    
    
    