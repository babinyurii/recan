# -*- coding: utf-8 -*-

import plotly.graph_objs as go
import pandas as pd
from plotly.offline import init_notebook_mode, iplot
from .rolling_window import RollingWindowOnAlignment
from .calc_pairwise_distance import calc_pairwise_distance


class Simgen():
    
    def __init__(self, in_file):
        
        self.alignment_roll_window = RollingWindowOnAlignment(in_file)
        self.align_sliced = None
        self.ticks_for_x_axis = []
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
            xaxis=dict(
                title="nucleotide position"),
            yaxis=dict(
                title="sequence identity"),
            legend=dict(x=-0.1, y=1.5, orientation="h"))
            #legend=dict(x=-0.1, y=1.5))
            
        fig = go.Figure(data=data, layout=layout)
        iplot(fig)
        print("potential recombinant: ", self.pot_rec_id)
        
        
    def _get_ticks_for_x_axis(self):
        """
        prepares ticks for x axis
        the ticks correspond to the alignment nucleotide positions
        the ticks for x axis are in the middle of the rolling window
        """
        
        self.ticks_for_x_axis.clear()
        
        for start_stop_nuc_index in self.align_sliced.keys():
            # we take a nucleotide position in the middle of the window
            # it'll be a point we plot
            middle_nuc = (start_stop_nuc_index[1] - start_stop_nuc_index[0]) / 2
            self.ticks_for_x_axis.append(start_stop_nuc_index[0] + middle_nuc)
            
    
    def _get_pot_rec_id(self):
        """ get id of the potential recombinant from
        the alignment slices
        """
        
        for slice in self.align_sliced.values():
            slice_0 = slice
            break
        
        self.pot_rec_id = slice_0[self.pot_rec_index].id
        
    
    
    def _prepare_distance_data(self):
        
        """ makes dictionary for distance data collection
        keys are sequences names (ids from the alignment)
        values are lists of distances to the potential recombinant
        """
        self.distance_data = {}
        
        for slice in self.align_sliced.values():
            slice_0 = slice
            break
        
        for seq in slice_0:
            if seq.id == self.pot_rec_id:
                continue
            self.distance_data[seq.id] = []
            
    
    
    def simgen(self, pot_rec, window, shift, region=False, dist="pdist"):
        """
        Parameters
        ----------
        pot_rec : int
            index of the potential recombinant in the alignment.
        window : int
            sliding window size.
        shift : int
            sliding window shift along the alignment.
        region : list or tuple, optional
             The default is False. start and end of the region
             to analyze, f.e. (1000, 3000)
        dist : str, optional
            the distance calculation method. default is "pdist". 
            
            available methods:
            pdist - pairwise distance (default)
            jcd - Jukes-Cantor distance
            k2p -  Kimura 2-parameter distance
            td - Tamura distance
            
        Returns
        -------
        None.

        """
        
        if region:
            assert region[0] < region[1], "start of the region must be less than the region end"
            self.align_sliced = self.alignment_roll_window.roll_window_along_alignment_region(window_len=window, 
                                                                          window_step=shift,
                                                                          region=region)
        else:
            self.align_sliced = self.alignment_roll_window.roll_window_along_alignment(window_len=window, 
                                                                          window_step=shift)
        self.pot_rec_index = pot_rec
        self._get_pot_rec_id()
        self._get_ticks_for_x_axis()
        self._prepare_distance_data()
        
        
        for window_borders, alignment_slice in self.align_sliced.items():
            for seq in alignment_slice:
                if seq.id == self.pot_rec_id:
                    continue
                seq1 = alignment_slice[self.pot_rec_index].seq
                seq2 = seq.seq
                distance = calc_pairwise_distance(seq1=seq1, seq2=seq2, 
                                                 dist_method=dist) 
                self.distance_data[seq.id].append(distance)

        
        self._plot_distance_by_plotly()
        
        
        
        
    def save_data(self, path=False, out="csv", out_name="distance_data",
                  data_cols="plot_ticks"):
        """saves the data spreadsheet as a csv file
        Parameters
        ---------
        path: str
            output destination
        out: str
            output file format: "csv" or "excel"
        out_name: str
            output file name
        columns: str
            denoting columns in the spreadsheet: "plot_ticks" by default,
            or "window_pos" for sliding window start and end nucleotide 
            position in alignment
        """
        
        if data_cols == "plot_ticks":
            columns = self.ticks_for_x_axis
        elif data_cols == "window_pos":
            columns = [x for x in self.align_sliced.keys()]
       
        
        df = pd.DataFrame.from_dict(self.distance_data, 
                                    orient='index',
                                    columns=columns)
       
        if path:
            if out == "csv":
                df.to_csv(out_name + ".csv")
            else:
                print("invalid output file")
        else:
            if out == "csv":
                df.to_csv(out_name + ".csv")
            else:
                print("invalid output file format")
        
   
        
        
    def get_data(self, df=True):
        """returns distance data
        Parameters
        ---------
        df: bool
            True: returns pandas DataFrame object
            False: returns a dictionary where keys are the sequence ids and 
            values are distance data
        
        """
        if df:
            return pd.DataFrame(data=self.distance_data, index=self.ticks_for_x_axis).T
        else:
            return self.ticks_for_x_axis, self.distance_data
        
        
    def get_info(self):
        """outputs information about the alignment: 
        index (which is the row number), 
        sequence names, and alignment length"""
        
        print("index:", "sequence id:", sep="\t")
        for seq_index, seq in enumerate(self.alignment_roll_window.align):
            print(seq_index, seq.id, sep="\t")
        print("alignment length: ", self.alignment_roll_window.align.get_alignment_length())
      
        
        
    