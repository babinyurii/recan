""" Simgen realizes the interface to manipulate alignment in Jupyter notebook,
and to explore recombination events using similarity plots
"""

from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
#import matplotlib.pyplot as plt
import pandas as pd
from plotly.offline import init_notebook_mode, iplot
import plotly.graph_objs as go



class Simgen(MultipleSeqAlignment):
    """
    A class used to represent an alignment object
    ---
    Methods
    -------
    simgen(pot_rec, window=500, shift=250, region=False, dist="pdist")
        outputs similarity plots using slidign window approach
    get_data(df=True)
        returns distance data as a pandas DataFrame object or a dictionary of lists 
        with key as a sequence id and distances in corresponding lists
    get_info()
        outputs information about the alignment, index (which is the row number), 
        sequence names, and alignment length
    save_data(path=False, out="csv")
        saves the data spreadsheet as a csv or excel file
    
    """
      
    def __init__(self, path):
        """
        initializing Simgen
        ---
        Parameters
        ----------
        path: str
            path to the fasta file containing alignment
        """
        from Bio import AlignIO
        recs_prepared = (x for x in AlignIO.read(path, "fasta")) # it works without it, but i believe it's formally right
        super(Simgen, self).__init__(recs_prepared)  # you need explicitly call the __init__ of the upperclass
        self._distance = {}  # empty, until call of the simgen function
        self._align = None  # current slice or whole MultipleSeqRecord for plotting
        self._ticks = None  # ticks for plot
        
    
    def _draw_simplot(self):
        """draws similarity plot using plotly"""
        init_notebook_mode()
        data = []
        for key in self._distance.keys():
            trace = go.Scatter(y=self._distance[key], x=self._ticks, name=key)
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
        
    #TODO
    # tick label merge
    def _draw_simplot_mpl(self):
        """draws similartiy plot using  matplotlib """
        
        data = list(self._distance.values())
        ticks = self._ticks[1:]
        #labels = list(self._distance.keys())
        #print(data, ticks, labels, sep="\n")
        plt.figure(figsize=(15, 8))
        for i in data:
            plt.plot(range(1, len(i) + 1), i)
        plt.xticks(list(range(1, len(ticks))), ticks, rotation='vertical')
        plt.show()
    
    
    def _get_x_labels(self, left_border, right_border, shift):
        """creates tick labels"""

        tick_container = []      
        tick_container.append(left_border)

        while tick_container[-1] < right_border:
            tick_container.append(tick_container[-1] + shift)
            if tick_container[-1] > right_border:
                tick_container[-1] = right_border
        #self._ticks = tick_container
        return tick_container
    
    def _move_window(self, window, pot_rec, shift, dist):
        """moves window"""
        distance_data = {}
        parents = list(range(0, len(self._align)))
        parents.remove(pot_rec)
        align_length = len(self._align[0, :])        
        
        for par in parents:
            dist_container = []
            start = 0
            finish = shift

            while start < align_length:
                seq1 = self._align[pot_rec, start:finish].seq # here is a potential recombinant sequence slice
                seq2 = self._align[par, start:finish].seq  # here's a parent's slice
                
                if dist == "pdist":
                    distance = self._pdistance(seq1, seq2)
                    dist_container.append(distance) #calculate pdistance, append to container
                elif dist == "k2p":
                    distance = self._K2Pdistance(seq1, seq2)
                    dist_container.append(distance)
                    #dist_container.append(self._K2Pdistance(seq1, seq2)) #calculate pdistance, append to container
                start += shift
                finish = start + window
            distance_data[self._align[par].id] = dist_container
    
        #self._distance = distance_data  
        return distance_data
    
    def _K2Pdistance(self, seq1, seq2):
        """
        Kimura 2-Parameter distance = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
        where:
        p = transition frequency
        q = transversion frequency
        """
        from math import log, sqrt
        pairs = []
        for x in zip(seq1, seq2):
            if '-' not in x: 
                pairs.append(x)
        ts_count=0
        tv_count=0
        length = len(pairs)

        transitions = [ "AG", "GA", "CT", "TC"]
        transversions = [ "AC", "CA", "AT", "TA",
                          "GC", "CG", "GT", "TG" ]

        for (x, y) in pairs:
            if x + y in transitions: 
                ts_count += 1 
            elif x + y in transversions: 
                tv_count += 1

        p = float(ts_count) / length
        q = float(tv_count) / length
        try: 
            d = -0.5 * log( (1 - 2*p - q) * sqrt( 1 - 2*q ) )
        except ValueError: 
            print ("Tried to take log of a negative number")
            return None
        return 1 - d
   
        
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
        
    def _get_collect_sliced_left_right_borders(self, region=False):
        
        if region:
            assert region[0] < region[1], "the value of the first nucleotide position should be less than the second one"
            collect_sliced = []
            for rec in self._records:  # access to seq of the SeqRecord obj inside MultipleSeqAlignment
                sliced_seq = rec.seq[region[0]:region[1]]
                collect_sliced.append(SeqRecord(sliced_seq, id=rec.id, name=rec.name, description=rec.description))
            left_border = region[0]   # border for the first tick
            right_border = region[1]  # if region, 'right_border' is actual position

        else:
            collect_sliced = []
            for rec in self._records:  # access to seq of the SeqRecord obj inside MultipleSeqAlignment
                sliced_seq = rec.seq[:]
                collect_sliced.append(SeqRecord(sliced_seq, id=rec.id, name=rec.name, description=rec.description))
            
            left_border = 1  # border for the first tick
            right_border = self.get_alignment_length()
            
        return collect_sliced, left_border, right_border
            
            
    def simgen(self, pot_rec, window=500, shift=250, region=False, dist='pdist'):
        """slices the alignment, collects the distance data, outputs the plot

        Parameters:
        -----------
        pot_rec: int
            the index of the sequence under study. use get_info() to find out.
            it is the index of the row in the alignment: starts with 0,
            like the 'x' dimension in the numpy array
        window: int
            sliding window size. 500 by default
        shift: int
            the step window slides downstream the alignment. 250 by default
        region: a tuple or a list of two integers or False by default
            the region of the alignment to analyze: the start
            and the end nucleotide positions. 
            by default takes the whole alignment length
        dist: str
            'pdist' or 'k2p' 
            pairwise or Kimura methods to calculate distance
       
            """

        assert window >=1, "window  parameter can't be a negative or zero"
        assert shift >= 1, "shift parameter can't be a negative or zero" 

        collect_sliced, left_border, right_border = self._get_collect_sliced_left_right_borders(region)
        self._align = MultipleSeqAlignment(collect_sliced)
        # creating tick labels for the plot
        
        #self._get_x_labels(left_border, right_border, shift)
        self._ticks = self._get_x_labels(left_border, right_border, shift)
        
        # calculating pairwise distance
        
        #self._move_window(window, pot_rec, shift, dist)
        self._distance = self._move_window(window, pot_rec, shift, dist)
        
        self._draw_simplot()
        #TODO
        # inter=True, when tick merge will be solved
        # add this parameter to choose between plotly and mpl
        #if inter:
        #    self._draw_simplot()
        #else:
        #    self._draw_simplot_mpl()
            

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
            return pd.DataFrame(data=self._distance, index=self._ticks[1:]).T
        else:
            return self._ticks[1:], self._distance
    
    
    def get_info(self):
        """outputs information about the alignment: 
        index (which is the row number), 
        sequence names, and alignment lengths"""
        
        print("index:", "sequence id:", sep="\t")
        for counter, value in enumerate(self):
            print(counter, value.id, sep="\t")
        print("alignment length: ", self.get_alignment_length())
        #print(self)
    
        
    def save_data(self, path=False, out="csv", out_name="distance_data"):
        """saves the data spreadsheet as a csv or excel file
        Parameters
        ---------
        path: str
            output destination
        out: str
            output file format: "csv" or "excel"
        out_name: str
            output file name
        """
        
        df = pd.DataFrame(data=self._distance, index=self._ticks[1:]).T
        if path:
            if out == "csv":
                df.to_csv(out_name + ".csv")
            elif out == "excel":
                writer = pd.ExcelWriter(out_name + ".xlsx")
                df.to_excel(path)
            else:
                print("invalid output file")
        else:
            if out == "csv":
                df.to_csv(out_name + ".csv")
            elif out == "excel":
                writer = pd.ExcelWriter(out_name + ".xlsx")
                df.to_excel(writer)
            else:
                print("invalid output file format")
        