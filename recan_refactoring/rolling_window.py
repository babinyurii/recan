from Bio import AlignIO

class RollingWindowOnAlignment():
    
    def __init__(self, in_file):
        self.align = AlignIO.read(in_file, "fasta")
        
    def roll_window_along_alignment(self, window_len, window_step):
    
        
        window_start = 0
        window_end = window_len
        window_step = window_step

        window_counter = 0
        sliced_alignment = {}
        while window_start < self.align.get_alignment_length():
            sliced_alignment[(window_start, window_end)] = self.align[:, window_start:window_end]
            window_start += window_step
            window_end += window_step

            window_counter += 1

        return sliced_alignment 
    
    def roll_window_along_alignment_region(self, window_len, window_step, region):
    
        
        window_start = region[0]
        window_end = region[0] + window_len
        window_step = window_step

        window_counter = 0
        sliced_alignment = {}
        while window_start < region[1]:
            sliced_alignment[(window_start, window_end)] = self.align[:, window_start:window_end]
            window_start += window_step
            window_end += window_step

            window_counter += 1

        return sliced_alignment
    
    
    
    
    