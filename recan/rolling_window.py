from Bio import AlignIO

class RollingWindowOnAlignment():
    """
    alignment obj as the biopython multiple alignment
    has two sliding window methods that slice the alignment into sections
    """
    
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
            
            ###############################################
            # redo
            #new_window_end = window_end + window_step
            if window_end + window_step < self.align.get_alignment_length():
                window_end += window_step
            else:
                window_end = self.align.get_alignment_length()
            # redo
            ############################################3
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
            
            if window_end + window_step < region[1]:
                window_end += window_step
            else:
                window_end = region[1]
            

            window_counter += 1
            

        return sliced_alignment
    
    
    
    
    