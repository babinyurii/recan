from Bio import AlignIO


def roll_window_along_alignment(in_file, window_len, window_step):
    
    align = AlignIO.read(in_file, "fasta")
    window_start = 0
    window_end = window_len
    window_step = window_step

    window_counter = 0
    sliced_alignment = {}
    while window_start < align.get_alignment_length():
        sliced_alignment[(window_start, window_end)] = align[:, window_start:window_end]
        window_start += window_step
        window_end += window_step

        window_counter += 1

    return sliced_alignment


def roll_window_along_alignment_region(in_file, window_len, window_step, region):
    
    align = AlignIO.read(in_file, "fasta")
    window_start = region[0]
    window_end = region[0] + window_len
    window_step = window_step

    window_counter = 0
    sliced_alignment = {}
    while window_start < region[1]:
        sliced_alignment[(window_start, window_end)] = align[:, window_start:window_end]
        window_start += window_step
        window_end += window_step

        window_counter += 1

    return sliced_alignment
    