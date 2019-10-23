# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 15:57:30 2019

@author: babin
"""
from Bio.Align import MultipleSeqAlignment

import unittest
from test_results import dist_whole_align_ref, dist_win_250_shift_100_ref, dist_whole_align_def_params_k2p
import sys
sys.path.append("..")
from recan.simgen import Simgen


class DistTestCase(unittest.TestCase):
    
    
    def setUp(self):
        
        self.sequences = [("AAAAAGGGGG", "AAAAAGGGGG"),
                            ("AAAAAAAAAT", "AAAAAAAAAA"),
                            ("AAAAATTTTT", "AAAAAAAAAA"),
                            ("ATTTTTTTTT", "AGGGGGGGGGG")]
    
        self.pdistances = [1.0, 0.9, 0.5, 0.1 ]
        
        self.ids = ["AB048704.1_genotype_C_", "AB033555.1_Ba", "AB010291.1_Bj"]
        
        self.default_shift = 250
        self.default_window = 500
        self.ticks_whole_align_ref = [1, 251, 501, 751, 1001, 1251, 1501, 1751, 2001, 2251, 2501, 2751, 3001, 3215]
        self.ticks_region_ref = [250, 500, 700]
        self.pot_rec = 1
        self.dist_reg_ref = {'AB048704.1_genotype_C_': [0.9319999999999999, 0.945], 'AB010291.1_Bj': [0.992, 0.975]}
        
        
        # whole alignment slice
        self.sim_obj = Simgen("./hbv_C_Bj_Ba.fasta")
        self.collect_sliced, self.left_border, self.right_border = self.sim_obj._get_collect_sliced_left_right_borders(region=False)
        self.whole_align_obj = MultipleSeqAlignment(self.collect_sliced)
        self.ticks_whole_align = self.sim_obj._get_x_labels(self.left_border, self.right_border, self.default_shift)
        self.sim_obj._align = MultipleSeqAlignment(self.collect_sliced) # store alignment in _align attr of the sim_obj, otherwise move window will throw error
        self.distance_whole_align = self.sim_obj._move_window(self.default_window, self.pot_rec, self.default_shift, "pdist")
        # changin window and shift params
        self.dist_whole_align_win_250_shift_100 = self.sim_obj._move_window(window=250, pot_rec=self.pot_rec, shift=100, dist="pdist")
        #k2p distance with the same params
        self.distance_whole_align_k2p = self.sim_obj._move_window(self.default_window, self.pot_rec, self.default_shift, "k2p")
        
        # region slice
        self.sim_obj_reg = Simgen("./hbv_C_Bj_Ba.fasta")
        self.collect_sliced_reg, self.left_border_reg, self.right_border_reg = self.sim_obj_reg._get_collect_sliced_left_right_borders(region=(250, 700))
        self.sliced_align_obj = MultipleSeqAlignment(self.collect_sliced_reg)
        self.ticks_region = self.sim_obj_reg._get_x_labels(self.left_border_reg, self.right_border_reg, self.default_shift)
        self.sim_obj_reg._align = MultipleSeqAlignment(self.collect_sliced_reg)
        self.dist_reg = self.sim_obj_reg._move_window(self.default_window, self.pot_rec, self.default_shift, "pdist")
        
        
        
    def test_get_collect_sliced_left_right_borders_whole_alignment_ids(self):
        for counter, value in enumerate(self.collect_sliced):
            self.assertEqual(value.id, self.ids[counter])
            
    def test_get_collect_sliced_left_right_borders_region_ids(self):
        for counter, value in enumerate(self.collect_sliced_reg):
            self.assertEqual(value.id, self.ids[counter])
            
   
    def test_get_collect_sliced_left_right_borders_whole_alignment_len(self):
        for i in self.collect_sliced:
            self.assertEqual(len(i.seq), 3215)
            
    def test_get_collect_sliced_left_right_borders_region_len(self):
        for i in self.collect_sliced_reg:
            self.assertEqual(len(i.seq), 450)
            
    
    def test_get_collect_sliced_left_right_borders_whole_alignment(self):
        self.assertEqual(self.left_border, 1)
        self.assertEqual(self.right_border, 3215)
        
    
    def test_get_collect_sliced_left_right_borders_region(self):
        self.assertEqual(self.left_border_reg, 250)
        self.assertEqual(self.right_border_reg, 700)
        
    def test_get_x_labels_whole_align_ticks(self):
        for counter, value in enumerate(self.ticks_whole_align):
            self.assertEqual(value, self.ticks_whole_align_ref[counter])
            
    def test_get_x_labels_region_ticks(self):
        for counter, value in enumerate(self.ticks_region):
            self.assertEqual(value, self.ticks_region_ref[counter])
            
        
    def test_move_window_whole_align_default_params(self):
        for key in self.distance_whole_align.keys():
            for counter, value in enumerate(self.distance_whole_align[key]):
                self.assertEqual(value, dist_whole_align_ref[key][counter] )
                
    def test_move_window_whole_align_win_250_shift_100(self):
        for key in self.dist_whole_align_win_250_shift_100.keys():
            for counter, value in enumerate(self.dist_whole_align_win_250_shift_100[key]):
                self.assertEqual(value, dist_win_250_shift_100_ref[key][counter] )
    
    def test_move_window_whole_align_def_params_k2p(self):
        for key in self.distance_whole_align_k2p.keys():
            for counter, value in enumerate(self.distance_whole_align_k2p[key]):
                self.assertEqual(value, dist_whole_align_def_params_k2p[key][counter] )
                
        
    def test_move_window_region(self):
        for key in self.dist_reg.keys():
            for counter, value in enumerate(self.dist_reg[key]):
                self.assertEqual(value, self.dist_reg_ref[key][counter])
            

    def test_pdistance_method(self):
        for counter, value in enumerate(self.pdistances):
            dist = self.sim_obj._pdistance(self.sequences[counter][0], self.sequences[counter][1])
            self.assertEqual(round(dist, 1), value)
    


if __name__ == '__main__':
    unittest.main()
        