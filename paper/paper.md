---
title: 'Recan: Python tool for analysis of recombination events in viral genomes'
tags:
  - Python
  - virology
  - recombination
authors:
- name: Yuriy Babin
  orcid: 0000-0002-7524-5921
  affiliation: "1"
affiliations: 
- name: National Medical Research Center for Tuberculosis and Infectious Diseases, Moscow, Russia
  index: 1
date: 29 November 2019
bibliography: paper.bib
---


# Summary

Recombination drives virus evolution in response to selective forces in a host environment and adaptation to new abiotic factors [@Perez-Losada2015].
Gaining insights into recombination events is important for a better understanding of viral biology. Analysis of recombination events can be performed through construction and exploration of similarity plots based on  genetic distances between nucleotide sequences.
Python package named "recan" (recombination analyzer) provides the means to construct genetic distance plots and explore them interactively. The package has been designed to operate with Jupyter notebooks. Compared to the previously designed desktop software [@Lole1999; @Etherington2005a; @Martin2015] recan has the ability to insert or delete sequences from the output without reconstructing the plots and recalculating the distance values. Finally, recan enables simultaneous analysis of several datasets in a single session. Recan is based on Biopython, Pandas, and Plotly libraries. The package requires a sequence alignment in fasta format as an input. The user can adjust the sliding window size, the window shift, method of distance calculation, sequence of interest (a sequence where breakpoints occur), and the length of alignment region which will be included into the distance calculation. The two methods of genetic distance calculations implemented in recan are the pairwise and Kimura 2-parameter models. The distance data can be saved in csv or excel file, or directly used to reconstruct the plot in the Jupyter notebook using the plotting library to obtain a final report. 

# Testing and verification
To test the package, we used four previously reported recombinant viral genomes representing different genuses: human immunodeficiency virus (HIV) [@Liitsola2000a], hepatitis C virus (HCV) [@Smith2014], norovirus [@Jiang1999], and lumpy skin disease virus (LSDV) [@Sprygin2018].
Each dataset included a recombinant virus sequence, its putative parental sequences and a set of sequences of the same virus closely related to the recombinant virus. HIV, HCV and Norovirus sequences were aligned using ClustalW [@Larkin2007], and LSDV genomes were aligned using MAFFT [@Katoh2002] as part of Ugene software [@Okonechnikov2012].
The HIV alignment contained twenty five 3135 bp sequences; the HCV alignment contained twenty three 9431 bp sequences; the norovirus alignment included nineteen 3366 bp sequences, and the LSDV alignment had a total of 150511 bp sequences. The resulting `simgen` method execution time with the default window size and shift parameters was the following: 437 ms ± 7.74 ms for HIV, 579 ms ± 58.7 ms for Norovirus, 648 ms ± 44.2 ms for HCV, and 3.55 s ± 239 ms for LSDV dataset. Time execution test was performed using a desktop PC with 4 CPU cores and 4 Gb RAM. LSDV has one of the largest genomes of all viruses (about 150 000 bp). Ultimately, recan can potentially be used to identify and analyze recombination events in a large subset of sequences regardless of the length of the viral genome. The distance plots with recombination events detected by recan are shown in Figures 1-4.

# Availability and implementation
Recan is supported on Linux and Windows. The package can be installed by `pip` Python package manager using `pip install recan` command. The source code, guide and datasets are available on the GitHub repository (https://github.com/babinyurii/recan). 

![](https://raw.githubusercontent.com/babinyurii/recan/master/paper_plots/hiv_rec_kal153.png)
_Figure 1. HIV recombinant strain AF193276 between sequences AF193275 and AF193278._


![](https://raw.githubusercontent.com/babinyurii/recan/master/paper_plots/hcv_2k_1b_rec.png)
_Figure 2. HCV intergenotype recombinant 2k/1b._


![](https://raw.githubusercontent.com/babinyurii/recan/master/paper_plots/norovirus_rec.png)
_Figure 3. Norovirus recombinant AF190817 between parental sequences U22498 and X86557._


![](https://raw.githubusercontent.com/babinyurii/recan/master/paper_plots/lsdv_rec_sar.png)
_Figure 4. LSDV recombinant vaccine-like strain LSDV RUSSIA/Saratov/2017 between sequences AF193275 and KY829023._

# Acknowledgement 
The author thanks Alexander Sprygin for editing the manuscript and providing LSDV data.



# References
