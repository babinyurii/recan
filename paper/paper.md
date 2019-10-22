---
title: 'Recan: A Python package for recombination events analysis by distance plotting'
tags:
  - Python
  - biology
authors:
  - name: Yuriy Babin
    orcid: 0000-0002-7524-5921
    affiliation: 
affiliations:
date:  November 2019
bibliography: paper.bib
---
# Summary

The similarity plot based on genetic distance between nucleotide sequences is intuitive and simple way to explore and discover recombination events. We have developed Python package named ‘recan’ (recombination analyzer) which allows to construct genetic distance plots and explore them interactively. 
Recan is based upon Biopython, Pandas, Matplotlib, and Plotly Python libraries. The package is intended to be used via Jupyter notebook, which makes it flexible and easy to use.
Recan code is organized using OOP paradigm. Simgen class is inherited from the MultipleSequenceAlignment class of Biopython library. The package requires sequence alignment in the fasta format as an input. User interface includes importing Simgen class, initializing an instance of this class, and applying one of its methods. The main method constructs and outputs a distance plot. Using various arguments a user can control the sliding window size, the window shift, a method of distance calculation, a sequence of interest, and the genome region, which will be included into distance calculation. Two methods of genetic distance calculation are implemented by now: pairwise and Kimura 2-parameter distance. After creating the plot, a user can also return distance data in the form of Pandas DataFrame object (to manipulate them, construct a plot using another library, or save the data), and data concerning the input alignment: sequence names and alignment length. The output plots are constructed using Plotly. It allows to interactively explore the distance data, and easily control the graphical output.
To test the package we used the nucleotide sequences which contain recombination events detected using the same method implemented in other tools. We’ve found that recan works fast both using short virus genome alignments, as hepatitis B virus, and long ones, as the lumpy skin disease virus genome which is about 150 000 bp in length. Recan produces the same output when compared to available tools (RAT, Simplot, RD4) which have the implementation of the distance methods. 
Recan can be installed by pip Python package manager using ‘pip install recan’ command. The source code and guide are available at https://github.com/babinyurii/recan.


# Citations

# Figures


# Acknowledgements



# References
