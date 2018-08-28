# RECAN
RECAN is a Python library to test DNA sequences for recombination events using distance plots. It generates the same distance plots as the RAT[1] and Simplot[2]. It is intended to be used in the Jupyter notebook.

## Requirements
To use RECAN, you will need:
- Python 3.6
- Biopython 1.7
- plotly 3.1
- pandas
- Jupyter notebook

## Intallation
Download the repository, unzip it and run from the folder 'RECAN*' to install it into your Python environment:

```
$ pip install .
```

## Usage example

Import 'simgen' function from the recan package:
```
from recan.simgen import simgen
```




































1. Recombination Analysis Tool (RAT): a program for the high-throughput detection of recombination. Bioinformatics, Volume 21, Issue 3, 1 February 2005, Pages 278–281, https://doi.org/10.1093/bioinformatics/bth500
2. https://sray.med.som.jhmi.edu/SCRoftware/simplot/ 
