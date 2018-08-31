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
