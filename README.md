# RECAN
RECAN is a Python library to test DNA sequences for recombination events using distance plots. It generates the same distance plots as the RAT[1] and Simplot[2]. It is intended to be used in the Jupyter notebook.

## Requirements
To use RECAN, you will need:
- Python 3.6
- Biopython
- plotly 
- pandas
- Jupyter notebook

## Intallation
Download the repository and unzip it. To install the package  into your Python environment, run from the folder 'RECAN*' :

```
$ pip install .
```

## Usage example

Import 'simgen' function from the recan package:
```python
from recan.simgen import simgen
```

create an object of the Simgen class. To initialize the object pass your alignment in 'fasta' format as an argument:
```python
sim = Simgen("./data/hbv_C_Bj_Ba.fasta")
```
The input data are taken from the article by Sugauchi et al.(2002). This paper describes recombination event observed in hepatitis B virus isolates.

The object of the Simgen class has method 'get_info()' which shows information about the alignment. 
```python
sim.get_info()
```
```
index:	sequence id:
0	AB048704.1_genotype_C_
1	AB033555.1_Ba
2	AB010291.1_Bj
alignment length:  3215
```


We have three sequences in our alignment. 'Simgen' class is based upon the 'MultipleSequenceAlignment' class of the Biopython library.  So, we treat our alignment as the array with n_samples and n_features, where 'samples' are sequences themselves, and the features are columns of nucleotides in the alignment. Index corresponds to the sequence. Note, that indices start with 0.


After you've created the object you can draw the similarity plot. 
Call the method 'simgen' of the Simgen object to draw the plot. Pass the following parameters to the method:
- 'window': sliding window size. The number of nucleotides the sliding window will span. It has the value of 500 by default.
- 'shift': this is the step our window slides downstream the alignment. It's value is set to 250 by default
- 'pot_rec': the index of the potential recombinant: use method 'get_info()' to get the indices, especially if your alignment has many sequences. All the other sequences will be plotted as function of distance to that sequence.

The isolate of Ba genotype is the recombinant between the virus of C genotype and genotype Bj. Let's plot it. We set genotype Ba as the potential recombinant : 

```python
sim.simgen(window=200, shift=50, pot_rec=1)
```
![image](https://user-images.githubusercontent.com/28758465/44982136-7eec7a80-af7d-11e8-9bb7-e76935821eaf.png)

1. Recombination Analysis Tool (RAT): a program for the high-throughput detection of recombination. Bioinformatics, Volume 21, Issue 3, 1 February 2005, Pages 278–281, https://doi.org/10.1093/bioinformatics/bth500
2. https://sray.med.som.jhmi.edu/SCRoftware/simplot/ 
3.  Hepatitis B Virus of Genotype B with or without Recombination with Genotype C over the Precore Region plus the Core Gene. Fuminaka Sugauchi et al. JOURNAL OF VIROLOGY, June 2002, p. 5985–5992. 10.1128/JVI.76.12.5985-5992.2002 https://jvi.asm.org/content/76/12/5985 
