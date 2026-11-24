# Metagenomic-Community-Complexity-Estimator

This project analyzes real metagenomic sequencing data to estimate the complexity of a microbial community purely from raw DNA sequences, without relying on reference databases. By measuring the Shannon entropy of short DNA words (k-mers), the script infers the diversity of organisms present in the sample. A higher k-mer entropy indicates a more diverse and complex microbial community.

# Project Overview
Instead of isolating individual species, metagenomic studies sequence all DNA fragments in an environmental sample together, producing millions of short reads. This script measures community diversity by analyzing k-mer usage patterns across these reads.   
Key features of this implementation include:
+ Efficient Memory Usage: Stream-processes large, compressed FASTQ files (```.fq.gz```) line-by-line using generators, without loading the entire file into memory.
+ Data Parsing: Extracts only the DNA sequence lines from the FASTQ format (skipping headers and quality scores) and ignores k-mers containing ambiguous bases ('N').
+ Dynamic Convergence: Computes the Shannon entropy of the k-mer distribution every 1,000 reads. It automatically stops processing when the entropy converges (i.e., the difference between current and previous entropy is strictly less than 0.00001).
+ Visualization: Generates a richness plot showing the number of unique k-mers discovered over the number of reads processed, visually demonstrating community richness and redundancy.

# Setup and Data 
To run the analysis, you need a metagenomic FASTQ file:   
1. Go to the [NCBI Trace Archive](https://www.ncbi.nlm.nih.gov/Traces/index.html?view=study&acc=ERP000108).
2. Download a file ending in ```_clean.1.fq.gz``` (for example, ```MH0026_081224.clean.1.fq.gz```).
3. Do not decompress the file. The Python script uses the ```gzip``` module to read the compressed ```.fq.gz``` archive directly.

# Usage
The main script is ```Project_4.py```. It accepts a required FASTQ file and an optional k-mer size (an integer between 4 and 7).   
## 1. Single-k Mode
If you specify a k-mer size, the script will run the analysis for that specific length and output a single plot.   
```
python Project_4.py MH0026_081224.clean.1.fq.gz 6
```
## 2. Multi-k Mode
If you omit the k-mer size, the script automatically defaults to "Multi-k mode". It will run the convergence algorithm sequentially for k = 4, 5, 6, and 7.   
```
python Project_4.py MH0026_081224.clean.1.fq.gz
```
