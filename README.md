# gcBB - Space-Efficient Genome Comparison using BOSS representation and the BWSD similarity measure 


## Install
```sh
git clone https://github.com/lucaspr98/gcBB --recursive
cd gcBB
```

## Pre-requisites
* A relatively recent version of *gcc*
* Python 3.X

### eGap
To construct the BOSS representation we will use the method described in [*External memory BWT and LCP computation for sequence collections with 
applications*](https://doi.org/10.1186/s13015-019-0140-0). 
The eGap repository comes within clone with flag `--recursive`

## Compile
If you do not want the BWSD to use coverage to weight the comparison, use the command:
```sh
make
```
or 
```sh
make COVERAGE=0
```
else, use the command:
```sh
make COVERAGE=1
```

## Run
The code of gcBB provides the possibility of comparing a pair of genomes or all pairs of genomes in a collection. After running the algorithm a directory named `results/` will be created containing:
* For each pair of genome, a file containing the BOSS representation constructed based on the merge the pair;
* A file containing the BWSD matrixes with the expectation and shannon's entropy between all pair of genomes;

### Pair of genomes comparison
To compute the BOSS representation and the BWSD between a pair of genomes run gcBB using the command:
```sh
./gcBB <path_to_dir> <input1.fastq> <input2.fastq>
```
**Example:**
```sh
./gcBB dataset/ reads1.fastq reads2.fastq 
```
In directory results, there will be eight files: 
-`reads1-reads2.boss-info`
-`reads1-reads2.2.last`
-`reads1-reads2.2.W`
-`reads1-reads2.2.Wm`
-`reads1-reads2.2.colors`
-`reads1-reads2.4.coverage`
-`reads1-reads2.2.reduced_lcp`
-`reads1-reads2_distance_matrixes_coverage_0.txt` or `reads1-reads2_distance_matrixes_coverage_1.txt`(depending on the compilation flag). 

### Genome collection comparison
To compute the BOSS representation and the BWSD between all pair of genomes from a directory run gcBB using the command:
```sh
./gcBB <path_to_dir>
```
**Example:**
```sh
./gcBB influenza_dataset/
```
In directory results, there will be _8*((N-1)*N/2)*_ files containing all possible pair of genomes in the directory BOSS representations, where **N** is the number of genomes in the directory, and `influenza_dataset_distance_matrixes_coverage_0.txt` or `influenza_dataset_distance_matrixes_coverage_1.txt`(depending on the compilation flag).

### Command line options
*-k*    
    specify the size of k-mers used in the BOSS construction. The default value is k=30.

*-m*    
    specify the size of the blocks read from the files constructed by the eGap. The default value is m=10000.

