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
To compile use the command:
```sh
make all
```
Options:
* `COVERAGE=1` to apply coverage to weight the comparison on BWSD. The default value is COVERAGE=0.
* `BOSS_ALL=1` to make eGap compute and merge all genomes instead of pairwise, constructs only one BOSS struct for all genomes in the collection. The default value is BOSS_ALL=0.

Example:
```sh
make all BOSS_ALL=1 COVERAGE=1
```
**Obs**: use `make clean` command before `make all` with new options. 
## Run
The code of gcBB provides the possibility of comparing a pair of genomes or all pairs of genomes in a collection. After running the algorithm a directory named `results/` will be created containing:
* For each pair of genome, a file containing the BOSS representation constructed based on the merge of the pair;
* Two files containing the BWSD matrixes with the expectation and shannon's entropy between all pair of genomes;
* Two files containing the newick files using expectation and shannon's entropy between all pair of genomes to reconstruct the phylogeny;

### Pair of genomes comparison
To compute the BOSS representation and the BWSD between a pair of genomes run gcBB using the command:
```sh
./gcBB <path_to_dir> <input1.fastq> <input2.fastq>
```
**Example:**
```sh
./gcBB dataset/ -k 16 reads1.fastq reads2.fastq
```
In directory results, there will be following files: 

-`reads1-reads2_k_16.info`;

-`reads1-reads2_expectation_k_16.dmat` and  `reads1-reads2_entropy_k_16.dmat`.

If the coverage flag was used in make, then we will have the following files:

-`reads1-reads2_expectation_k_16_coverage.dmat` and `reads1-reads2_entropy_k_16_coverage.dmat`. 

### Genome collection comparison
To compute the BOSS representation and the BWSD between all pair of genomes from a directory run gcBB using the command:
```sh
./gcBB <path_to_dir>
```
**Example:**
```sh
./gcBB influenza_dataset/ -k 16 -p
```
In directory results, there will be _8*((N-1)*N/2)*_ files containing all possible pair of genomes in the directory BOSS representations, where **N** is the number of genomes in the directory, and `influenza_dataset_expectation_k_16.dmat` and  `influenza_dataset_entropy_k_16.dmat`.

### Command line options
*-k*    
    specify the size of k-mers used in the BOSS construction. The default value is k=32.

*-m*    
    specify the maximum usage of ram in MB provided to eGap and gcBB. The default value is m=2048.

*-p*    
    used to print BOSS files (last, w, wm, colors, coverage, summarized\_LCP, summarized\_SL) in results directory.

## Funding

Supported by the Brazilian agencies CAPES, CNPq (grant number 406418/2021-7) and FAPEMIG (grant number APQ-01217-22).
