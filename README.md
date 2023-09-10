# gcBB: Space-Efficient Genome Comparison using BOSS representation and the BWSD similarity measure 
This software is an implementation of the gcBB algorithm described in Genome Comparison on Succinct Colored de Bruijn Graphs by Lucas P. Ramos, Felipe A. Louza and Guilherme P. Telles, String Processing and Information Retrieval: 29th International Symposium (2022).

Given a collection of _N_ genomes, gcBB uses the BOSS representation and BWSD computation to compare all genomes in the collection, outputting a distance matrix and a newick file that can be used to visualize the collection phylogeny. It's divided in three phases:
* Phase 1: Uses eGap algorithm to construct needed arrays (BWT,LCP,DA,CL);
* Phase 2: Constructs the BOSS representation for each pair of collection or the entire collection at once, depending on the option compiled;
* Phase 3: Computes the BWSD between all pairs of genomes.

## Install
```sh
git clone https://github.com/lucaspr98/gcBB --recursive
cd gcBB
```

## Pre-requisites
* A relatively recent version of *gcc*
* Python 3.X

### eGap
To construct the BOSS representation we will use the method described in [[1](https://doi.org/10.1186/s13015-019-0140-0)]. 
The eGap repository comes within clone with flag `--recursive`.

## Compile
To compile use the command:
```sh
make all
```
Options:
* `COVERAGE=1` to apply coverage to weight the comparison on BWSD. The default value is COVERAGE=0.
* `ALL_VS_ALL=1` to make eGap compute and merge all genomes instead of pairwise, constructs BOSS for one merge and computes BWSD using bitvector approach described in [[2](https://doi.org/10.1016/j.tcs.2019.03.012)]. The default value is ALL_VS_ALL=0.
* `DEBUG=1` to print information over the BOSS construction and BWSD computation. The default value is DEBUG=0.


Example:
```sh
make all ALL_VS_ALL=1 COVERAGE=1 DEBUG=1
```
**Obs**: use `make clean` command before `make all` with new options. 
## Run
The code of gcBB provides the possibility of comparing a pair of genomes or all pairs of genomes in a collection. After running the algorithm a directory named `results/` will be created containing:
* Two files containing the BWSD matrixes with the expectation and shannon's entropy between all pair of genomes;
* Two files containing the newick files using expectation and shannon's entropy between all pair of genomes to reconstruct the phylogeny;
* For each pair of genome, a file containing the BOSS and BWSD information. That is, _8*((N-1)*N/2)*_ files, where **N** is the number of genomes. (**ALL_VS_ALL=0**);
* One file containing the BOSS and BWSD information for the entire collection (**ALL_VS_ALL=1**);

### Pair of genomes comparison
To construct the BOSS representation and compute the BWSD between a pair of genomes run gcBB using the command:
```sh
./gcBB <path_to_dir> <input1.fastq> <input2.fastq>
```
**Example:**
```sh
./gcBB dataset/ -k 3 reads1.fastq reads2.fastq
```
In directory results, there will be the following files: 
* `reads1-reads2_expectation_k_16.dmat` and  `reads1-reads2_entropy_k_16.dmat`;
* `reads1-reads2_k_16.info`.

If the COVERAGE flag was used in make, `coverage` will be in the suffix of the outputted files.
### Genome collection comparison
To construct the BOSS representation and compute the BWSD between all pair of genomes from a directory run gcBB using the command:
```sh
./gcBB <path_to_dir>
```
**Example:**
```sh
./gcBB dataset/ -k 3 -p
```
Consider that `dataset/` contains the following genomes `reads1.fastq`, `reads2.fastq`, `reads3.fastq`.\
In directory results, there will be the following files: 
* `dataset_expectation_k_3.dmat` and  `dataset_entropy_k_3.dmat`;
* `reads1-reads2_k_3.info`, `reads1-reads3_k_3.info`, `reads2-reads3_k_3.info` (**ALL_VS_ALL=0**);
* `dataset_k_3_all.info` (**ALL_VS_ALL=1**).

Note that if the `ALL_VS_ALL` flag was used in make, `all` will be in the suffix of the outputted files.
### Command line options
*-k*, specify the size of k-mers used in the BOSS construction. The default value is k=32.

*-m*, specify the maximum usage of ram in MB provided to eGap and gcBB. The default value is m=2048.

*-p*, used to print BOSS files (last, w, wm, colors, coverage, summarized\_LCP, summarized\_SL) in results directory.

## References
[1] [*External memory BWT and LCP computation for sequence collections with applications*](https://doi.org/10.1186/s13015-019-0140-0);\
[2] [*Algorithms to compute the Burrows-Wheeler Similarity Distribution*](https://doi.org/10.1016/j.tcs.2019.03.012);

## Funding

Supported by the Brazilian agencies CAPES, CNPq (grant number 406418/2021-7) and FAPEMIG (grant number APQ-01217-22).
