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
applications*](https://doi.org/10.1186/s13015-019-0140-0). To compile eGap run the following commands:
```sh
cd egap
make
cd ../
```
The eGap repository comes within clone with flag `--recursive`

## Compile
```sh
gcc main.c -o gcBB
```

## Run
The code of gcBB provides the possibility of comparing a pair of genomes or all pairs of genomes in a collection. After running the algorithm a directory named `results/` will be created containing:
* For each pair of genome, a file containing the BOSS representation constructed based on the merge the pair;
* A file containing the BWSD matrixes with the expectation and shannon's entropy between all pair of genomes;

### Pair of genomes comparison
To compute the BOSS representation and the BWSD between a pair of genomes run gcBB using the command:
```sh
./gcBB <path_to_dir> <input1.fastq> <input2.fastq> <k>
```
**Example:**
```sh
./gcBB dataset/ reads1.fastq reads2.fastq 3
```
In directory results, there will be two files, `1-2.boss` and `bwsd_matrixes.txt`.

### Genome collection comparison
To compute the BOSS representation and the BWSD between all pair of genomes from a directory run gcBB using the command:
```sh
./gcBB <path_to_dir> <k>
```
**Example:**
```sh
./gcBB influenza_dataset/ 3
```
In directory results, there will be _((N-1)*N/2)_ files containing all possible pair of genomes in the directory BOSS representation, where **N** is the number of genomes in the directory, and `bwsd_matrixes.txt`.