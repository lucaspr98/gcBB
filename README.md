# gcBB - Space-Efficient Genome Comparison using BOSS representation and the BWSD similarity measure 

## Install
```sh
git clone https://github.com/lucaspr98/gcBB --recursive
cd gcBB
```

## Pre-requisites
```sh
cd egap
make
cd ../
```

## Compile
```sh
gcc main.c -o gcBB
```

## Run
```sh
./gcBB <path/input1.fastq> <path/input2.fastq> <k>
```
Example:
```sh
./gcBB dataset/reads1.fastq dataset/reads2.fastq 3
```
