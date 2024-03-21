![RabbitTrim logo](RabbitTrim_logo.png "RabbitTrim logo")

# `RabbitTrim version 1.0.0`
RabbitTrim is an efficient adapter removal tool based on Trimmomatic, which implemented in C++ and deeply optimized in performance. 
RabbitTrim refactors and improves the efficiency of Trimmomatic in processing plain FASTQ data and gzip-compressed data.
RabbitTrim supports all function modules available in Trimmomatic and maintains identical results.  

## Installation
`RabbitTrim` is written in `C++` for GNU Linux/Unix platforms.`

### Dependancy
* cmake version 3.20 or later
* gcc version 9.4.0 or later
* [intel isa-l](https://github.com/intel/isa-l)
* [zlib](http://zlib.net/)

### Compile and install 
```bash
git clone git@github.com:RabbitBio/RabbitTrim.git 
cd RabbitTrim
mkdir build && cd build
cmake .. -DIGZIP_PATH=[your installation path for intel isa-l]
make -j
```
## Usage
```bash
Usage: ./RabbitTrim trimmomatic [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -P,--PE,--no-PE{false} Needs: --forward --reverse Excludes: --SE
                              specify whether to be pair end data
  -S,--SE,--no-SE{false} Needs: --forward Excludes: --PE --reverse
                              specify whether to be single end data
  -t,--threads INT:POSITIVE   specify the number of threads
  -p,--phred INT:{0,33,64}    specify the baseline of the phred score
  -f,--forward TEXT:PATH(existing) ...
                              specify the path to the forward read file
  -r,--reverse TEXT:PATH(existing) ... Excludes: --SE
                              specify the path to the reverse read file
  -o,--output TEXT REQUIRED   specify the path to the output file
  --log TEXT                  specify the path to the trim log file
  --stats TEXT REQUIRED       specify the path to the trim statistical data file
  -s,--steps TEXT ... REQUIRED
                              specify the steps to be performed
  -c,--compressLevel INT:compression level is limited to be between 0 and 9
                              specify the compression level for the output file
  --quiet,--no-quiet{false}   specify whether to print program runtime information
  --validatePair,--no-validatePair{false}
                              specify whether to validate pair data
  --pigz,--no-pigz{false}     specify whether to use pigz
  -g,--pigzThreadsNum INT:POSITIVE:INT bounded to [2 - 64] Needs: --pigz
                              specify the max thread number of pigz
```
## Example
1. Process plain FASTQ data in simple mode (for Single-End data)
```bash
./RabbitTrim trimmomatic --SE -f /data/demo.read1.fastq -o /data/output/out.fastq --stats /data/output/trim_stat -s ILLUMINACLIP:../adapter/TruSeq-DNA-Free-SE.fa:2:30:12:1:true MINLEN:36 -t 1
```
2. Process plain FASTQ data in palindrome mode (for Paired-End data)
```bash
./RabbitTrim trimmomatic --PE -f /data/demo.read1.fastq -r /data/demo.read2.fastq -o /data/output/out.fastq --stats /data/output/trim_stat -s ILLUMINACLIP:../adapter/TruSeq-DNA-Free-PE.fa:2:30:12:1:true MINLEN:36 -t 1
```
```
3. Process gzip-compressed data in simple mode (for Single-End data)
```bash
./RabbitTrim trimmomatic --SE -f /data/demo.read1.fastq.gz -o /data/output/out.fastq.gz --stats /data/output/trim_stat -s ILLUMINACLIP:../adapter/TruSeq-DNA-Free-SE.fa:2:30:12:1:true MINLEN:36 -t 1
```
4. Process gzip-compressed data in palindrome mode (for Paired-End data)
```bash
./RabbitTrim trimmomatic --PE -f /data/demo.read1.fastq.gz -r /data/demo.read2.fastq.gz -o /data/output/out.fastq.gz --stats /data/output/trim_stat -s ILLUMINACLIP:../adapter/TruSeq-DNA-Free-PE.fa:2:30:12:1:true MINLEN:36 -t 1
```
## Output
The output of RabbitTrim is identical to Trimmomatic.
For single-end data, it typically includes one processed FASTQ file (either plain or gzip-compressed format) and one trimming information statistics file.
For paired-end data, it typically includes four FASTQ files and one trimming information statistics file.

# Bug Report
We highly appreciate all bug reports, comments, and suggestions from our users.  
Please feel free to raise any concerns or feedback with us without hesitation by `issue`. 

## Cite
