![RabbitTrim logo](RabbitTrim_logo.png "RabbitTrim logo")

# `RabbitTrim version 2.0.0`
RabbitTrim is an efficient adapter removal tool based on Trimmomatic and Ktrim, which implemented in C++ and deeply optimized in performance. 
RabbitTrim refactors and improves the efficiency of Trimmomatic and Ktrim in processing plain FASTQ data and gzip-compressed data.
RabbitTrim supports all function modules available in Trimmomatic and Ktrim and maintains identical results.  

## Installation
`RabbitTrim` is written in `C++` for GNU Linux/Unix platforms.

### Dependancy
* cmake version 3.20 or later [REQUIRED]
* gcc version 9.4.0 or later [REQUIRED]
* [zlib](http://zlib.net/) [REQUIRED]
* [intel isa-l](https://github.com/intel/isa-l) [OPTIONAL]

### Compile and install 
```bash
git clone https://github.com/RabbitBio/RabbitTrim.git
cd RabbitTrim
mkdir build && cd build
cmake .. -DIGZIP_PATH=[your installation path for intel isa-l]
make -j
```
## Trimmomatic Usage
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
                              specify the max thread number of pigz
```

## Ktrim Usage
```bash
Usage: ./RabbitTrim ktrim [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -P,--PE,--no-PE{false} Needs: --forward --reverse Excludes: --SE
                              specify whether to be pair end data
  -S,--SE,--no-SE{false} Needs: --forward Excludes: --PE --reverse
                              specify whether to be single end data
  -f,--forward TEXT:PATH(existing) ...
                              specify the path to the forward read file
  -r,--reverse TEXT:PATH(existing) ... Excludes: --SE
                              specify the path to the reverse read file
  -o,--output TEXT REQUIRED   specify the path to the output file
  -t,--threads INT:POSITIVE   specify the number of threads
  -p,--phred INT:{0,33,64}    specify the baseline of the phred score
  -m,--mismatch FLOAT:FLOAT in [0 - 1]
                              Set the proportion of mismatches allowed during index and sequence comparison, Default: 0.125
  -l,--minlen INT:POSITIVE    specify the minimum read size to be kept
  -q,--quality INT:POSITIVE   specify the minimum quality score
  -w,--window INT:POSITIVE    specify the window size for quality check
  -a,--seqA TEXT              specify the adapter sequence of read1
  -b,--seqB TEXT              specify the adapter sequence of read2
  -k,--seqKit TEXT:{BGI,Illumina,Nextera,Transposase}
                              specify the sequencing kit to use built-in adapters
  -c,--compressLevel INT:compression level is limited to be between 0 and 9
                              specify the compression level for the output file
```

## Output
The output of RabbitTrim is identical to Trimmomatic and Ktrim.


# Bug Report
We highly appreciate all bug reports, comments, and suggestions from our users.  
Please feel free to raise any concerns or feedback with us without hesitation by `issue`. 

