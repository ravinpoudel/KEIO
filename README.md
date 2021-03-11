#KEIO
## KEIO: A python software to process illumina reads for keio-collection type project.


## Authors

* Ravin Poudel, PhD, Department of Microbiology and Cell Science, University of Florida
* Adam R. Rivers, PhD , US Department of Agriculture, Agricultural Research Service
* Christopher Reisch, PhD, Department of Microbiology and Cell Science, University of Florida


## Introduction



## Installation

KEIO can be installed from:

1. The Github repository: https://github.com/ravinpoudel/KEIO.git

```{bash}

    git clone https://github.com/ravinpoudel/KEIO.git
    
```


## Dependencies

Following are the required softwares/programs.

* ``pybedtools``
* ``NMSLib``
* ``Biopython``
* ``Pandas``


## Usage

```
KEIO: A python software to process illumina reads for keio-collection type project.

optional arguments:
  -h, --help            show this help message and exit
  --fastq FASTQ [FASTQ ...], -f FASTQ [FASTQ ...]
                        input fastq file
  --upstreamFasta UPSTREAMFASTA, -uf UPSTREAMFASTA
                        A upstreamFasta file
  --downstreamrcFasta DOWNSTREAMRCFASTA, -drcf DOWNSTREAMRCFASTA
                        A downstreamFasta file
  --threads THREADS     The number of cpu threads to use
  --tempdir TEMPDIR     The temp file directory
  --keeptemp            Should intermediate files be kept?

```

## Examples


Use case: Retrieving target sequence for a given PAM motif in the forward and reverse strands, where length of guide sequence is 20 base pair.
12 base pair close to PAM motif is conserved ,i.e. unique and the full sequence has a hamming distance of more than 2.
Here the number of used threads is 2
Return a table of pam sites and associated data, at in current folder.

```
		keio -f test/test_data/sample.fq.gz -uf test/test_data/upstream.fasta -drcf test/test_data/downstream_rc.fasta --threads 8

```

## API documentation


## License information

