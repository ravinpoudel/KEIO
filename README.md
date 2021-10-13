## KEIO: A python software to process illumina reads for keio-collection type project.

The Keio collection in Escherichia coli K-12 represents a collection of single-gene deleted mutants. This collection was created manually one by one by replacing the predicted ORF with a kanamycin cassette to inactivate chromosomal genes. Then primers were designed to create in-frame deletions upon excision of the resistance cassette. Of 4288 genes targeted, 3985 mutants were obtained. Majority of these mutants represents mutation of non-essential genes. [More @](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1681482/pdf/msb4100050.pdf)

Primarily, Keio collection provides a molecular tool/resource to understand functional and physiological aspects of gene at the system levels. Although, creating such molecular collection / tools takes lot of resources and time. Thus, here we explored TnSeq methods to create a single gene mutant type collection at high-throughput scale. Once the Randomly Barcode Transposons are created, the constructs are randomly inserted into bacterial genome to create mutated clone. Each clone theoretically should represent a single gene mutation. [More @](https://mbio.asm.org/content/6/3/e00306-15)

Our methods involves mapping the location of random barcode sequence (about 20 base pair in lengths) as well as inline barcodes. Mapping-positions are then mapped to the genome coordinates to identity targeted features, whereas in-line barcodes mapping allows to map the clone to specific plates ( we use 384 wells plate).


Barcode Insertion Mapping             |  Plate Mapping
:-------------------------:|:-------------------------:
![](https://raw.githubusercontent.com/ravinpoudel/KEIO/master/keio/data/instruction.png)  |  ![](https://raw.githubusercontent.com/ravinpoudel/KEIO/master/keio/data/platemapping.png)

## Authors

* Ravin Poudel, PhD, Department of Microbiology and Cell Science, University of Florida
* Adam R. Rivers, PhD , US Department of Agriculture, Agricultural Research Service
* Christopher Reisch, PhD, Department of Microbiology and Cell Science, University of Florida


## Introduction



## Installation

KEIO can be installed from GitHub:

```bash
# Create a conda environment and install and pybedtools
conda create -n keio python=3.7 pybedtools=0.8.2
conda activate keio

git clone https://github.com/ravinpoudel/KEIO.git
cd keio
pip install .

# check if the installation works
keio -h
    
```


## Dependencies

Following are the required softwares/programs.


* ``VSEARCH``
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


```{bash}
keio -f test/test_data/sample.fq.gz \
     -uf test/test_data/upstream.fasta \
     -drcf test/test_data/downstream_rc.fasta \
     --threads 8

```




[![IMAGE ALT TEXT HERE](https://raw.githubusercontent.com/ravinpoudel/KEIO/master/keio/data/applogo.png)](https://ravinpoudel.shinyapps.io/keioplatemapper/)


[![IMAGE ALT TEXT HERE](https://raw.githubusercontent.com/ravinpoudel/KEIO/master/keio/data/KeioMapper.png)](https://ravinpoudel.shinyapps.io/keioplatemapper/)




## API documentation

[Available Here](https://ravinpoudel.github.io/KEIO/)


## License information

