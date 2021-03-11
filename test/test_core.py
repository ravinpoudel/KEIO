import os
import sys
import logging
import argparse
import tempfile
import shutil
import textdistance
import hashlib
import subprocess
import Bio
import pickle
import nmslib
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
import gzip
import keio

#######

tempdir='/var/folders/ct/2nl65wtn3978tk0rq17mg2fw0000gn/T/keio_vpjbfm4f'
mapping_fasta = "test/test_data/upstream.fasta"
reads_fasta ='/var/folders/ct/2nl65wtn3978tk0rq17mg2fw0000gn/T/keio_vpjbfm4f/forward.fasta'
cluster_id=0.95
minseq_length=8
threads = 2

keio.run_vsearch(mapping_fasta, reads_fasta, tempdir=tempdir)



with open(file, 'r') as f:
    for line in f:
        print(line)


with open(file, 'r') as f:
    sdict = {}
    for line in f:
        ll = line.strip().split()
        qname = ll[0]
        tname = ll[1]
        pmatch = float(ll[4])
        tcov = float(ll[5])
        gaps = int(ll[8])
        spos = int(ll[12])
        epos = int(ll[13])
        qstrand = ll[14]
        tstrand = ll[15]
        td = {'spos': spos, 'epos': epos, 'pmatch': pmatch, 'tcov': tcov,
                'gaps': gaps, 'qstrand': qstrand, 'tstrand': tstrand}
        if tcov >= 50:
            if gaps < 5:
                if qname in sdict:
                    if tname in sdict[qname]:
                        if sdict[qname][tname]['pmatch'] < pmatch:
                            sdict[qname][tname] = td
                    else:
                        sdict[qname][tname] = td
                else:
                    sdict[qname] = {}
                    sdict[qname][tname] = td







tempdir = '/var/folders/ct/2nl65wtn3978tk0rq17mg2fw0000gn/T/keio_byqni62i'
fastq = "test/test_data/sample.fq.gz"

def test_fq2fa():
    keio.fq2fa(fastq, tempdir)
    

keio.fq2fa(fastq, tempdir)



tempdir = '/var/folders/ct/2nl65wtn3978tk0rq17mg2fw0000gn/T/keio_u0dumd2j'
reads_fasta ='/var/folders/ct/2nl65wtn3978tk0rq17mg2fw0000gn/T/keio_u0dumd2j/forward.fasta'
mapping_fasta = "test/test_data/upstream.fasta"
cluster_id=0.95
minseq_length=8
threads = 2
run_vsearch(args.upstreamFasta, fastapath, cluster_id=0.95, minseq_length=8,  threads=2)