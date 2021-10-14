#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 12:31:42 2020

@author: Ravin Poudel
"""


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

tempdir = "test/test_data/"
fq="test/test_data/sample.fq.gz"
filelist = [fq]


def test_is_zip(fq):
	keio.is_gzip(fq)

def test_fq2fa():
	keio.fq2fa(filelist)

def test_run_vsearch():
	keio.run_vsearch("test/test_data/upstream.fasta", "test/test_data/forward.fasta", cluster_id=0.75, minseq_length=5, tempdir=tempdir, threads=2)

# is_gzip(fq)

# fq2fa(filelist, tempdir)

# fastpath = os.path.join(tempdir, "forward.fasta")

# with open(fastpath, "w") as f1:
# 	for file in filelist:
# 		if is_gzip(file):
# 			with gzip.open(file, 'rt') as f:
# 				records = SeqIO.parse(f, "fastq")
# 				SeqIO.write(records, f1, "fasta")
# 		else:
# 			with open(file, 'r') as f:
# 				records = SeqIO.parse(f, "fastq")
# 				SeqIO.write(records, f1, "fasta")
# test_run_vsearch()
# for file in filelist:
# 	print(file)