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
